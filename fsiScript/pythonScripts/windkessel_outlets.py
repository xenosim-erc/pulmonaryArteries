#!/usr/bin/env python3
"""Split a total three-element Windkessel across the outlets by outlet area.

The whole distal bed is described by one proximal resistance, one distal
resistance and one compliance. Each outlet is given the fraction

    f_i = A_i^(k/2) / sum_j A_j^(k/2)

of the total flow, where A_i is the outlet area and k the exponent of the flow
split Q_i ~ r_i^k: k = 2 is the area rule (Q proportional to area), k = 3 is
Murray's law. Outlet i then gets

    Rch_i = Rp / f_i,    R_i = Rd / f_i,    C_i = C * f_i

The two resistances are inversely proportional to the flow fraction and the
compliance directly proportional to it, so a larger outlet gets less resistance
and more compliance and draws more of the flow. The outlets in parallel
reproduce the totals exactly, and every outlet keeps the time constant
R_i C_i = Rd C of the whole bed. Only area ratios are used, so the areas may be
in any unit.

As a check, the values as written are recombined in parallel (resistances
1 / sum(1/R_i), compliance sum(C_i)) and compared with the totals; a mismatch
beyond rounding stops the script.

With --settled-start WAVEFORM the outlets start at the pressure the
Windkessel settles to at the start of a beat when driven by that inflow, not
at --initial-pressure: the compliance charge is the slowest part of the
model to settle (time constant Rd C, so a starting error shrinks only by
exp(-T / (Rd C)) per beat), and starting it settled leaves little to settle.
When the 3D wall's compliance has been taken off C (--wall-compliance, with
the pressure the wall starts at, --wall-pressure), the model is two
compliances: the wall's at the inlet, then Rp, then the outlets' and Rd. The
wall then starts at its own pressure (the imaging pressure, where the
prestressed wall is in equilibrium), and the outlets start where the total
charge on both compliances equals that of the settled state; the small
imbalance between the two evens out within a fraction of a beat (time
constant about Rp Cwall Cout / (Cwall + Cout)).

The result is written as a dictionary fragment of per-outlet patch entries,
which the template's 0/fluid/p includes. That include is the only definition
of the outlet conditions; the old uniform "outlet.*" entry is commented out.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path


# Values are written to 6 significant figures, so recombining them can drift
# from the totals by about 1e-6; anything larger means the split is wrong.
CHECK_TOLERANCE = 1e-4

RHO = re.compile(r"(?m)^[ \t]*rho\b[^;]*?(?P<value>[0-9.eE+-]+)[ \t]*;")


def positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise argparse.ArgumentTypeError("must be finite and greater than 0")
    return number


def read_outlet_areas(filename: Path) -> list[tuple[str, float]]:
    with filename.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = {"area", "name"} - set(reader.fieldnames or ())
        if missing:
            raise RuntimeError(
                f"profile CSV is missing columns: {', '.join(sorted(missing))}; "
                "regenerate it with the current pulmonary_centerlines.py"
            )
        outlets = []
        for row in reader:
            if not row["name"].startswith("outlet"):
                continue
            area = float(row["area"])
            if not math.isfinite(area) or area <= 0.0:
                raise RuntimeError(f"{row['name']} has a non-positive area: {area}")
            outlets.append((row["name"], area))
    if not outlets:
        raise RuntimeError(f"no outlet rows found in {filename}")
    # Natural order (outlet2 before outlet10) keeps the file easy to read.
    outlets.sort(key=lambda item: int(re.sub(r"\D", "", item[0]) or 0))
    return outlets


def read_density(filename: Path | None, default: float) -> float:
    if filename is None or not filename.is_file():
        return default
    match = RHO.search(filename.read_text(encoding="utf-8"))
    return float(match.group("value")) if match else default


def read_waveform(filename: Path) -> list[tuple[float, float]]:
    """(time, flow) rows of an OpenFOAM table: lines of the form ( t Q )."""
    rows = []
    for line in filename.read_text(encoding="utf-8").splitlines():
        parts = line.replace("(", " ").replace(")", " ").split()
        if len(parts) >= 2:
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    if len(rows) < 2:
        raise RuntimeError(f"no (time flow) rows in {filename}")
    return rows


def settled_start_pressure(rows, distal_resistance: float, compliance: float, steps: int = 20000) -> float:
    """Compliance pressure at the start of a beat once the Windkessel is
    periodic, for the total model C dp/dt = Q - p/Rd driven by the waveform
    (one period, the table's span). The equation is linear, so the pressure
    after one beat is a p0 + b with a = exp(-T/(Rd C)): integrate one beat
    from p0 = 0 to get b, then p* = b / (1 - a). Exponential integrator with
    the flow held at its midpoint value over each of the steps."""
    t0, t1 = rows[0][0], rows[-1][0]
    period = t1 - t0
    tau = distal_resistance * compliance
    dt = period / steps
    decay = math.exp(-dt / tau)

    def flow(t):
        for (ta, qa), (tb, qb) in zip(rows, rows[1:]):
            if ta <= t <= tb:
                return qa + (qb - qa) * (t - ta) / (tb - ta) if tb > ta else qa
        return rows[-1][1]

    p = 0.0
    for i in range(steps):
        q = flow(t0 + (i + 0.5) * dt)
        p = p * decay + q * distal_resistance * (1.0 - decay)
    return p / (1.0 - math.exp(-period / tau))


def settled_start_two_compliances(rows, rp: float, rd: float, c_out: float, c_wall: float,
                                  p_wall: float, steps: int = 20000) -> tuple[float, float, float]:
    """Outlet start pressure for the wall + outlet compliance model
        Cwall dp1/dt = Q - (p1 - p2)/Rp,   Cout dp2/dt = (p1 - p2)/Rp - p2/Rd,
    so that, with the wall starting at p_wall, the total charge
    Cwall p1 + Cout p2 equals that of the periodic state. Returns
    (outlet start, periodic p1, periodic p2) at the start of a beat. The
    equations are linear: one beat maps the state x to A x + b, found by
    integrating from zero and from the two unit states (RK4), and the
    periodic state solves (I - A) x = b."""
    import bisect
    times = [r[0] for r in rows]
    t0, period = rows[0][0], rows[-1][0] - rows[0][0]
    dt = period / steps

    def flow(t):
        i = min(max(bisect.bisect_right(times, t) - 1, 0), len(rows) - 2)
        (ta, qa), (tb, qb) = rows[i], rows[i + 1]
        return qa + (qb - qa) * (t - ta) / (tb - ta) if tb > ta else qa

    def beat(p1, p2, forced):
        def f(t, a, b):
            q = flow(t) if forced else 0.0
            qo = (a - b) / rp
            return (q - qo) / c_wall, (qo - b / rd) / c_out
        for i in range(steps):
            t = t0 + i * dt
            k1 = f(t, p1, p2)
            k2 = f(t + dt / 2, p1 + dt / 2 * k1[0], p2 + dt / 2 * k1[1])
            k3 = f(t + dt / 2, p1 + dt / 2 * k2[0], p2 + dt / 2 * k2[1])
            k4 = f(t + dt, p1 + dt * k3[0], p2 + dt * k3[1])
            p1 += dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            p2 += dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        return p1, p2

    b = beat(0.0, 0.0, True)
    a1 = beat(1.0, 0.0, False)
    a2 = beat(0.0, 1.0, False)
    # (I - A) x = b, A = [[a1[0], a2[0]], [a1[1], a2[1]]]
    m11, m12, m21, m22 = 1 - a1[0], -a2[0], -a1[1], 1 - a2[1]
    det = m11 * m22 - m12 * m21
    x1 = (b[0] * m22 - m12 * b[1]) / det
    x2 = (m11 * b[1] - m21 * b[0]) / det
    outlet = (c_wall * x1 + c_out * x2 - c_wall * p_wall) / c_out
    return max(outlet, 0.0), x1, x2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiles_csv", type=Path)
    parser.add_argument("output", type=Path, help="dictionary fragment to write")
    parser.add_argument(
        "--proximal-resistance", type=positive_float, default=3.3e6,
        help="total proximal (characteristic) resistance Rp [Pa s/m^3] (default: 3.3e6)",
    )
    parser.add_argument(
        "--distal-resistance", type=positive_float, default=2.0e7,
        help="total distal resistance Rd [Pa s/m^3] (default: 2e7)",
    )
    parser.add_argument(
        "--compliance", type=positive_float, default=3.5e-8,
        help="total compliance C [m^3/Pa] (default: 3.5e-8)",
    )
    parser.add_argument(
        "--exponent", type=positive_float, default=2.0,
        help="flow-split exponent k in Q ~ r^k: 2 area rule, 3 Murray (default: 2)",
    )
    parser.add_argument(
        "--transport-properties", type=Path,
        help="read rho from this transportProperties instead of --rho",
    )
    parser.add_argument(
        "--rho", type=positive_float, default=1060.0,
        help="fluid density used by windkesselPressure (default: 1060)",
    )
    parser.add_argument(
        "--initial-pressure", type=float, default=0.0,
        help=(
            "pressure every outlet starts at [Pa] (default: 0). windkesselPressure "
            "takes its starting pressure from 'value', which holds the kinematic "
            "pressure p/rho, so this starts the compliance charged to it"
        ),
    )
    parser.add_argument(
        "--settled-start", type=Path, metavar="WAVEFORM",
        help=(
            "start the outlets at the pressure the Windkessel settles to at the "
            "start of a beat when driven by this inflow table (constant/fluid/"
            "inletFlow format), instead of --initial-pressure"
        ),
    )
    parser.add_argument(
        "--wall-compliance", type=positive_float,
        help=(
            "with --settled-start: compliance of the 3D wall [m^3/Pa], already taken "
            "off --compliance; the settled start then accounts for it"
        ),
    )
    parser.add_argument(
        "--wall-pressure", type=float,
        help="with --wall-compliance: pressure the 3D wall starts at [Pa] (the imaging pressure)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    profiles_csv = args.profiles_csv.expanduser().resolve()
    if not profiles_csv.is_file():
        raise FileNotFoundError(profiles_csv)

    outlets = read_outlet_areas(profiles_csv)
    rho = read_density(args.transport_properties, args.rho)
    if args.initial_pressure < 0.0:
        raise ValueError("--initial-pressure must not be negative")
    start_note = "given"
    if args.settled_start is not None:
        rows = read_waveform(args.settled_start.expanduser())
        if abs(rows[0][1]) > 1e-3 * max(abs(q) for _, q in rows):
            print(f"Warning: the inflow at t = {rows[0][0]:g} is not zero, so the outlet pressure at the "
                  "start also has an Rp Q part, which this does not include", file=sys.stderr)
        if args.wall_compliance is not None:
            if args.wall_pressure is None:
                raise ValueError("--wall-compliance needs --wall-pressure")
            settled, x1, x2 = settled_start_two_compliances(
                rows, args.proximal_resistance, args.distal_resistance, args.compliance,
                args.wall_compliance, args.wall_pressure)
            print(f"Settled start (3D wall {args.wall_compliance:.4g} m^3/Pa + outlets "
                  f"{args.compliance:.4g} m^3/Pa): periodic state {x1 / 133.322387415:.3g} mmHg "
                  f"in the wall, {x2 / 133.322387415:.3g} mmHg at the outlets; the wall starts at "
                  f"{args.wall_pressure / 133.322387415:.3g} mmHg, so the outlets start at "
                  f"{settled:.6g} Pa ({settled / 133.322387415:.3g} mmHg) for the same total charge")
        else:
            settled = settled_start_pressure(rows, args.distal_resistance, args.compliance)
            print(f"Settled start: {settled:.6g} Pa ({settled / 133.322387415:.3g} mmHg), from "
                  f"{args.settled_start.name} (period {rows[-1][0] - rows[0][0]:g} s, "
                  f"Rd C = {args.distal_resistance * args.compliance:.3g} s), instead of "
                  f"{args.initial_pressure:g} Pa")
        args.initial_pressure = settled
        start_note = f"settled start from {args.settled_start.name}"
    initial_value = args.initial_pressure / rho
    weights = [area ** (args.exponent / 2.0) for _, area in outlets]
    total_weight = sum(weights)

    lines = [
        "// Generated by pythonScripts/windkessel_outlets.py; edit the totals",
        "// and rerun rather than editing these entries by hand.",
        f"//   Rp = {args.proximal_resistance:g} Pa s/m^3, "
        f"Rd = {args.distal_resistance:g} Pa s/m^3, C = {args.compliance:g} m^3/Pa",
        f"//   flow split Q ~ r^{args.exponent:g}, rho = {rho:g}",
        f"//   initial pressure {args.initial_pressure:.6g} Pa ({start_note}; "
        f"value = p/rho = {initial_value:.6g})",
        "//   Rch = Rp/f and R = Rd/f (inversely proportional to the flow",
        "//   fraction f); C = C*f (directly proportional to it).",
        "",
    ]
    print(f"{'Outlet':<10} {'flow fraction':>13} {'Rch':>12} {'R':>12} {'C':>12}")
    written = []
    for (name, _), weight in zip(outlets, weights):
        fraction = weight / total_weight
        # Round once here so the check below uses exactly what the file holds.
        rch = float(f"{args.proximal_resistance / fraction:.6g}")
        r = float(f"{args.distal_resistance / fraction:.6g}")
        c = float(f"{args.compliance * fraction:.6g}")
        written.append((fraction, rch, r, c))
        print(f"{name:<10} {fraction:>13.4f} {rch:>12.4g} {r:>12.4g} {c:>12.4g}")
        lines += [
            f"{name}",
            "{",
            "    type        windkesselPressure;",
            f"    // flow fraction {fraction:.6f}",
            f"    R           {r:.6g};",
            f"    Rch         {rch:.6g};",
            f"    C           {c:.6g};",
            f"    rho         {rho:g};",
            f"    value       uniform {initial_value:.6g};",
            "}",
            "",
        ]

    # Recombine the outlets in parallel and compare with the requested totals.
    checks = [
        ("Proximal resistance Rp", "Pa s/m^3", args.proximal_resistance,
         1.0 / sum(1.0 / rch for _, rch, _, _ in written)),
        ("Distal resistance Rd", "Pa s/m^3", args.distal_resistance,
         1.0 / sum(1.0 / r for _, _, r, _ in written)),
        ("Compliance C", "m^3/Pa", args.compliance,
         sum(c for _, _, _, c in written)),
        ("Sum of flow fractions", "", 1.0,
         sum(fraction for fraction, _, _, _ in written)),
    ]
    print()
    print("Check: totals recombined from the per-outlet values")
    print(f"{'Quantity':<24} {'requested':>12} {'recombined':>12} {'rel. error':>11}")
    failed = []
    for label, unit, requested, recombined in checks:
        error = abs(recombined - requested) / requested
        print(f"{label:<24} {requested:>12.6g} {recombined:>12.6g} {error:>11.1e}"
              + (f"  {unit}" if unit else ""))
        if error > CHECK_TOLERANCE:
            failed.append(label)
    if failed:
        raise RuntimeError(
            "recombined totals do not match the requested ones: " + ", ".join(failed)
        )
    print(f"All within {CHECK_TOLERANCE:g}: the outlets reproduce the requested totals.")
    print()

    output = args.output.expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {len(outlets)} Windkessel outlet(s) to {output}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
