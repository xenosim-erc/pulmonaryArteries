#!/usr/bin/env python3
"""Extract the results of a finished FSI run case into fsiScript/results/<case>/.

Reads the run case in place (decomposed or reconstructed; nothing is written
into it) and writes:

  metadata.json        provenance: settings, period, density, solver versions,
                       time window, how TAWSS/OSI were computed, warnings
  timeseries.csv       every time step of every cycle: inlet/lung/outlet flow
                       and pressure, fluid volume, kinetic energy, viscous
                       dissipation, helicity, FSI iterations, ...
  wall_cycleAvg.vtp    last-cycle wall fields on the reference (t = 0) wall:
                       TAWSS, OSI, RRT, ECAP, mean WSS vector, mean and
                       amplitude of pressure, displacement statistics
  wall_phases.pvd      the wall at every --phase-stride-th write time of the
  wall_phases/*.vtp    last cycle and at peak systole and end-diastole (WSS,
                       pressure, displacement), indexed by time in the cycle
  solid_{innerWall,outerWall}_{peakSystole,endDiastole}.vtp
                       the solid wall surfaces (D, sigmaEq, sigma, largest
                       principal Green strain)
  summary.json         the scalar results of this case
  setup/, logs/        copies of the input dictionaries and compressed logs

and adds or replaces this case's row in fsiScript/results/summary.csv.

Units: pressures in the solver are kinematic (p/rho) and so are the wall shear
stresses; both are multiplied by rho here. Column names carry their units.

Last-cycle TAWSS and OSI come from the OSI function object's running integrals
of WSS and |WSS| at the two ends of the last cycle, which is exact when a write
time falls on both ends (meshing.sh sets this up). Otherwise they are
integrated from the WSS snapshots with the trapezoidal rule and metadata.json
says so.

Usage: extract_results.py RUN_CASE [RUN_CASE ...] [options]
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import vtk
from vtk.util import numpy_support as vnp

SCRIPT_DIR = Path(__file__).resolve().parent
FSI_DIR = SCRIPT_DIR.parent
DEFAULT_RESULTS = FSI_DIR / "results"
MMHG = 133.322387415  # Pa
NUMBER = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")
trapz = getattr(np, "trapezoid", None) or np.trapz  # NumPy 2 renamed trapz


# --------------------------------------------------------------------------
# Case settings
# --------------------------------------------------------------------------

def read_params(case: Path) -> dict:
    """meshing.params as a dict: the '# Command:' and '# Repository commit:'
    header lines and every 'Name:   value' setting."""
    path = case / "meshing.params"
    params: dict = {}
    if not path.is_file():
        return params
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        header = re.match(r"^#\s*(Command|Repository commit|Written by meshing\.sh on)\s*:?\s*(.*)$", line)
        if header:
            params[header.group(1)] = header.group(2).strip()
            continue
        setting = re.match(r"^([^#\s][^:]*):\s+(.*)$", line)
        if setting:
            params[setting.group(1).strip()] = setting.group(2).strip()
    return params


def foam_entry(path: Path, key: str) -> str | None:
    """Value of a top-level 'key value;' entry in an OpenFOAM dictionary."""
    if not path.is_file():
        return None
    match = re.search(rf"^\s*{re.escape(key)}\s+([^;]*);", path.read_text(encoding="utf-8", errors="replace"), re.M)
    return match.group(1).strip() if match else None


def fluid_density(case: Path) -> float:
    value = foam_entry(case / "constant/fluid/transportProperties", "rho")
    if value is None:
        raise ValueError("no rho in constant/fluid/transportProperties")
    return float(value.split()[-1])


def inlet_waveform(case: Path) -> Path | None:
    text = (case / "0/fluid/U").read_text(encoding="utf-8", errors="replace") if (case / "0/fluid/U").is_file() else ""
    match = re.search(r'"<constant>/fluid/inletFlow/([^"]+)"', text)
    if match:
        return case / "constant/fluid/inletFlow" / match.group(1)
    return None


def cycle_period(case: Path) -> float:
    """Last time in the inlet waveform table, as in meshing.sh."""
    waveform = inlet_waveform(case)
    if waveform is None or not waveform.is_file():
        raise ValueError("cannot find the inlet waveform; give --period")
    last = None
    for line in waveform.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.replace("(", " ( ").replace(")", " ) ").split()
        if len(parts) >= 3 and parts[0] == "(" and NUMBER.match(parts[1]):
            last = float(parts[1])
    if last is None:
        raise ValueError(f"no times in {waveform}; give --period")
    return last


def solver_versions(case: Path) -> dict:
    """OpenFOAM build and solids4foam version from the head of the solver log."""
    log = case / "log.solids4Foam"
    versions = {}
    if log.is_file():
        with log.open(encoding="utf-8", errors="replace") as handle:
            head = "".join(next(handle, "") for _ in range(60))
        for key, pattern in (("openfoamBuild", r"Build\s*:\s*(.+)"),
                             ("openfoamVersion", r"Version:\s*(v?\d[\w.]*)"),
                             ("solids4foam", r"solids4foam.*?[Vv]ersion:\s*(\S+)")):
            match = re.search(pattern, head)
            if match:
                versions[key] = match.group(1).strip()
    return versions


def git_commit() -> str | None:
    try:
        commit = subprocess.run(["git", "-C", str(FSI_DIR), "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(FSI_DIR), "status", "--porcelain", "--", "pythonScripts"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return commit + (" (uncommitted changes in pythonScripts)" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return None


# --------------------------------------------------------------------------
# Function-object histories
# --------------------------------------------------------------------------

def read_dat(path: Path) -> tuple[list[str], np.ndarray]:
    """Column names (from the last header line) and numeric rows of a
    function-object .dat file. Rows that do not parse are skipped."""
    names: list[str] = []
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        fields = stripped.lstrip("#").split()
        if stripped.startswith("#") or not NUMBER.match(fields[0]):
            if fields and fields[0] == "Time":
                names = fields
            continue
        try:
            row = [float(x) for x in fields]
        except ValueError:
            continue
        if rows and len(row) != len(rows[0]):
            continue
        rows.append(row)
    return names, np.array(rows) if rows else np.empty((0, max(len(names), 1)))


def history(directory: Path, filename: str) -> tuple[list[str], np.ndarray] | None:
    """A function object's history over every restart: the files in its
    <startTime>/ subdirectories joined in time order, a later restart
    replacing any times it repeats."""
    if not directory.is_dir():
        return None
    starts = sorted((d for d in directory.iterdir() if d.is_dir() and NUMBER.match(d.name)),
                    key=lambda d: float(d.name))
    names, pieces = [], []
    stem, suffix = Path(filename).stem, Path(filename).suffix
    for start in starts:
        # A restart from the same start time writes <stem>_<restartTime>.dat
        # beside the original file; read them in time order, original first.
        later = [f for f in start.glob(f"{stem}_*{suffix}") if NUMBER.match(f.stem[len(stem) + 1:])]
        for path in [start / filename, *sorted(later, key=lambda f: float(f.stem[len(stem) + 1:]))]:
            if not path.is_file():
                continue
            names_here, data = read_dat(path)
            if len(data):
                names = names_here or names
                pieces.append(data)
    if not pieces:
        return None
    ncol = min(p.shape[1] for p in pieces)
    data = np.vstack([p[:, :ncol] for p in pieces])
    # Keep the last occurrence of each time (restarts rewrite times).
    keys = np.round(data[:, 0], 9)
    _, last = np.unique(keys[::-1], return_index=True)
    data = data[::-1][np.sort(last)][::-1]
    data = data[np.argsort(data[:, 0], kind="stable")]
    return names[:ncol], data


def restart_times(case: Path) -> list[float]:
    """Times the solver was (re)started at, from the <startTime>/ directories
    and <name>_<time>.dat files function objects write on each start."""
    post = case / "postProcessing"
    times: set[float] = set()
    if not post.is_dir():
        return []
    for directory in [post, *[d for d in post.rglob("*") if d.is_dir()]]:
        if NUMBER.match(directory.name):
            times.add(round(float(directory.name), 9))
            for f in directory.glob("*_*.dat"):
                tail = f.stem.rsplit("_", 1)[-1]
                if NUMBER.match(tail):
                    times.add(round(float(tail), 9))
    return sorted(times)


def time_series(case: Path, rho: float) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Every available history as name -> (time, value), in output units."""
    post = case / "postProcessing"
    series: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    fluid = post / "fluid"
    if fluid.is_dir():
        for fo in sorted(fluid.iterdir()):
            name = fo.name
            if name.endswith("Flow") or name.endswith("Pressure"):
                result = history(fo, "surfaceFieldValue.dat")
                if result is None:
                    continue
                _, data = result
                where = name[:-4] if name.endswith("Flow") else name[:-8]
                if name.endswith("Flow"):
                    series[f"Q_{where}_mLps"] = (data[:, 0], data[:, 1] * 1e6)
                else:
                    series[f"p_{where}_mmHg"] = (data[:, 0], data[:, 1] * rho / MMHG)
        result = history(fluid / "fluidMetrics", "fluidMetrics.dat")
        if result is not None:
            names, data = result
            scale = {"volume": ("volume_mL", 1e6), "kineticEnergy": ("kineticEnergy_mJ", rho * 1e3),
                     "dissipation": ("dissipation_mW", rho * 1e3),
                     "helicity": ("helicity_m4ps2", 1.0), "absHelicity": ("absHelicity_m4ps2", 1.0),
                     "maxMagU": ("maxMagU_mps", 1.0)}
            for i, col in enumerate(names[1:], start=1):
                if col in scale:
                    label, factor = scale[col]
                    series[label] = (data[:, 0], data[:, i] * factor)
        result = history(fluid / "fluidMetrics", "fluidEnergyFlux.dat")
        if result is not None:
            names, data = result
            for i, patch in enumerate(names[1:], start=1):
                series[f"energyOut_{patch}_mW"] = (data[:, 0], data[:, i] * rho * 1e3)
    # Region-less function objects (solid and FSI) write to postProcessing/<startTime>/.
    for filename, columns in (
        ("fsiConvergenceData.dat", {"nFsiCorrectors": ("nFsiCorrectors", 1.0)}),
        ("solidForcesinnerWall.dat", {"normalForce": ("innerWall_normalForce_N", 1.0)}),
    ):
        result = history(post, filename)
        if result is None:
            continue
        names, data = result
        for i, col in enumerate(names[1:], start=1):
            if col in columns:
                label, factor = columns[col]
                series[label] = (data[:, 0], data[:, i] * factor)
    result = history(post, "solidDisplacementsinnerWall.dat")
    if result is not None:
        names, data = result
        if data.shape[1] >= 7:
            extreme = np.abs(data[:, 1:7]).max(axis=1)
            series["innerWall_maxDispComponent_mm"] = (data[:, 0], extreme * 1e3)
    # Robin coupling: outer correctors per step from fsiResiduals.dat when the
    # convergence function object did not write.
    if "nFsiCorrectors" not in series and (post / "fsiResiduals.dat").is_file():
        _, data = read_dat(post / "fsiResiduals.dat")
        if len(data):
            keys = np.round(data[:, 0], 9)
            times, idx = np.unique(keys, return_inverse=True)
            counts = np.zeros(len(times))
            np.maximum.at(counts, idx, data[:, 1])
            series["nFsiCorrectors"] = (times, counts)
    return series


def tabulate(series: dict) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """All histories on the union of their times (NaN where one is missing)."""
    if not series:
        return np.empty(0), {}
    times = np.unique(np.concatenate([np.round(t, 9) for t, _ in series.values()]))
    table = {}
    for name, (t, v) in series.items():
        column = np.full(len(times), np.nan)
        column[np.searchsorted(times, np.round(t, 9))] = v
        table[name] = column
    return times, table


# --------------------------------------------------------------------------
# Cycle statistics
# --------------------------------------------------------------------------

def window(t: np.ndarray, v: np.ndarray, start: float, end: float, tol: float):
    keep = (t >= start - tol) & (t <= end + tol) & np.isfinite(v)
    return t[keep], v[keep]


def cycle_mean(t, v, period):
    return float(trapz(v, t) / period) if len(t) > 1 else math.nan


def periodicity(t, v, end, period, tol):
    """Largest difference between the last two cycles, relative to the range
    of the last one (0 = periodic)."""
    t1, v1 = window(t, v, end - period, end, tol)
    t0, v0 = window(t, v, end - 2 * period, end - period, tol)
    if len(t0) < 3 or len(t1) < 3 or t0[0] > end - 2 * period + 2 * tol:
        return math.nan
    previous = np.interp(t1 - period, t0, v0)
    spread = np.ptp(v1)
    return float(np.max(np.abs(v1 - previous)) / spread) if spread > 0 else math.nan


def cycle_summary(times, table, period, end, tol) -> tuple[dict, dict]:
    """Scalar results over the last cycle, and the times of peak systole and
    end-diastole."""
    out: dict = {}
    events: dict = {}

    def last(name):
        return window(times, table[name], end - period, end, tol) if name in table else (None, None)

    for name in sorted(table):
        if name.startswith("p_") and name.endswith("_mmHg"):
            where = name[2:-5]
            t, v = last(name)
            if len(t) < 2:
                continue
            out[f"p_{where}_mean_mmHg"] = cycle_mean(t, v, period)
            out[f"p_{where}_sys_mmHg"] = float(v.max())
            out[f"p_{where}_dia_mmHg"] = float(v.min())
            out[f"p_{where}_pulse_mmHg"] = float(v.max() - v.min())
        elif name.startswith("Q_") and name.endswith("_mLps"):
            where = name[2:-5]
            t, v = last(name)
            if len(t) < 2:
                continue
            if where == "inlet":
                v = -v  # report inflow as positive
            out[f"Q_{where}_mean_mLps"] = cycle_mean(t, v, period)
            out[f"Q_{where}_peak_mLps"] = float(v.max())
            total = trapz(np.abs(v), t)
            out[f"Q_{where}_reverseFraction"] = float(trapz(np.clip(-v, 0, None), t) / total) if total > 0 else math.nan

    # Cardiac events from the inlet: peak systole at peak inflow, end-diastole
    # at minimum inlet pressure (falling back on the inflow minimum).
    t, q = last("Q_inlet_mLps")
    if t is not None and len(t):
        events["peakSystole"] = float(t[np.argmin(q)])  # inflow is negative
    t, p = last("p_inlet_mmHg")
    if t is not None and len(t) > 1:
        events["endDiastole"] = float(t[np.argmin(p)])
        events.setdefault("peakSystole", float(t[np.argmax(p)]))
    elif "Q_inlet_mLps" in table:
        t, q = last("Q_inlet_mLps")
        if len(t):
            events["endDiastole"] = float(t[np.argmax(q)])

    if {"Q_LPA_mean_mLps", "Q_RPA_mean_mLps"} <= out.keys():
        total = out["Q_LPA_mean_mLps"] + out["Q_RPA_mean_mLps"]
        out["LPA_fraction_mean"] = out["Q_LPA_mean_mLps"] / total if total else math.nan
        if "peakSystole" in events:
            i = int(np.argmin(np.abs(times - events["peakSystole"])))
            left, right = table["Q_LPA_mLps"][i], table["Q_RPA_mLps"][i]
            out["LPA_fraction_peakSystole"] = float(left / (left + right)) if left + right else math.nan
    if {"Q_inlet_mLps", "Q_LPA_mLps", "Q_RPA_mLps"} <= table.keys() and out.get("Q_inlet_mean_mLps"):
        # With a compliant wall, inflow - outflow is the rate the lumen fills
        # (dV/dt), not an error: in + out + dV/dt = 0 (inflow is negative).
        net = table["Q_inlet_mLps"] + table["Q_LPA_mLps"] + table["Q_RPA_mLps"]
        t, stored = window(times, net, end - period, end, tol)
        if len(t) > 1:
            out["netStorage_fractionOfInflow"] = -cycle_mean(t, stored, period) / out["Q_inlet_mean_mLps"]
        if "volume_mL" in table:
            ok = np.isfinite(table["volume_mL"])
            if ok.sum() > 2:
                dvdt = np.interp(times, times[ok], np.gradient(table["volume_mL"][ok], times[ok]))
                t, imbalance = window(times, net + dvdt, end - period, end, tol)
                if len(t) > 1:
                    out["massBalanceError"] = cycle_mean(t, np.abs(imbalance), period) / out["Q_inlet_mean_mLps"]
    if "p_inlet_mean_mmHg" in out and out.get("Q_inlet_mean_mLps"):
        out["resistance_inlet_Pasm3"] = out["p_inlet_mean_mmHg"] * MMHG / (out["Q_inlet_mean_mLps"] * 1e-6)
    for side in ("LPA", "RPA"):
        if f"p_{side}_mean_mmHg" in out and "p_inlet_mean_mmHg" in out:
            out[f"dp_inlet_{side}_mean_mmHg"] = out["p_inlet_mean_mmHg"] - out[f"p_{side}_mean_mmHg"]

    t, v = last("dissipation_mW")
    if t is not None and len(t) > 1:
        out["energyLoss_mJ"] = float(trapz(v, t))  # mW s = mJ
        out["dissipation_mean_mW"] = cycle_mean(t, v, period)
    for name, label in (("kineticEnergy_mJ", "kineticEnergy"), ("helicity_m4ps2", "helicity"),
                        ("absHelicity_m4ps2", "absHelicity")):
        t, v = last(name)
        if t is not None and len(t) > 1:
            unit = name.split("_", 1)[1]
            out[f"{label}_mean_{unit}"] = cycle_mean(t, v, period)
            out[f"{label}_max_{unit}"] = float(v.max())
    t, v = last("volume_mL")
    if t is not None and len(t) > 1:
        out["volume_min_mL"] = float(v.min())
        out["volume_max_mL"] = float(v.max())
        out["volume_stroke_mL"] = float(v.max() - v.min())
    t, v = last("maxMagU_mps")
    if t is not None and len(t):
        out["maxMagU_mps"] = float(v.max())
    t, v = last("nFsiCorrectors")
    if t is not None and len(t):
        out["nFsiCorrectors_mean"] = float(v.mean())
        out["nFsiCorrectors_max"] = float(v.max())

    out.update(power_summary(times, table, period, end, tol, out))

    for name in ("p_inlet_mmHg", "Q_LPA_mLps", "Q_RPA_mLps", "volume_mL"):
        if name in table:
            out[f"periodicity_{name.rsplit('_', 1)[0]}"] = periodicity(times, table[name], end, period, tol)
    return out, events


def power_summary(times, table, period, end, tol, cycle: dict) -> dict:
    """Hydraulic power the inflow delivers at the inlet, split into steady
    and oscillatory parts, and the energy budget of the fluid."""
    out: dict = {}
    if not {"Q_inlet_mLps", "p_inlet_mmHg"} <= table.keys():
        return out
    q = -table["Q_inlet_mLps"] * 1e-6          # inflow [m^3/s], positive
    p = table["p_inlet_mmHg"] * MMHG           # [Pa]
    t, qq = window(times, q, end - period, end, tol)
    _, pq = window(times, p * q, end - period, end, tol)
    _, pp = window(times, p, end - period, end, tol)
    if len(t) < 3 or len(pq) != len(t) or len(pp) != len(t):
        return out
    total = cycle_mean(t, pq, period)          # W
    steady = cycle_mean(t, pp, period) * cycle_mean(t, qq, period)
    out["power_total_mW"] = total * 1e3
    out["power_steady_mW"] = steady * 1e3
    out["power_oscillatory_mW"] = (total - steady) * 1e3
    out["oscillatoryPowerFraction"] = (total - steady) / total if total else math.nan
    stroke_volume = trapz(qq, t)               # m^3
    out["strokeVolume_mL"] = stroke_volume * 1e6
    out["strokeWork_mJ"] = total * period * 1e3
    if stroke_volume > 0:
        # Work per volume ejected: the flow-weighted mean inlet pressure.
        out["workPerStrokeVolume_mmHg"] = total * period / stroke_volume / MMHG
    if "energyLoss_mJ" in cycle and out["strokeWork_mJ"]:
        out["energyLoss_fractionOfStrokeWork"] = cycle["energyLoss_mJ"] / out["strokeWork_mJ"]

    # Energy budget from the patch fluxes (fluidEnergyFlux.dat): the inlet's
    # total (pressure + kinetic) energy inflow goes to the outlets, the wall,
    # viscous dissipation and the change in kinetic energy.
    flux = {n[len("energyOut_"):-3]: n for n in table if n.startswith("energyOut_")}
    if "inlet" in flux:
        def mean(name):
            tt, v = window(times, table[name], end - period, end, tol)
            return cycle_mean(tt, v, period) if len(tt) > 1 else math.nan
        inflow = -mean(flux["inlet"])
        outlets = sum(mean(n) for k, n in flux.items() if k != "inlet" and not k.startswith("wall"))
        wall = sum(mean(n) for k, n in flux.items() if k.startswith("wall"))
        out["power_totalWithKinetic_mW"] = inflow
        if inflow:
            out["kineticPowerFraction"] = (inflow - out["power_total_mW"]) / inflow
            out["energyToOutlets_fraction"] = outlets / inflow
            out["energyToWall_fraction"] = wall / inflow
        change = 0.0
        if "kineticEnergy_mJ" in table:
            tt, ke = window(times, table["kineticEnergy_mJ"], end - period, end, tol)
            change = (ke[-1] - ke[0]) / period if len(tt) > 1 else 0.0
        if "dissipation_mean_mW" in cycle and inflow:
            out["energyLoss_fractionOfInflow"] = cycle["dissipation_mean_mW"] / inflow
            out["energyBudgetResidual_fraction"] = (inflow - outlets - wall - cycle["dissipation_mean_mW"]
                                                    - change) / inflow
    return out


def windkessel_baseline(case: Path, times, table, period, end, tol) -> dict:
    """The same inlet power measures for the outlet Windkessels alone (no 3D
    arteries): the case's per-outlet RCR models in parallel, driven by its
    inflow. The FSI values minus these are the effect of the 3D arteries."""
    path = case / "0/fluid/windkesselOutlets"
    if not path.is_file():
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")
    outlets = []
    for block in re.finditer(r"(?:^|\n)\s*(\w+)\s*\n\s*\{(.*?)\}", text, re.S):
        body = block.group(2)
        values = {k: re.search(rf"\b{k}\s+([-+0-9.eE]+)\s*;", body) for k in ("R", "Rch", "C")}
        if all(values.values()):
            outlets.append({k: float(v.group(1)) for k, v in values.items()})
    if not outlets:
        return {}
    if "Q_inlet_mLps" in table:
        t, q = window(times, -table["Q_inlet_mLps"] * 1e-6, end - period, end, tol)
        source = "the case's inlet flow over the last cycle"
    else:
        t, q = np.empty(0), np.empty(0)
    if len(t) < 8:
        waveform = inlet_waveform(case)
        if waveform is None or not waveform.is_file():
            return {}
        rows = [line.replace("(", " ").replace(")", " ").split() for line in
                waveform.read_text(encoding="utf-8", errors="replace").splitlines()]
        rows = [r for r in rows if len(r) >= 2 and NUMBER.match(r[0]) and NUMBER.match(r[1])]
        t, q = np.array([float(r[0]) for r in rows]), np.array([float(r[1]) for r in rows])
        t = t - t[0]
        source = "the inlet waveform"
    n = 4096
    grid = t[0] + np.arange(n) * period / n
    flow = np.interp(grid, t, q, period=period)
    harmonics = np.fft.rfft(flow) / n
    omega = 2 * np.pi * np.arange(len(harmonics)) / period
    admittance = sum(1.0 / (o["Rch"] + o["R"] / (1 + 1j * omega * o["R"] * o["C"])) for o in outlets)
    pressure = np.fft.irfft(harmonics / admittance * n, n)
    total = float(np.mean(pressure * flow))
    steady = float(pressure.mean() * flow.mean())
    return {
        "baseline_power_total_mW": total * 1e3,
        "baseline_oscillatoryPowerFraction": (total - steady) / total if total else math.nan,
        "baseline_p_mean_mmHg": float(pressure.mean() / MMHG),
        "baseline_p_pulse_mmHg": float(np.ptp(pressure) / MMHG),
        "_source": source,
    }


def write_pv(path: Path, times, table, period, end, tol) -> bool:
    """Last-cycle inlet pressure against ejected volume (the ejection limb of
    a right-ventricular PV loop) and against the lumen volume of the 3D model."""
    if not {"Q_inlet_mLps", "p_inlet_mmHg"} <= table.keys():
        return False
    keep = (times >= end - period - tol) & (times <= end + tol)
    keep &= np.isfinite(table["Q_inlet_mLps"]) & np.isfinite(table["p_inlet_mmHg"])
    if keep.sum() < 3:
        return False
    t = times[keep]
    q = -table["Q_inlet_mLps"][keep]
    ejected = np.concatenate([[0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(t))])
    volume = table["volume_mL"][keep] if "volume_mL" in table else np.full(len(t), np.nan)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s", "phase_s", "p_inlet_mmHg", "Q_inflow_mLps", "V_ejected_mL",
                         "volume_mL", "dVolume_mL"])
        for i in range(len(t)):
            writer.writerow([f"{t[i]:.6g}", f"{t[i] - t[0]:.6g}", f"{table['p_inlet_mmHg'][keep][i]:.8g}",
                             f"{q[i]:.8g}", f"{ejected[i]:.8g}",
                             "" if not np.isfinite(volume[i]) else f"{volume[i]:.8g}",
                             "" if not np.isfinite(volume[i] - volume[0]) else f"{volume[i] - volume[0]:.8g}"])
    return True


# --------------------------------------------------------------------------
# Fields
# --------------------------------------------------------------------------

class FoamCase:
    """Thin wrapper round vtkPOpenFOAMReader for one patch or mesh at a time."""

    def __init__(self, case: Path):
        self.reader = vtk.vtkPOpenFOAMReader()
        self.reader.SetFileName(str(case / "results.foam"))  # need not exist
        self.reader.SetCaseType(0 if (case / "processor0").is_dir() else 1)
        self.reader.CreateCellToPointOff()
        self.reader.SetUse64BitFloats(1)
        self.reader.SetUse64BitLabels(0)
        self.reader.UpdateInformation()
        values = self.reader.GetTimeValues()
        self.times = np.array([values.GetValue(i) for i in range(values.GetNumberOfTuples())]) \
            if values is not None else np.empty(0)
        self.patches = [self.reader.GetPatchArrayName(i) for i in range(self.reader.GetNumberOfPatchArrays())]
        self.cell_arrays = [self.reader.GetCellArrayName(i) for i in range(self.reader.GetNumberOfCellArrays())]

    def read(self, block: str, time: float, arrays: list[str]):
        r = self.reader
        r.DisableAllPatchArrays()
        r.SetPatchArrayStatus(block, 1)
        r.DisableAllCellArrays()
        for name in arrays:
            r.SetCellArrayStatus(name, 1)
        r.DisableAllPointArrays()
        r.UpdateTimeStep(time)
        r.Update()
        leaves = []
        it = r.GetOutput().NewIterator()
        it.InitTraversal()
        while not it.IsDoneWithTraversal():
            leaf = it.GetCurrentDataObject()
            if leaf is not None and leaf.GetNumberOfCells() > 0:
                copy = leaf.NewInstance()
                copy.DeepCopy(leaf)
                leaves.append(copy)
            it.GoToNextItem()
        if not leaves:
            raise ValueError(f"{block} is empty at t = {time}")
        if len(leaves) == 1:
            return leaves[0]
        append = vtk.vtkAppendPolyData() if leaves[0].IsA("vtkPolyData") else vtk.vtkAppendFilter()
        for leaf in leaves:
            append.AddInputData(leaf)
        append.Update()
        return append.GetOutput()


def cell_array(data, name):
    array = data.GetCellData().GetArray(name)
    return None if array is None else vnp.vtk_to_numpy(array).astype(float)


def points_of(data) -> np.ndarray:
    return vnp.vtk_to_numpy(data.GetPoints().GetData()).astype(float)


def face_areas(poly) -> np.ndarray:
    sizes = vtk.vtkCellSizeFilter()
    sizes.SetInputData(poly)
    sizes.ComputeAreaOn()
    sizes.ComputeLengthOff()
    sizes.ComputeVolumeOff()
    sizes.ComputeVertexCountOff()
    sizes.Update()
    return vnp.vtk_to_numpy(sizes.GetOutput().GetCellData().GetArray("Area")).astype(float)


def add_array(data, name, values, point=False):
    # float32 on disk: everything is computed in float64 first.
    values = np.ascontiguousarray(values, dtype=np.float32)
    array = vnp.numpy_to_vtk(values, deep=True)
    array.SetName(name)
    (data.GetPointData() if point else data.GetCellData()).AddArray(array)


def write_xml(data, path: Path):
    writer = vtk.vtkXMLPolyDataWriter() if data.IsA("vtkPolyData") else vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(data)
    writer.SetDataModeToAppended()
    writer.SetCompressorTypeToZLib()
    writer.Write()


def surface(reference_poly):
    """Copy of a surface's geometry (float32 points) with no arrays, to fill
    with results."""
    points = vtk.vtkPoints()
    points.SetData(vnp.numpy_to_vtk(points_of(reference_poly).astype(np.float32), deep=True))
    poly = vtk.vtkPolyData()
    poly.SetPoints(points)
    poly.SetPolys(reference_poly.GetPolys())
    return poly


def point_to_face(poly, values: np.ndarray) -> np.ndarray:
    """Mean of each face's point values: lets point data be area-weighted
    without counting the points duplicated along processor boundaries."""
    polys = poly.GetPolys()
    conn = vnp.vtk_to_numpy(polys.GetConnectivityArray()).astype(np.int64)
    offsets = vnp.vtk_to_numpy(polys.GetOffsetsArray()).astype(np.int64)
    return np.add.reduceat(values[conn], offsets[:-1]) / np.diff(offsets)


def weighted_percentile(values, weights, q):
    order = np.argsort(values)
    cum = np.cumsum(weights[order])
    return float(np.interp(q / 100 * cum[-1], cum, values[order]))


def find_array(names: list[str], suffix: str) -> str | None:
    hits = [n for n in names if n.endswith(suffix)]
    return hits[0] if hits else None


def extract_wall(foam: FoamCase, case: Path, out: Path, rho: float, period: float, end: float,
                 tol: float, events: dict, args, meta: dict) -> dict:
    """Wall fields over the last cycle; returns summary scalars."""
    summary: dict = {}
    wall = "/fluid/patch/wall"
    if wall not in foam.patches:
        meta["warnings"].append("no fluid wall patch: wall results skipped")
        return summary
    acc_vec = find_array(foam.cell_arrays, "accumulatedWallShearStress")
    acc_mag = find_array(foam.cell_arrays, "accumulatedMagWallShearStress")
    snapshot_times = foam.times[(foam.times >= end - period - tol) & (foam.times <= end + tol)]
    if len(snapshot_times) < 2:
        meta["warnings"].append("fewer than two write times in the last cycle: wall results skipped")
        return summary
    t0, t1 = snapshot_times[0], snapshot_times[-1]
    full_cycle = abs((t1 - t0) - period) <= tol

    reference = foam.read(wall, float(foam.times[0]), [])
    x_ref = points_of(reference)
    area = face_areas(reference)

    wss, pressure, disp = [], [], []
    for t in snapshot_times:
        data = foam.read(wall, float(t), ["wallShearStress", "p"])
        w = cell_array(data, "wallShearStress")
        p = cell_array(data, "p")
        if w is None or p is None:
            raise ValueError(f"wallShearStress or p missing on the wall at t = {t}")
        wss.append(w * rho)
        pressure.append(p * rho)
        disp.append(points_of(data) - x_ref)
    wss, pressure, disp = np.array(wss), np.array(pressure), np.array(disp)

    # The OSI function object restarts its integrals from zero whenever the
    # solver restarts. A restart before the cycle is harmless; one exactly at
    # its start means the integral at the end covers the cycle alone; one
    # inside it leaves only the snapshots.
    restarts = [r for r in restart_times(case) if r > tol]
    inside = [r for r in restarts if t0 + tol < r < t1 - tol]
    at_start = any(abs(r - t0) <= tol for r in restarts)
    if inside:
        meta["warnings"].append(f"solver restart(s) at t = {inside} inside the last cycle reset the "
                                "OSI running integrals")
    method = None
    if full_cycle and acc_vec and acc_mag and not inside:
        stop = foam.read(wall, float(t1), [acc_vec, acc_mag])
        e_vec, e_mag = cell_array(stop, acc_vec), cell_array(stop, acc_mag)
        if at_start:
            s_vec, s_mag = np.zeros_like(e_vec), np.zeros_like(e_mag)
        else:
            start = foam.read(wall, float(t0), [acc_vec, acc_mag])
            s_vec, s_mag = cell_array(start, acc_vec), cell_array(start, acc_mag)
        if all(a is not None for a in (s_vec, s_mag, e_vec, e_mag)):
            mean_wss = (e_vec - s_vec) * rho / period
            tawss = (e_mag - s_mag) * rho / period
            method = ("running integrals of the OSI function object at both ends of the cycle (exact)"
                      + ("; restarted at the cycle start" if at_start else ""))
    if method is None:
        span = t1 - t0
        mean_wss = trapz(wss, snapshot_times, axis=0) / span
        tawss = trapz(np.linalg.norm(wss, axis=2), snapshot_times, axis=0) / span
        method = (f"trapezoidal rule over {len(snapshot_times)} WSS snapshots spanning {span:g} s "
                  f"of a {period:g} s cycle (approximate)")
        meta["warnings"].append("TAWSS/OSI from snapshots, not the running integrals: " + method)
    meta["tawssMethod"] = method
    meta["wallWindow"] = {"start": float(t0), "end": float(t1), "snapshots": len(snapshot_times),
                          "fullCycle": bool(full_cycle)}

    mag_mean = np.linalg.norm(mean_wss, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        osi = np.clip(0.5 * (1 - mag_mean / tawss), 0, 0.5)
        rrt = 1.0 / mag_mean
        ecap = osi / tawss

    def nearest(name):
        if name not in events:
            return None
        return int(np.argmin(np.abs(snapshot_times - events[name])))

    i_ps, i_ed = nearest("peakSystole"), nearest("endDiastole")
    if i_ed is None:
        i_ed = 0
    disp_mag = np.linalg.norm(disp, axis=2)
    cyclic = np.linalg.norm(disp - disp[i_ed], axis=2).max(axis=0)

    poly = surface(reference)
    add_array(poly, "TAWSS_Pa", tawss)
    add_array(poly, "OSI", osi)
    add_array(poly, "RRT_1pPa", rrt)
    add_array(poly, "ECAP_1pPa", ecap)
    add_array(poly, "meanWSS_Pa", mean_wss)
    add_array(poly, "p_mean_mmHg", trapz(pressure, snapshot_times, axis=0) / (t1 - t0) / MMHG)
    add_array(poly, "p_amplitude_mmHg", np.ptp(pressure, axis=0) / MMHG)
    add_array(poly, "area_m2", area)
    add_array(poly, "dispMag_mean_mm", trapz(disp_mag, snapshot_times, axis=0) / (t1 - t0) * 1e3, point=True)
    add_array(poly, "dispMag_max_mm", disp_mag.max(axis=0) * 1e3, point=True)
    add_array(poly, "dispCyclic_mm", cyclic * 1e3, point=True)
    add_array(poly, "D_endDiastole_m", disp[i_ed], point=True)
    if i_ps is not None:
        add_array(poly, "WSSmag_peakSystole_Pa", np.linalg.norm(wss[i_ps], axis=1))
        add_array(poly, "p_peakSystole_mmHg", pressure[i_ps] / MMHG)
        add_array(poly, "D_peakSystole_m", disp[i_ps], point=True)
    write_xml(poly, out / "wall_cycleAvg.vtp")

    phase_dir = out / "wall_phases"
    phase_dir.mkdir(exist_ok=True)
    entries = []
    # Every phase_stride-th write, plus peak systole and end-diastole; the
    # full cycle is in the archived run.
    keep = set(range(0, len(snapshot_times), max(args.phase_stride, 1))) | {i_ed, len(snapshot_times) - 1}
    if i_ps is not None:
        keep.add(i_ps)
    meta["wallPhaseTimes"] = [float(snapshot_times[k]) for k in sorted(keep)]
    for k in sorted(keep):
        t = snapshot_times[k]
        frame = surface(reference)
        add_array(frame, "WSS_Pa", wss[k])
        add_array(frame, "WSSmag_Pa", np.linalg.norm(wss[k], axis=1))
        add_array(frame, "p_mmHg", pressure[k] / MMHG)
        add_array(frame, "D_m", disp[k], point=True)
        name = f"wall_{k:03d}.vtp"
        write_xml(frame, phase_dir / name)
        entries.append(f'    <DataSet timestep="{t - t0:.6g}" file="wall_phases/{name}"/>')
    (out / "wall_phases.pvd").write_text(
        '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n  <Collection>\n'
        + "\n".join(entries) + "\n  </Collection>\n</VTKFile>\n")

    ok = np.isfinite(tawss) & np.isfinite(osi) & (area > 0)
    a = area[ok]
    summary["wallArea_cm2"] = float(area.sum() * 1e4)
    summary["TAWSS_mean_Pa"] = float(np.average(tawss[ok], weights=a))
    for q in (5, 50, 95):
        summary[f"TAWSS_p{q}_Pa"] = weighted_percentile(tawss[ok], a, q)
    summary["OSI_mean"] = float(np.average(osi[ok], weights=a))
    summary["OSI_p95"] = weighted_percentile(osi[ok], a, 95)
    summary[f"lowTAWSS_areaFraction_lt{args.low_tawss:g}Pa"] = float(a[tawss[ok] < args.low_tawss].sum() / a.sum())
    summary[f"highOSI_areaFraction_gt{args.high_osi:g}"] = float(a[osi[ok] > args.high_osi].sum() / a.sum())
    finite = ok & np.isfinite(ecap)
    summary[f"highECAP_areaFraction_gt{args.high_ecap:g}"] = float(area[finite][ecap[finite] > args.high_ecap].sum() / a.sum())
    summary["RRT_p50_1pPa"] = weighted_percentile(rrt[finite], area[finite], 50)
    if i_ps is not None:
        summary["WSSmag_peakSystole_mean_Pa"] = float(np.average(np.linalg.norm(wss[i_ps], axis=1)[ok], weights=a))
    summary["disp_max_mm"] = float(disp_mag.max() * 1e3)
    summary["dispCyclic_mean_mm"] = float(np.average(point_to_face(reference, cyclic), weights=area) * 1e3)
    summary["dispCyclic_max_mm"] = float(cyclic.max() * 1e3)
    return summary


def green_strain_max(F: np.ndarray) -> np.ndarray:
    """Largest principal Green-Lagrange strain of each deformation gradient."""
    F = F.reshape(-1, 3, 3)
    E = 0.5 * (np.einsum("nki,nkj->nij", F, F) - np.eye(3))
    return np.linalg.eigvalsh(E)[:, -1]


def extract_solid(foam: FoamCase, out: Path, events: dict, meta: dict) -> dict:
    """The solid's inner and outer wall surfaces at peak systole and
    end-diastole (the full volume stays in the archived run)."""
    summary: dict = {}
    wanted = [n for n in ("D", "sigmaEq", "sigma", "F") if n in foam.cell_arrays]
    grads: dict = {}
    for event in ("peakSystole", "endDiastole"):
        if event not in events or not len(foam.times):
            continue
        t = float(foam.times[np.argmin(np.abs(foam.times - events[event]))])
        meta.setdefault("solidTimes", {})[event] = t
        for patch in ("innerWall", "outerWall"):
            block = f"/solid/patch/{patch}"
            if block not in foam.patches:
                continue
            face = foam.read(block, t, wanted)
            poly = surface(face)
            a = face_areas(face)
            add_array(poly, "area_m2", a)
            for name, scale, label in (("D", 1.0, "D_m"), ("sigmaEq", 1e-3, "sigmaEq_kPa"),
                                       ("sigma", 1e-3, "sigma_kPa")):
                values = cell_array(face, name)
                if values is not None:
                    add_array(poly, label, values * scale)
            F = cell_array(face, "F")
            if F is not None and F.shape[1] == 9:
                grads[(event, patch)] = F
                add_array(poly, "greenStrainMaxPrincipal", green_strain_max(F))
            write_xml(poly, out / f"solid_{patch}_{event}.vtp")
            s = cell_array(face, "sigmaEq")
            if s is not None:
                summary[f"sigmaEq_{patch}_{event}_mean_kPa"] = float(np.average(s, weights=a) / 1e3)
                summary[f"sigmaEq_{patch}_{event}_p99_kPa"] = weighted_percentile(s, a, 99) / 1e3
                summary[f"sigmaEq_{patch}_{event}_max_kPa"] = float(s.max() / 1e3)
            if (("peakSystole", patch) in grads and ("endDiastole", patch) in grads):
                # Strain of peak systole relative to end-diastole: F_rel = F_ps F_ed^-1.
                f_ps = grads[("peakSystole", patch)].reshape(-1, 3, 3)
                f_ed = grads[("endDiastole", patch)].reshape(-1, 3, 3)
                cyclic = green_strain_max(np.einsum("nij,njk->nik", f_ps, np.linalg.inv(f_ed)))
                summary[f"cyclicGreenStrain_{patch}_p50"] = weighted_percentile(cyclic, a, 50)
                summary[f"cyclicGreenStrain_{patch}_p99"] = weighted_percentile(cyclic, a, 99)
    return summary


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

def copy_setup(case: Path, out: Path):
    setup = out / "setup"
    setup.mkdir(exist_ok=True)
    for rel in ("meshing.params", "lungSides.csv", "log.windkesselOutlets", "system/controlDict",
                "system/lungFlows", "constant/fluid/transportProperties",
                "constant/solid/mechanicalProperties", "constant/fsiProperties",
                "0/fluid/windkesselOutlets"):
        src = case / rel
        if src.is_file():
            shutil.copy2(src, setup / rel.replace("/", "_"))
    waveform = inlet_waveform(case)
    if waveform and waveform.is_file():
        shutil.copy2(waveform, setup / f"inletFlow_{waveform.name}")
    logs = out / "logs"
    logs.mkdir(exist_ok=True)
    for src in sorted(case.glob("log.*")) + sorted((case / "prestress").glob("log.*")):
        if src.is_file():
            prefix = "prestress_" if src.parent.name == "prestress" else ""
            with src.open("rb") as fin, gzip.open(logs / f"{prefix}{src.name}.gz", "wb") as fout:
                shutil.copyfileobj(fin, fout)


def write_timeseries(path: Path, times, table, period):
    names = sorted(table)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s", "cycle", "phase_s", *names])
        for i, t in enumerate(times):
            cycle = int(math.floor(t / period - 1e-9))
            row = [f"{t:.6g}", cycle + 1, f"{t - cycle * period:.6g}"]
            row += ["" if not np.isfinite(table[n][i]) else f"{table[n][i]:.8g}" for n in names]
            writer.writerow(row)


def update_summary(path: Path, case_name: str, row: dict):
    rows = []
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            rows = [r for r in csv.DictReader(handle) if r.get("case") != case_name]
    rows.append({"case": case_name, **{k: v for k, v in row.items()}})
    columns = ["case"] + sorted({k for r in rows for k in r} - {"case"})
    rows.sort(key=lambda r: r["case"])
    tmp = path.with_suffix(".csv.tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: fmt(r.get(k, "")) for k in columns})
    tmp.replace(path)


def fmt(value):
    if isinstance(value, float):
        return "" if not math.isfinite(value) else f"{value:.6g}"
    return value


def jsonable(value):
    if isinstance(value, dict):
        return {k: jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, (np.floating, float)):
        return None if not math.isfinite(float(value)) else float(value)
    if isinstance(value, np.integer):
        return int(value)
    return value


def extract(case: Path, args) -> dict:
    case = case.resolve()
    name = args.name or case.name
    out = (args.results_dir / name).resolve()
    if out.exists() and not args.force:
        raise FileExistsError(f"{out} exists; pass --force to replace it")
    previous_archive = None
    if out.exists():
        try:  # keep the record archive_run.py added
            previous_archive = json.loads((out / "metadata.json").read_text(encoding="utf-8")).get("archive")
        except (OSError, ValueError):
            pass
    # Build in <name>.partial and swap it in only once everything succeeded,
    # so a failed --force leaves the earlier results untouched.
    final, out = out, out.with_name(out.name + ".partial")
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    try:
        return _extract(case, name, out, final, previous_archive, args)
    except BaseException:
        shutil.rmtree(out, ignore_errors=True)
        raise


def _extract(case: Path, name: str, out: Path, final: Path, previous_archive, args) -> dict:
    params = read_params(case)
    rho = args.rho or fluid_density(case)
    period = args.period or cycle_period(case)
    delta_t = float(foam_entry(case / "system/controlDict", "deltaT") or 0) or None
    end_time = float(foam_entry(case / "system/controlDict", "endTime") or 0) or None
    tol = 0.5 * (delta_t or 1e-6)
    meta = {
        "case": name,
        "runCase": str(case),
        "extracted": dt.datetime.now().isoformat(timespec="seconds"),
        "extractorCommit": git_commit(),
        "period_s": period,
        "rho_kgpm3": rho,
        "deltaT_s": delta_t,
        "endTime_s": end_time,
        "writeInterval_s": foam_entry(case / "system/controlDict", "writeInterval"),
        "purgeWrite": foam_entry(case / "system/controlDict", "purgeWrite"),
        "inletWaveform": str(inlet_waveform(case)) if inlet_waveform(case) else None,
        "solver": solver_versions(case),
        "meshingParams": params,
        "thresholds": {"lowTAWSS_Pa": args.low_tawss, "highOSI": args.high_osi, "highECAP_1pPa": args.high_ecap},
        "warnings": [],
    }
    if previous_archive:
        meta["archive"] = previous_archive

    series = time_series(case, rho)
    times, table = tabulate(series)
    has_mesh = any((case / d / r / "polyMesh").is_dir() for d in ("constant", "processor0/constant")
                   for r in ("fluid", "solid"))
    foam = FoamCase(case) if has_mesh else None
    if foam is None:
        meta["warnings"].append("no mesh found: wall and solid results skipped")
    field_times = foam.times if foam is not None else np.empty(0)
    last_time = float(max(times.max() if len(times) else 0, field_times.max() if len(field_times) else 0))
    meta["lastTime_s"] = last_time
    meta["lastWriteTime_s"] = float(field_times.max()) if len(field_times) else None
    meta["complete"] = bool(end_time is not None and last_time >= end_time - tol)
    meta["cyclesRun"] = last_time / period
    if not meta["complete"]:
        meta["warnings"].append(f"run stopped at t = {last_time:g} s before endTime {end_time} s")
    if last_time < period - tol:
        meta["warnings"].append("less than one full cycle was run")
    end = float(field_times.max()) if len(field_times) and not args.series_end else last_time
    meta["lastCycle_s"] = [end - period, end]

    summary: dict = {}
    events: dict = {}
    if len(times):
        write_timeseries(out / "timeseries.csv", times, table, period)
        cycle, events = cycle_summary(times, table, period, end, tol)
        summary.update(cycle)
        baseline = windkessel_baseline(case, times, table, period, end, tol)
        if baseline:
            meta["windkesselBaselineSource"] = baseline.pop("_source")
            summary.update(baseline)
            if "oscillatoryPowerFraction" in summary:
                summary["oscillatoryPowerFraction_minusBaseline"] = (
                    summary["oscillatoryPowerFraction"] - summary["baseline_oscillatoryPowerFraction"])
        write_pv(out / "pv_inlet.csv", times, table, period, end, tol)
    else:
        meta["warnings"].append("no function-object histories found in postProcessing/")
    if "peakSystole" not in events and len(field_times):
        meta["warnings"].append("no inlet history: cardiac events not known, wall/solid phases use the cycle start")
        events["endDiastole"] = end - period
    meta["events_s"] = events

    if foam is not None:
        summary.update(extract_wall(foam, case, out, rho, period, end, tol, events, args, meta))
        summary.update(extract_solid(foam, out, events, meta))
    copy_setup(case, out)

    (out / "summary.json").write_text(json.dumps(jsonable(summary), indent=2) + "\n", encoding="utf-8")
    (out / "metadata.json").write_text(json.dumps(jsonable(meta), indent=2) + "\n", encoding="utf-8")
    row = {"complete": meta["complete"], "cyclesRun": round(meta["cyclesRun"], 3),
           "animal": params.get("Animal", "").split(" ")[0], "commit": params.get("Repository commit", ""),
           "nWarnings": len(meta["warnings"]), **summary}
    if final.exists():
        shutil.rmtree(final)
    out.rename(final)
    update_summary(args.results_dir / "summary.csv", name, row)
    for warning in meta["warnings"]:
        print(f"  warning: {warning}", file=sys.stderr)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=__doc__.split("\n\n", 1)[1])
    parser.add_argument("cases", nargs="+", type=Path, help="run case directories")
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS,
                        help=f"where results/<case>/ and summary.csv go (default {DEFAULT_RESULTS})")
    parser.add_argument("--name", help="results name, if not the case directory name (one case only)")
    parser.add_argument("--period", type=float, help="cycle length [s] (default: last time of the inlet waveform)")
    parser.add_argument("--rho", type=float, help="fluid density [kg/m^3] (default: transportProperties)")
    parser.add_argument("--low-tawss", type=float, default=0.4, help="low-TAWSS threshold [Pa] (default 0.4)")
    parser.add_argument("--high-osi", type=float, default=0.2, help="high-OSI threshold (default 0.2)")
    parser.add_argument("--high-ecap", type=float, default=1.4, help="high-ECAP threshold [1/Pa] (default 1.4)")
    parser.add_argument("--phase-stride", type=int, default=4,
                        help="keep every Nth write time of the last cycle in wall_phases (default 4)")
    parser.add_argument("--series-end", action="store_true",
                        help="end the last cycle at the last history time rather than the last write time")
    parser.add_argument("--force", action="store_true", help="replace existing results for the case")
    args = parser.parse_args()
    if args.name and len(args.cases) > 1:
        parser.error("--name needs a single case")
    args.results_dir = args.results_dir.resolve()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    failures = 0
    for case in args.cases:
        print(f"Extracting {case}")
        try:
            summary = extract(case, args)
            print(f"  wrote {args.results_dir / (args.name or case.resolve().name)} ({len(summary)} summary values)")
        except Exception as error:  # report and go on to the next case
            print(f"  error: {type(error).__name__}: {error}", file=sys.stderr)
            failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
