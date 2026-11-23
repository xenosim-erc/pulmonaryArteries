#!/usr/bin/env python3
"""Archive the last cycle of a finished FSI run case to the HDD, then remove it.

For each RUN_CASE this

1. checks that the run reached endTime, that no Slurm job is running in it,
   and that results/<case>/ has been extracted from it (extract_results.py);
2. picks the time directories of the last cardiac cycle (both ends) plus 0/,
   running reconstructPar for any of them that are not reconstructed;
3. writes <archive-dir>/<case>.tar.zst holding the reconstructed case:
   setup files, logs, constant/, system/, 0/, the last-cycle time directories
   (without restart-only fields, unless --keep-all-fields), postProcessing/
   (every cycle), the prestress/ and windkessel/ setup and logs, and the
   top-level files of the meshing directory (STL, centerlines, profiles) under
   meshing/. Processor directories, dynamicCode/ and earlier cycles are left
   out;
4. checks the archive by listing it back and comparing every file's size,
   and writes <case>.manifest.json beside it (files, kept times, sha256);
5. records the archive in results/<case>/metadata.json; and
6. deletes the run case, after asking (or straight away with --delete;
   never with --keep).

The meshing directory itself (caseFiles/<case>, intermediate meshes) is not
deleted.

Usage: archive_run.py RUN_CASE [RUN_CASE ...] [options]
"""

from __future__ import annotations

import argparse
import datetime as dt
import fnmatch
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
FSI_DIR = SCRIPT_DIR.parent
DEFAULT_ARCHIVE = Path("/mnt/backup") / os.environ.get("USER", "") / "fsiRuns"
DEFAULT_RESULTS = FSI_DIR / "results"
FOAM_BASHRC = os.environ.get("FOAM_BASHRC", "/usr/lib/openfoam/openfoam2412/etc/bashrc")
NUMBER = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")

# Fields needed only to restart the solver, or recoverable from what is kept
# (mesh motion from polyMesh/points, sigmaHyd from sigma).
RESTART_ONLY = ["Uf", "meshPhi", "cellMotionU", "pointMotionU", "F_", "Ff_", "sigmaHyd", "*_0"]
# Subdirectories of the case that hold set-up and logs; their time and
# processor directories are left out (prestress's result is 0/solid/Fpre).
SIDE_CASES = ("prestress", "windkessel")


def foam_entry(path: Path, key: str) -> str | None:
    if not path.is_file():
        return None
    match = re.search(rf"^\s*{re.escape(key)}\s+([^;]*);",
                      path.read_text(encoding="utf-8", errors="replace"), re.M)
    return match.group(1).strip() if match else None


def cycle_period(case: Path) -> float:
    """Last time of the inlet waveform named in 0/fluid/U (as meshing.sh)."""
    text = (case / "0/fluid/U").read_text(encoding="utf-8", errors="replace")
    match = re.search(r'"<constant>/fluid/inletFlow/([^"]+)"', text)
    if not match:
        raise ValueError("cannot find the inlet waveform in 0/fluid/U; give --period")
    last = None
    waveform = case / "constant/fluid/inletFlow" / match.group(1)
    for line in waveform.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.replace("(", " ( ").replace(")", " ) ").split()
        if len(parts) >= 3 and parts[0] == "(" and NUMBER.match(parts[1]):
            last = float(parts[1])
    if last is None:
        raise ValueError(f"no times in {waveform}; give --period")
    return last


def meshing_dir(case: Path) -> Path | None:
    params = case / "meshing.params"
    if not params.is_file():
        return None
    match = re.search(r"^Meshing directory:\s+(.+)$",
                      params.read_text(encoding="utf-8", errors="replace"), re.M)
    return Path(match.group(1).strip()) if match else None


def time_dirs(directory: Path) -> dict[float, Path]:
    if not directory.is_dir():
        return {}
    return {float(d.name): d for d in directory.iterdir() if d.is_dir() and NUMBER.match(d.name)}


def running_jobs(case: Path) -> list[str]:
    """Slurm jobs of this user whose working directory is the case."""
    try:
        out = subprocess.run(["squeue", "-h", "-u", os.environ.get("USER", ""), "-o", "%i %Z"],
                             capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return []
    jobs = []
    for line in out.splitlines():
        parts = line.split(maxsplit=1)
        if len(parts) == 2 and Path(parts[1]).resolve() == case:
            jobs.append(parts[0])
    return jobs


def reconstruct(case: Path, start: float, end: float):
    """reconstructPar both regions over [start, end] (needs OpenFOAM)."""
    for region in ("fluid", "solid"):
        log = case / f"log.reconstructPar.{region}.archive"
        command = (f'source "{FOAM_BASHRC}" >/dev/null 2>&1; cd "{case}" && '
                   f'reconstructPar -region {region} -time "{start:.10g}:{end:.10g}" > "{log}" 2>&1')
        print(f"  reconstructPar -region {region} for t = {start:g} to {end:g} (log: {log.name})")
        if subprocess.run(["bash", "-c", command]).returncode != 0:
            raise RuntimeError(f"reconstructPar -region {region} failed; see {log}")


def excluded(name: str, patterns: list[str]) -> bool:
    return any(fnmatch.fnmatch(name, p) for p in patterns)


def inside(path: Path, case: Path) -> bool:
    try:
        path.resolve().relative_to(case.resolve())
        return True
    except ValueError:
        return False


def collect(case: Path, keep_times: list[Path], drop_fields: list[str],
            meshing: Path | None) -> list[tuple[Path, str]]:
    """(source file, name in the archive) for everything to archive."""
    name = case.name
    members: list[tuple[Path, str]] = []

    # Symlinks that point inside the case (solids4foam's fsiProperties ->
    # fsiProperties.pimpleFluid, 0/fluid -> fluid.pimpleFluid, ...) are
    # archived as the files they point to; links leaving the case are not.
    def add_tree(root: Path, arc_root: str, skip_dir=lambda d: False, skip_file=lambda f: False):
        for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
            here = Path(dirpath)
            dirnames[:] = sorted(d for d in dirnames if not skip_dir(here / d)
                                 and (not (here / d).is_symlink() or inside(here / d, case)))
            for f in sorted(filenames):
                path = here / f
                if path.is_file() and (not path.is_symlink() or inside(path, case)) and not skip_file(path):
                    members.append((path, f"{arc_root}/{path.relative_to(root)}"))

    kept = {p.resolve() for p in keep_times}
    for entry in sorted(case.iterdir()):
        arc = f"{name}/{entry.name}"
        if entry.is_symlink() and not inside(entry, case):
            continue
        if entry.is_file():
            if not entry.name.endswith(".foam"):
                members.append((entry, arc))
        elif entry.name.startswith("processor") or entry.name == "dynamicCode":
            continue
        elif NUMBER.match(entry.name):
            if entry.resolve() not in kept:
                continue
            if entry.name.strip("0.") == "":  # the 0 directory: initial fields, keep all
                add_tree(entry, arc)
            else:
                add_tree(entry, arc, skip_file=lambda f: f.parent.name in ("fluid", "solid")
                         and excluded(f.name, drop_fields))
        elif entry.name in SIDE_CASES:
            add_tree(entry, arc, skip_dir=lambda d: d.name.startswith("processor")
                     or d.name == "dynamicCode"
                     or (NUMBER.match(d.name) is not None and d.parent == entry and d.name.strip("0.") != ""))
        else:
            add_tree(entry, arc)
    if meshing is not None and meshing.is_dir():
        for f in sorted(meshing.iterdir()):
            if f.is_file() and not f.is_symlink():
                members.append((f, f"{name}/meshing/{f.name}"))
    return members


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_archive(members: list[tuple[Path, str]], target: Path, level: int):
    """tar the members through zstd into target. The archive is written to a
    temporary name and checked there, so an existing target is only replaced
    by an archive that has passed verify()."""
    tmp = target.with_name(target.name + ".partial")
    file_list = target.with_name(target.name + ".files")
    # GNU tar reads (source -> archive name) pairs only through --transform,
    # so stage the names with a symlink farm instead: simpler and exact.
    stage = target.with_name(target.name + ".stage")
    if stage.exists():
        shutil.rmtree(stage)
    try:
        for src, arc in members:
            link = stage / arc
            link.parent.mkdir(parents=True, exist_ok=True)
            link.symlink_to(src.resolve())
        file_list.write_text("\n".join(arc for _, arc in members) + "\n", encoding="utf-8")
        # --hard-dereference: a symlink and its target are the same file, which
        # tar would otherwise store once plus a hard link; store both in full.
        command = ["tar", "--create", "--dereference", "--hard-dereference", "--file", str(tmp),
                   "--use-compress-program", f"zstd -T0 -{level}",
                   "-C", str(stage), "--verbatim-files-from", "--files-from", str(file_list)]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"tar failed: {result.stderr.strip()}")
        verify(tmp, members)
        tmp.replace(target)
    finally:
        shutil.rmtree(stage, ignore_errors=True)
        file_list.unlink(missing_ok=True)
        tmp.unlink(missing_ok=True)


def verify(target: Path, members: list[tuple[Path, str]]):
    """List the archive back through zstd and compare names and sizes."""
    listing = subprocess.run(["tar", "--list", "--verbose", "--use-compress-program", "zstd -d",
                              "--file", str(target), "--quoting-style=literal"],
                             capture_output=True, text=True)
    if listing.returncode != 0:
        raise RuntimeError(f"cannot read back {target}: {listing.stderr.strip()}")
    sizes = {}
    for line in listing.stdout.splitlines():
        # -rw-r--r-- user/group SIZE YYYY-MM-DD HH:MM name
        parts = line.split(None, 5)
        if len(parts) == 6 and parts[0].startswith("-"):
            sizes[parts[5]] = int(parts[2])
    problems = []
    for src, arc in members:
        expected = src.stat().st_size
        if arc not in sizes:
            problems.append(f"missing: {arc}")
        elif sizes[arc] != expected:
            problems.append(f"size differs: {arc} ({sizes[arc]} != {expected})")
    if len(sizes) != len(members):
        problems.append(f"{len(sizes)} files in the archive, {len(members)} expected")
    if problems:
        raise RuntimeError("archive check failed:\n    " + "\n    ".join(problems[:20]))


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def archive(case: Path, args) -> None:
    case = case.resolve()
    name = case.name
    if not (case / "system/controlDict").is_file():
        raise ValueError(f"{case} is not a run case (no system/controlDict)")
    jobs = running_jobs(case)
    if jobs:
        raise ValueError(f"Slurm job(s) {', '.join(jobs)} are running in {case}")
    # A job submitted from another directory is not matched above, so also
    # refuse while the solver log is still being written.
    recent = [f for f in (case / "log.solids4Foam", case / "log.Allrun")
              if f.is_file() and dt.datetime.now().timestamp() - f.stat().st_mtime < args.quiet_minutes * 60]
    if recent and not args.force_active:
        raise ValueError(f"{recent[0].name} changed in the last {args.quiet_minutes:g} min: is the run still "
                         "going? (--force-active to archive anyway)")

    period = args.period or cycle_period(case)
    delta_t = float(foam_entry(case / "system/controlDict", "deltaT") or 1e-6)
    end_time = float(foam_entry(case / "system/controlDict", "endTime") or 0)
    tol = 0.5 * delta_t
    root_times = time_dirs(case)
    proc_times = time_dirs(case / "processor0")
    last = max([*root_times, *proc_times, 0.0])
    if last < end_time - tol and not args.allow_incomplete:
        raise ValueError(f"the run stopped at t = {last:g} s, before endTime {end_time:g} s "
                         "(--allow-incomplete to archive it anyway)")
    start = last - period
    wanted = sorted(t for t in {*root_times, *proc_times} if start - tol <= t <= last + tol)
    if not wanted:
        raise ValueError("no written times in the last cycle")
    if wanted[0] > start + tol:
        print(f"  warning: the first write of the last cycle is t = {wanted[0]:g}, not {start:g} "
              "(written before purgeWrite was set up?)", file=sys.stderr)

    results = args.results_dir / name / "metadata.json"
    if not args.skip_results_check:
        if not results.is_file():
            raise ValueError(f"no extracted results at {results.parent}; run extract_results.py first "
                             "(or --skip-results-check)")
        meta = json.loads(results.read_text(encoding="utf-8"))
        if Path(meta.get("runCase", "")).resolve() != case:
            raise ValueError(f"{results} was extracted from {meta.get('runCase')}, not {case}")
        extracted_last = meta.get("lastWriteTime_s", meta.get("lastTime_s"))
        if extracted_last is None or abs(extracted_last - last) > max(tol, 1e-9):
            raise ValueError(f"{results} was extracted at t = {extracted_last}, but the run is now at "
                             f"t = {last:g}: extract it again (--force)")

    def reconstructed(t):
        d = next((p for tt, p in root_times.items() if abs(tt - t) <= tol), None)
        return d is not None and all((d / r).is_dir() for r in ("fluid", "solid"))

    missing = [t for t in wanted if t > tol and not reconstructed(t)]
    if missing and args.dry_run:
        print(f"  dry run: would reconstructPar t = {min(missing):g} to {max(missing):g}; the size below "
              "leaves those times out")
    elif missing:
        if not proc_times:
            raise ValueError(f"times {missing} are neither reconstructed nor decomposed")
        reconstruct(case, min(missing), max(missing))
        root_times = time_dirs(case)
        still = [t for t in wanted if t > tol and not reconstructed(t)]
        if still:
            raise RuntimeError(f"times {still} still not reconstructed")
    keep = [p for t, p in root_times.items() if any(abs(t - w) <= tol for w in wanted) or t <= tol]

    drop = [] if args.keep_all_fields else RESTART_ONLY
    meshing = meshing_dir(case)
    if meshing is None or not meshing.is_dir():
        print(f"  warning: meshing directory {meshing} not found; not archived", file=sys.stderr)
    members = collect(case, keep, drop, meshing)
    total = sum(src.stat().st_size for src, _ in members)

    target = args.archive_dir / f"{name}.tar.zst"
    manifest_path = args.archive_dir / f"{name}.manifest.json"
    if target.exists() and not args.force:
        raise FileExistsError(f"{target} exists (--force to replace it)")
    if target.exists() and manifest_path.is_file():
        source = json.loads(manifest_path.read_text(encoding="utf-8")).get("source")
        if source and Path(source).resolve() != case:
            raise FileExistsError(f"{target} holds a different case ({source}); not replacing it")
    existing = args.archive_dir
    while not existing.exists():
        existing = existing.parent
    free = shutil.disk_usage(existing).free
    kept_times = sorted(float(p.name) for p in keep)
    print(f"  {len(members)} files, {human(total)} before compression; "
          f"times {', '.join(f'{t:g}' for t in kept_times[:3])}"
          f"{' ... ' + f'{kept_times[-1]:g}' if len(kept_times) > 3 else ''} ({len(kept_times)}); "
          f"{human(free)} free on {args.archive_dir}")
    if total > free:
        raise OSError(f"not enough space on {args.archive_dir}")
    if args.dry_run:
        print("  dry run: nothing written")
        return

    args.archive_dir.mkdir(parents=True, exist_ok=True)
    print(f"  writing {target}")
    write_archive(members, target, args.level)
    size = target.stat().st_size
    digest = sha256(target)
    manifest = {
        "case": name,
        "source": str(case),
        "archive": str(target),
        "created": dt.datetime.now().isoformat(timespec="seconds"),
        "sha256": digest,
        "compressedBytes": size,
        "uncompressedBytes": total,
        "period_s": period,
        "lastTime_s": last,
        "cycle_s": [start, last],
        "keptTimes_s": kept_times,  # includes the 0 directory
        "droppedFields": drop,
        "meshingDirectory": str(meshing) if meshing else None,
        "files": [{"name": arc, "bytes": src.stat().st_size} for src, arc in members],
    }
    manifest_path.write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")
    print(f"  checked: {len(members)} files, {human(size)} compressed ({size / total:.0%}); "
          f"manifest {manifest_path.name}")

    meta = json.loads(results.read_text(encoding="utf-8")) if results.is_file() else {}
    if Path(meta.get("runCase", "/nonexistent")).resolve() == case:
        meta["archive"] = {"path": str(target), "sha256": digest, "created": manifest["created"],
                           "keptTimes_s": [kept_times[0], kept_times[-1]], "count": len(kept_times)}
        results.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")

    if args.keep:
        print(f"  kept {case} (--keep)")
        return
    if not args.delete:
        if not sys.stdin.isatty():
            print(f"  kept {case}: pass --delete to remove it when not run interactively")
            return
        answer = input(f"  Delete {case}? [y/N] ").strip().lower()
        if answer not in ("y", "yes"):
            print(f"  kept {case}")
            return
    shutil.rmtree(case)
    print(f"  deleted {case}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=__doc__.split("\n\n", 1)[1])
    parser.add_argument("cases", nargs="+", type=Path, help="run case directories")
    parser.add_argument("--archive-dir", type=Path, default=DEFAULT_ARCHIVE,
                        help=f"where <case>.tar.zst goes (default {DEFAULT_ARCHIVE})")
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS,
                        help=f"extracted results (default {DEFAULT_RESULTS})")
    parser.add_argument("--period", type=float, help="cycle length [s] (default: from the inlet waveform)")
    parser.add_argument("--level", type=int, default=3, help="zstd compression level (default 3)")
    parser.add_argument("--keep-all-fields", action="store_true",
                        help=f"also archive the restart-only fields ({', '.join(RESTART_ONLY)})")
    parser.add_argument("--allow-incomplete", action="store_true", help="archive a run that stopped before endTime")
    parser.add_argument("--skip-results-check", action="store_true",
                        help="archive even if results/<case>/ has not been extracted")
    parser.add_argument("--dry-run", action="store_true", help="report what would be archived, write nothing")
    parser.add_argument("--force", action="store_true", help="replace an existing archive")
    delete = parser.add_mutually_exclusive_group()
    delete.add_argument("--delete", action="store_true", help="delete the run case after a checked archive, without asking")
    delete.add_argument("--keep", action="store_true", help="never delete the run case")
    parser.add_argument("--quiet-minutes", type=float, default=15,
                        help="refuse if the solver log changed this recently (default 15)")
    parser.add_argument("--force-active", action="store_true",
                        help="archive even if the solver log changed recently")
    args = parser.parse_args()
    args.archive_dir = args.archive_dir.resolve()
    args.results_dir = args.results_dir.resolve()
    failures = 0
    for case in args.cases:
        print(f"Archiving {case}")
        try:
            archive(case, args)
        except (OSError, ValueError, RuntimeError) as error:
            print(f"  error: {error}", file=sys.stderr)
            failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
