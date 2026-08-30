#!/usr/bin/env python3
"""Insert a zero-face patch as the first entry in an OpenFOAM boundary file."""

from __future__ import annotations

import argparse
import os
import re
import stat
import sys
import tempfile
from pathlib import Path


BOUNDARY_LIST = re.compile(
    r"(?m)^(?P<indent>[ \t]*)(?P<count>[0-9]+)[ \t]*\r?\n"
    r"(?P<paren_indent>[ \t]*)\([ \t]*\r?\n"
)


def foam_word(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", value):
        raise argparse.ArgumentTypeError(f"invalid OpenFOAM word: {value!r}")
    return value


def add_patch(boundary_file: Path, patch_name: str, patch_type: str) -> None:
    text = boundary_file.read_text(encoding="utf-8")
    match = BOUNDARY_LIST.search(text)
    if match is None:
        raise RuntimeError("could not find the boundary patch count and list")

    # Refuse to create a duplicate patch. Patch names appear alone immediately
    # before their opening brace in a standard polyMesh/boundary file.
    existing_name = re.compile(
        rf"(?m)^[ \t]*{re.escape(patch_name)}[ \t]*\r?\n[ \t]*\{{"
    )
    if existing_name.search(text, match.end()):
        raise RuntimeError(f"patch {patch_name!r} already exists")

    old_count = int(match.group("count"))
    first_start_match = re.search(
        r"(?m)^[ \t]*startFace[ \t]+(?P<start>[0-9]+)[ \t]*;",
        text[match.end() :],
    )
    if first_start_match is None:
        raise RuntimeError("could not determine the first boundary-face index")
    first_start_face = int(first_start_match.group("start"))
    prefix = (
        f"{match.group('indent')}{old_count + 1}\n"
        f"{match.group('paren_indent')}(\n"
        f"    {patch_name}\n"
        "    {\n"
        f"        type            {patch_type};\n"
        "        nFaces          0;\n"
        f"        startFace       {first_start_face};\n"
        "    }\n"
    )
    updated = text[: match.start()] + prefix + text[match.end() :]

    original_mode = stat.S_IMODE(boundary_file.stat().st_mode)
    temporary_name = ""
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=boundary_file.parent,
            prefix=f".{boundary_file.name}.",
            delete=False,
        ) as temporary:
            temporary.write(updated)
            temporary.flush()
            os.fsync(temporary.fileno())
            temporary_name = temporary.name
        os.chmod(temporary_name, original_mode)
        os.replace(temporary_name, boundary_file)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("boundary_file", type=Path)
    parser.add_argument("patch_name", type=foam_word)
    parser.add_argument("patch_type", type=foam_word)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    boundary_file = args.boundary_file.expanduser().resolve()
    if not boundary_file.is_file():
        raise FileNotFoundError(boundary_file)
    add_patch(boundary_file, args.patch_name, args.patch_type)
    print(
        f"Inserted zero-face patch {args.patch_name!r} "
        f"of type {args.patch_type!r} into {boundary_file}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, OSError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
