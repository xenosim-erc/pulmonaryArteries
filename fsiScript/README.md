# Pulmonary-artery FSI workflow

This directory prepares an open pulmonary-artery STL for a coupled
fluid-structure interaction (FSI) simulation. The workflow extracts vessel
centerlines, estimates a diameter-dependent wall thickness, creates fluid and
solid meshes, and assembles an OpenFOAM/solids4foam run case.

## What is version-controlled

- `geometries/`: input STL surfaces
- `pythonScripts/`: VTK/VMTK geometry-processing utilities
- `templateMesh/`: reusable cfMesh and wall-extrusion case files
- `templateCase/`: reusable solids4foam case files
- `meshing.sh`: the end-to-end meshing driver
- `environment.yml`: the Python dependency specification

The local Conda installation and environment (`.env/`), intermediate meshing
files (`caseFiles/`), and assembled solver cases (`run/`) are generated locally
and ignored by Git. These directories are large, machine-specific, or
reproducible from the tracked inputs.

## Requirements

The Python tools use Python 3.10, VMTK 1.5, VTK 9.2.6, ITK 5.3, NumPy 1.26,
and SciPy. They are installed from conda-forge using `environment.yml`.

The complete workflow also requires a compatible OpenFOAM installation,
solids4foam, cfMesh (`surfaceFeatureEdges` and `cartesianMesh`), and the local
`mapExtrudeDistance` and `varExtrudeMesh` utilities from this repository. The
included Slurm template currently loads OpenFOAM v2412 and should be adapted to
the target cluster.

## Create the Python environment

Install Miniforge, Mambaforge, or Miniconda outside this repository, then run
from `fsiScript`:

```bash
conda env create --prefix .env/fsi-env --file environment.yml
conda activate ./.env/fsi-env
```

The distribution's own version is not a project dependency: Miniforge and
Miniconda are installers for the Conda package manager. Reproducibility comes
from the dependencies in `environment.yml`, not from committing an installer
or a complete environment directory. Compatibility-sensitive packages are
pinned here; Conda resolves their platform-specific transitive dependencies.

To update the environment after changing the file:

```bash
conda env update --prefix .env/fsi-env --file environment.yml --prune
```

For a bit-for-bit platform lock, generate a separate lock file with a tool such
as `conda-lock`; keep `environment.yml` as the readable source specification.

## Run the workflow

First source the required OpenFOAM/solids4foam environment so their commands
are on `PATH`. Then, from this directory, run for example:

```bash
./meshing.sh geometries/h08.stl
```

The script uses `.env/fsi-env/bin/python` by default. An environment installed
elsewhere can be selected explicitly:

```bash
FSI_PYTHON=/path/to/fsi-env/bin/python ./meshing.sh geometries/h08.stl
```

Use `./meshing.sh` without arguments to see the command-line options. Generated
intermediate geometry is written beneath `caseFiles/`, and the assembled FSI
case is written beneath `run/`.

Before running a generated case, review the material properties, boundary
conditions, decomposition settings, processor count, and cluster-specific
paths in `templateCase/`.
