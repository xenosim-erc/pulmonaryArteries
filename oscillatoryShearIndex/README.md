# oscillatoryShearIndex

OpenFOAM function object that computes the oscillatory shear index (OSI) on
wall patches. The FSI case templates in `fsiScript/template/` load it with
`libs (oscillatoryShearIndex);` and `fsiScript/pythonScripts/extract_results.py`
uses its output for TAWSS and OSI.

Taken from <https://github.com/solids4foam/oscillatoryShearIndex> (commit
57cd043; author Philip Cardiff, philip.cardiff@ucd.ie), GPL v3 or later like
the rest of this repository.

## What it computes

At every time step, on the selected patches, it takes the wall shear stress
WSS = -n · devReff (kinematic, per unit density, for an incompressible solver)
and accumulates

- `<name>_accumulatedWallShearStress` = ∫ WSS dt (vector)
- `<name>_accumulatedMagWallShearStress` = ∫ |WSS| dt

then sets OSI = 0.5 (1 - |∫ WSS dt| / ∫ |WSS| dt). With `writeFields yes`
all three fields are written at write times; the minimum and maximum OSI per
patch go to `postProcessing/<name>/<startTime>/oscillatoryShearIndex.dat`.
`<name>` is `oscillatoryShearIndex` (the function object's type name).

TAWSS and OSI over any interval [t0, t1] whose two ends were written follow
from the differences of the two integrals: TAWSS = Δ∫|WSS| dt / (t1 - t0),
mean WSS = Δ∫WSS dt / (t1 - t0).

## Change from upstream

Upstream creates the two integrals with `NO_READ`, so every solver restart
starts them again from zero and the integrals then cover only the time since
the restart. Here they are `READ_IF_PRESENT`: on a restart they are read back
from the start-time directory and carry on; a run from t = 0 has no such files
and starts from zero as before. Uninterrupted runs give identical results.

Checked on upstream's cavity case: running to t = 0.3 and restarting to 0.5,
the |WSS| integral at 0.5 is 0.39 of the uninterrupted run's with the upstream
code (only 0.3-0.5 s counted) and equal to it with this change.

To start a deliberately fresh count from a non-zero time, delete the two
`*_accumulated*` files from that time directory first.

## Building

`./Allwmake` at the repository root builds it, or

    cd oscillatoryShearIndex && wmake libso

with OpenFOAM sourced. It installs `liboscillatoryShearIndex.so` into
`$FOAM_USER_LIBBIN`, replacing any upstream build of the same name there.
Tested with OpenFOAM v2412.

## Usage

```
oscillatoryShearIndex1
{
    type            oscillatoryShearIndex;
    libs            (oscillatoryShearIndex);
    region          fluid;          // multi-region (solids4foam) cases
    patches         (wall);
    writeFields     yes;
    writeControl    writeTime;
}
```
