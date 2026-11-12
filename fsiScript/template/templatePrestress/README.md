# Prestress case

This solid-only case computes the prestress `Fpre` of the imaged, loaded
arterial wall, and copies it to `../0/solid/Fpre` for the fluid-solid
interaction run. `meshing.sh` places it in every run case as `prestress/`, and
the run case's `Allrun` runs it before the fluid-solid interaction run. It is
adapted from the prestress case of the solids4foam `pulmonaryArteries`
tutorial.

The imaged geometry is loaded by the blood pressure at the time of imaging, so
it is not stress-free. `Allrun` applies that pressure to `innerWall` and
repeats a modified updated Lagrangian (MULF) pass, `solids4Foam` then
`updatePrestress`, until the imaged wall is in equilibrium with it. The
deformation gradient the wall law sees in the fluid-solid interaction run is
then `F & Fpre`.

Set by `meshing.sh`:

- the pressure, in `constant/timeVsPressure.ramp` (pass 1) and
  `constant/timeVsPressure.full` (later passes), in Pa;
- `tolerance` in `Allrun`: passes stop once a run's maximum displacement,
  `maxMag(D)` in `prestressHistory`, is below it, in metres;
- `maxPasses` in `Allrun`: the cap on passes;
- `numberOfSubdomains` in `system/decomposeParDict`.

The wall material and the solid mesh are not stored here: `Allrun` copies them
from `../constant/solid/` at the start of every run, so `Fpre` is always
computed for what the fluid-solid interaction run solves. Each pass appends
its convergence measures to `prestressHistory`. If the tolerance is not met,
the latest `Fpre` is still copied, and `PRESTRESS_NOT_CONVERGED` is written
here.

To rerun it on its own from the run case:

```bash
cd prestress
./Allclean
./Allrun            # parallel; ./Allrun serial to run in serial
```
