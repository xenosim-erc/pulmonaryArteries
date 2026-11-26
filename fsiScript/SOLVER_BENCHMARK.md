# FSI solver-cost benchmark (November 2026)

Why the Robin-Neumann FSI cases had become slow, what was tested, and what the
templates now use. The benchmark cases themselves have been deleted; this file
keeps their results.

## Practical upshots (now in `template/templateCaseRobin`)

- **The solid solve is the cost.** Profiling puts all fluid linear solves
  (pressure, velocity, mesh motion) at about 1.5% of the run time and the
  interface interpolation below 0.1%; about 98% is the solid's PETSc solve. A
  finer fluid mesh costs little; the solid mesh size and solver settings set
  the cost.
- **Solid preconditioner: BoomerAMG with tuned options** (HMIS coarsening,
  ext+i interpolation, aggressive coarsening; the same as the Dirichlet
  template) instead of block-Jacobi with LU: 2.7x faster per step.
- **Solid tolerances: `snes_rtol 1e-3`, `snes_atol 1e-6`** (was `snes_rtol
  1e-6`, no `atol`): half the Newton iterations, a further 1.7x, with results
  identical to the tight tolerances. The FSI outer iterations converge the
  coupled step anyway.
- **32 cores** instead of 16: 1.3x on their own, 66% parallel efficiency;
  with the looser tolerances, 2.85x over the 16-core tight-tolerance run.
- **Time step: keep 2e-4 s.** 2.5e-4 s saves only 10% per beat (more FSI
  iterations); at 5e-4 s the coupling barely converges (10 FSI iterations per
  step) and a beat costs 3x more.
- **The Fung law is not the problem.** Neo-Hookean at the same small-strain
  stiffness needs exactly the same iterations. Prestress costs about 25% per
  step at equal iterations.

Overall, on this case, 100 s per step (old set-up) became 13 s per step
(current template, 32 cores).

## Case and method

- h08 human case from `meshing.sh` defaults at the time: 366,741 fluid cells,
  179,643 solid cells (3 wall layers; the default is now 5, which makes the
  solid about 1.7x larger), maxCellSize 0.8 mm, both meshes pass `checkMesh`.
- Robin-Neumann coupling, isotropic Fung wall (c 7.87 kPa, a 7.05, K 1.38 MPa)
  with prestress at 10 mmHg, human waveform and Windkessel, time step 2e-4 s
  unless stated, MPI ranks bound to cores on xenosim.
- Each variant: 100 steps from t = 0, field writes off. Rates are over steps
  51 to 100, still in the start-up, so they overstate the cost over a whole
  beat. "Hours per beat" = that rate times the steps in a 0.86 s beat.
- Cost split by OpenFOAM's `profiling` dictionary in `system/controlDict`
  (`active true`), written to `processor*/<time>/uniform/profiling`. The
  solid's PETSc solve is not instrumented, so its share is the remainder.

## Results

All on 16 cores unless stated. Its: iterations per time step.

| Variant | s/step | Hours per beat | FSI its | Solid Newton its | Solid linear its |
| --- | --- | --- | --- | --- | --- |
| Original: block-Jacobi LU, `rtol 1e-6` | 100 | 120 | 3.5 | 8.9 | 1057 |
| Time step 2.5e-4 s | 113 | 108 | 4.3 | 10.5 | 1413 |
| Time step 5e-4 s (capped at 47 steps) | 332 | 158 | 10.4 | 22.8 | 4095 |
| No prestress | 77 | 92 | 3.7 | 9.4 | 1149 |
| Neo-Hookean E 82.5 kPa (Fung small-strain stiffness), no prestress | 88 | 105 | 3.7 | 9.4 | 1148 |
| Plain BoomerAMG | 55 | 66 | 3.5 | 8.9 | 488 |
| Plain BoomerAMG, neo-Hookean E 82.5 kPa, no prestress | 46 | 55 | 3.7 | 9.4 | 515 |
| Plain BoomerAMG, neo-Hookean E 300 kPa, no prestress | 73 | 87 | 5.6 | 13.1 | 885 |
| Tuned BoomerAMG (reference for the last block) | 37 | 45 | 3.5 | 8.9 | 716 |
| Tuned BoomerAMG, 32 cores | 28 | 34 | 3.5 | 8.9 | 729 |
| Tuned BoomerAMG, `rtol 1e-3` + `atol 1e-6` | 22 | 26 | 3.5 | 4.5 | 312 |
| **Tuned BoomerAMG, `rtol 1e-3` + `atol 1e-6`, 32 cores (current template)** | **13** | **16** | 3.5 | 4.5 | 314 |

**Accuracy of the looser tolerances and 32 cores**, against the 16-core
tight-tolerance run at the same times over the 100 steps: inlet pressure within
1.6e-5 mmHg, LPA and RPA flow within 0.0007 mL/s, lumen volume within 3e-6 mL,
largest wall displacement within 0.004 um (of about 740 um).

**Profile of the plain BoomerAMG run** (7126 s): fluid pressure solve 79 s,
mesh motion 9 s, fluid velocity 8 s, interface interpolation 5 s, the rest
(the solid solve) about 98%.

## Comparison with the October case

The October case (`~/OpenFOAM/conormoore-v2412/run/prestressTutorial/
pulmonaryArteriesNoPrestress`) ran at 3.8 s per step over 5 s of flow: 87,459
fluid and 68,781 solid cells, neo-Hookean E 300 kPa, no prestress, time step
2.5e-4 s, block-Jacobi LU. Over the same steps 51 to 100 it needed similar
iterations (4.2 FSI, 10.4 Newton, 476 linear per step, at 7.0 s per step):
the difference is the cost of each solid linear iteration, 15 ms then against
52 ms with tuned BoomerAMG here, which follows the 2.6x larger solid mesh and
the prestress. After the start-up it fell to 2.0 FSI iterations and 3.6 s
per step.

## Wall stiffness

A stiffer wall does not cut the FSI iterations: neo-Hookean at 300 kPa needed
5.6 per step against 3.7 at 82.5 kPa. For reference, the Fung wall's
equibiaxial tangent stiffness against the 300 kPa neo-Hookean wall, at the
prestretch read from the case (largest principal stretch median 1.08, 90th
percentile 1.13):

| Stretch | Fung tangent | Neo-Hookean 300 kPa | Ratio |
| --- | --- | --- | --- |
| 1.00 | 166 kPa | 604 kPa | 0.28 |
| 1.08 | 229 kPa | 531 kPa | 0.43 |
| 1.13 | 447 kPa | 504 kPa | 0.89 |
| 1.16 | 718 kPa | 493 kPa | 1.46 |

The Fung wall reaches the 300 kPa wall's stiffness at a stretch of about 1.14.
