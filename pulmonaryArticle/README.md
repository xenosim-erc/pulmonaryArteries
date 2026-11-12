# pulmonaryArticle

Draft of the human vs porcine pulmonary-artery FSI article.

- `main.tex` — top level; one file per section in `sections/`
- `macros.tex` — shared notation and `\todo{}`
- `figures/` — generated figures only (PDF); never edit by hand
- `scripts/` — post-processing scripts that produce `figures/` and `tables/`
  from the simulation results
- `tables/` — generated LaTeX tables, `\input` from the sections
- `data/` — small extracted results (CSV) that the scripts read

Build with `latexmk` (output in `build/`).
