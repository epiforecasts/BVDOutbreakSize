# Scripts

Helper scripts for the analysis.
Most are Julia and run against the `scripts/` project:

```sh
julia --project=scripts scripts/<name>.jl
```

`scripts/download_sitreps.jl` fetches the INSP SitRep PDFs into
`data/sitrep_pdfs/` (see `data/README.md`).

## `digitize_onset_curve.py` (Python)

This one step is Python, not Julia: it digitises the symptom-onset
epidemic-curve figure out of the analytique SitRep PDFs into
`data/onset_curve_scanned.csv` (image analysis, so Pillow + numpy).
It is kept separate because the automated data-updater that runs it has
Python, not Julia, access.
The digitised CSV is validated on the Julia side by
`test/test_onset_curve.jl`, so the model repo still gets a Julia check.

### Getting the Python dependencies with uv

The script carries [PEP 723](https://peps.python.org/pep-0723/) inline
metadata, so [uv](https://docs.astral.sh/uv/) fetches Pillow and numpy into a
throwaway environment automatically:

```sh
uv run scripts/digitize_onset_curve.py
```

Pass a PDF directory and output path if they differ from the defaults
(`data/sitrep_pdfs` and `data/onset_curve_scanned.csv`):

```sh
uv run scripts/digitize_onset_curve.py data/sitrep_pdfs data/onset_curve_scanned.csv
```

Without uv, install into an isolated venv instead:

```sh
python3 -m venv .venv
.venv/bin/pip install Pillow numpy
.venv/bin/python scripts/digitize_onset_curve.py
```

Either way, poppler must be on `PATH` for `pdfimages` / `pdftotext` /
`pdfinfo` (`apt install poppler-utils`, or `brew install poppler`).
Download the SitRep PDFs first with `scripts/download_sitreps.jl`.
