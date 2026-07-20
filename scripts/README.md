# Scripts

Helper scripts for the analysis.
Most are Julia and run against the `scripts/` project:

```sh
julia --project=scripts scripts/<name>.jl
```

`scripts/download_sitreps.jl` fetches the INSP SitRep PDFs into
`data/sitrep_pdfs/` (see `data/README.md`).

## `digitize_onset_curve` — onset epidemic-curve digitiser

Digitises the symptom-onset epidemic-curve figure out of the analytique
SitRep PDFs into `data/onset_curve_scanned.csv` (the figure is a raster bar
chart with no data table; see `data/README.md`).

There are two implementations that produce a byte-identical CSV:

- `digitize_onset_curve.jl` — the **Julia reference**, dependency-free (it
  parses the PPM that poppler's `pdfimages` emits by default, so it needs no
  image library). This is the in-repo version:

  ```sh
  julia scripts/digitize_onset_curve.jl
  # optional: julia scripts/digitize_onset_curve.jl <pdf_dir> <out_csv>
  ```

- `digitize_onset_curve.py` — a Python port kept for the automated
  data-updater, which has Python (not Julia) access. Same algorithm, same
  output.

Poppler (`pdfimages` / `pdftotext` / `pdfinfo`) must be on `PATH` for both
(`apt install poppler-utils`, or `brew install poppler`). Download the SitRep
PDFs first with `scripts/download_sitreps.jl`.

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
