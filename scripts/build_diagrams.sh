#!/usr/bin/env bash
# Compile every TikZ diagram in docs/diagrams/ (except the shared style.tex)
# with pdflatex and write docs/src/public/diagrams/<name>.svg.
# Needs pdflatex (TeX Live with the pgf, standalone and helvet packages) and
# poppler's pdftocairo on PATH. Run from anywhere:
#
#   scripts/build_diagrams.sh          # every diagram
#   scripts/build_diagrams.sh --png    # also write 150 dpi PNG previews to
#                                      # output/diagrams/ or $BVD_DIAGRAM_PNG_DIR
#
# The SVGs are committed, since CI has no LaTeX.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
src_dir="$repo_root/docs/diagrams"
out_dir="$repo_root/docs/src/public/diagrams"
png=false
[[ "${1:-}" == "--png" ]] && png=true

for tool in pdflatex pdftocairo; do
  command -v "$tool" >/dev/null || {
    echo "build_diagrams: $tool is not on PATH" >&2
    exit 1
  }
done

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$out_dir"

shopt -s nullglob
for tex in "$src_dir"/*.tex; do
  name="$(basename "$tex" .tex)"
  [[ "$name" == "style" ]] && continue
  echo "build_diagrams: $name"
  (
    cd "$src_dir"
    pdflatex -interaction=nonstopmode -halt-on-error \
      -output-directory "$tmp" "$tex" >"$tmp/$name.out" 2>&1
  ) || {
    tail -n 40 "$tmp/$name.out" >&2
    echo "build_diagrams: $name failed to compile" >&2
    exit 1
  }
  pdftocairo -svg "$tmp/$name.pdf" "$out_dir/$name.svg"
  if $png; then
    pdftocairo -png -r 150 -singlefile "$tmp/$name.pdf" "$tmp/$name"
    keep="${BVD_DIAGRAM_PNG_DIR:-$repo_root/output/diagrams}"
    mkdir -p "$keep"
    cp "$tmp/$name.png" "$keep/$name.png"
  fi
done
echo "build_diagrams: wrote $out_dir"
