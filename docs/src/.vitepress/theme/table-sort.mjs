// Cell parsing and ordering for the large-table enhancement in tables.ts.
// Plain JavaScript with no DOM access, so docs/test/table_sort.test.mjs can
// run it under `node --test` without a build step or packages.

// The leading number of a cell: "1.2 (0.8–1.9)", "12%", "-2.21", "−2.21"
// (U+2212), "286,637" and "1e-3" all qualify.
const NUMBER = /^[^\d\p{L}+\-−.]*([+\-−]?(?:\d{1,3}(?:,\d{3})+|\d+)?(?:\.\d+)?(?:e[+\-−]?\d+)?)/iu
const DATE = /^\s*\d{4}-\d{2}-\d{2}/
const BLANK = /^\s*(|-|—|–|missing|NaN|n\/a)\s*$/i

// Share of a column's filled cells that must parse as numbers for the
// column to sort and align as numbers.
export const NUMERIC_SHARE = 0.8

/**
 * The first number in a cell, or NaN when it does not start with one.
 * ISO dates are NaN so that they sort as text.
 * @param {string} text
 * @returns {number}
 */
export function cellNumber(text) {
  if (DATE.test(text)) return NaN
  const m = NUMBER.exec(text)
  if (!m || !/\d/.test(m[1])) return NaN
  return Number(m[1].replace(/,/g, '').replace(/−/g, '-'))
}

/**
 * Whether a cell holds no value: empty, a dash, "missing", "NaN" or "n/a".
 * @param {string} text
 * @returns {boolean}
 */
export function isBlank(text) {
  return BLANK.test(text)
}

/**
 * Whether a column's cells are numbers, ignoring blank cells.
 * @param {string[]} texts
 * @returns {boolean}
 */
export function isNumericColumn(texts) {
  let filled = 0
  let numeric = 0
  for (const text of texts) {
    if (isBlank(text)) continue
    filled++
    if (!Number.isNaN(cellNumber(text))) numeric++
  }
  return filled > 0 && numeric / filled >= NUMERIC_SHARE
}

const COLLATOR = new Intl.Collator(undefined, { numeric: true, sensitivity: 'base' })

/**
 * Sort order of two cells in one column. Blank cells come last in both
 * directions, and in a numeric column so do cells without a number.
 * @param {string} x
 * @param {string} y
 * @param {boolean} numeric whether the column is numeric
 * @param {boolean} descending
 * @returns {number}
 */
export function compareCells(x, y, numeric, descending) {
  const bx = isBlank(x)
  const by = isBlank(y)
  if (bx || by) return bx === by ? 0 : bx ? 1 : -1
  const sign = descending ? -1 : 1
  if (numeric) {
    const nx = cellNumber(x)
    const ny = cellNumber(y)
    const ix = Number.isNaN(nx)
    const iy = Number.isNaN(ny)
    if (ix && iy) return sign * COLLATOR.compare(x, y)
    if (ix || iy) return ix ? 1 : -1
    return sign * (nx - ny)
  }
  return sign * COLLATOR.compare(x, y)
}
