// Tests of the large-table cell parsing and ordering,
// docs/src/.vitepress/theme/table-sort.mjs.
// Run with `node --test docs/test/table_sort.test.mjs` (no packages needed).
import test from 'node:test'
import assert from 'node:assert/strict'
import {
  cellNumber,
  compareCells,
  isBlank,
  isNumericColumn,
  NUMERIC_SHARE,
} from '../src/.vitepress/theme/table-sort.mjs'

const sortColumn = (texts, numeric, descending) =>
  texts.slice().sort((x, y) => compareCells(x, y, numeric, descending))

test('the first number of a cell', () => {
  assert.equal(cellNumber('1.2 (0.8–1.9)'), 1.2)
  assert.equal(cellNumber('0.73–1.53'), 0.73)
  assert.equal(cellNumber('4-5'), 4)
  assert.equal(cellNumber('12%'), 12)
  assert.equal(cellNumber('-2.21'), -2.21)
  assert.equal(cellNumber('.5'), 0.5)
  assert.equal(cellNumber('< 0.01'), 0.01)
})

test('the U+2212 minus sign', () => {
  assert.equal(cellNumber('−2.21'), -2.21)
  assert.equal(cellNumber('−0.5 (−1.0–0.1)'), -0.5)
  assert.equal(cellNumber('2e−3'), 0.002)
})

test('thousands separators', () => {
  assert.equal(cellNumber('286,637'), 286637)
  assert.equal(cellNumber('1,234,567.5'), 1234567.5)
  assert.equal(cellNumber('331 (156–680)'), 331)
  // A comma not followed by three digits is not a separator.
  assert.equal(cellNumber('1,2'), 1)
})

test('exponents', () => {
  assert.equal(cellNumber('1e-3'), 0.001)
  assert.equal(cellNumber('2.5E+4'), 25000)
})

test('ISO dates and text are not numbers', () => {
  assert.ok(Number.isNaN(cellNumber('2026-10-02')))
  assert.ok(Number.isNaN(cellNumber('2026-10-02T12:00')))
  assert.ok(Number.isNaN(cellNumber('Mandima')))
  assert.ok(Number.isNaN(cellNumber('P(R > 1)')))
  assert.ok(Number.isNaN(cellNumber('Week 12')))
  assert.ok(Number.isNaN(cellNumber('')))
})

test('blank cells', () => {
  for (const t of ['', '  ', '-', '—', '–', 'missing', 'NaN', 'n/a']) assert.ok(isBlank(t), t)
  for (const t of ['0', 'no', '1.2', 'Beni']) assert.ok(!isBlank(t), t)
})

test('numeric sort uses the first number, not the text', () => {
  const cells = ['9.1 (2.0–20.0)', '50.2 (17.8–119.7)', '12%', '-2.21', '−3', '286,637', '1e-3']
  assert.deepEqual(sortColumn(cells, true, false), [
    '−3', '-2.21', '1e-3', '9.1 (2.0–20.0)', '12%', '50.2 (17.8–119.7)', '286,637',
  ])
  assert.deepEqual(sortColumn(cells, true, true), [
    '286,637', '50.2 (17.8–119.7)', '12%', '9.1 (2.0–20.0)', '1e-3', '-2.21', '−3',
  ])
})

test('blank and missing cells sort last in both directions', () => {
  const cells = ['missing', '3', '', '1', '-', '2', 'NaN']
  for (const descending of [false, true]) {
    const out = sortColumn(cells, true, descending)
    assert.deepEqual(out.slice(0, 3), descending ? ['3', '2', '1'] : ['1', '2', '3'])
    assert.ok(out.slice(3).every(isBlank), out.join(','))
  }
  const text = sortColumn(['Beni', 'missing', 'Aru', ''], false, true)
  assert.deepEqual(text.slice(0, 2), ['Beni', 'Aru'])
})

test('text cells in a numeric column sort after the numbers', () => {
  for (const descending of [false, true]) {
    const out = sortColumn(['n.a. yet', '2', 'pending', '10'], true, descending)
    assert.deepEqual(out.slice(0, 2), descending ? ['10', '2'] : ['2', '10'])
  }
})

test('text and date columns sort as text', () => {
  assert.deepEqual(sortColumn(['zone 10', 'Zone 2', 'adja'], false, false), ['adja', 'Zone 2', 'zone 10'])
  assert.deepEqual(
    sortColumn(['2026-10-10', '2026-09-27', '2026-10-02'], false, false),
    ['2026-09-27', '2026-10-02', '2026-10-10']
  )
})

test('a column is numeric when 80% of its filled cells are numbers', () => {
  assert.equal(NUMERIC_SHARE, 0.8)
  const col = (n, t) => [...Array(n).fill('1.5'), ...Array(t).fill('yes')]
  assert.ok(isNumericColumn(col(8, 2)))
  assert.ok(!isNumericColumn(col(7, 3)))
  assert.ok(isNumericColumn([...col(4, 1), '', 'missing', '-']), 'blanks do not count')
  assert.ok(!isNumericColumn(['', 'missing']), 'an all-blank column is text')
  assert.ok(!isNumericColumn(['2026-10-02', '2026-10-06']), 'dates are text')
})
