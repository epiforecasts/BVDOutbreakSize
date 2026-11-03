// Large report tables get a compact style, a scroll box with a sticky
// header, a row filter and sortable columns. Tables with no more than
// `MIN_ROWS` body rows are left alone, as is every table without
// JavaScript. The enhancement only adds a wrapper and a toolbar around the
// table and reorders or hides its rows, so the server-rendered markup is
// unchanged until the page has hydrated. Cell parsing and ordering live in
// table-sort.mjs, which docs/test/table_sort.test.mjs tests.

import { compareCells, isNumericColumn } from './table-sort.mjs'

// About a screen of rows at the compact size. The province tables (eight
// provinces) and the ten-row summaries stay as they are.
const MIN_ROWS = 10

const DONE = 'data-bvd-table'

interface Parts {
  header: HTMLTableRowElement
  rows: HTMLTableRowElement[]
}

// The header row whose cells label the columns, and the body rows; <tfoot>
// rows are neither. A raw HTML table without a <thead> has its leading
// all-<th> row as the header. A <thead> with more than one row is skipped,
// except a DataFrame display's second row, told apart by every cell after the
// row-number stub carrying its element type as its title.
function tableParts(table: HTMLTableElement): Parts | null {
  const rows: HTMLTableRowElement[] = []
  for (const body of Array.from(table.tBodies)) {
    rows.push(...Array.from(body.rows))
  }
  const head = table.tHead
  if (head && head.rows.length > 1 && !isDataFrameHead(table, head)) return null
  const header = head ? head.rows[0] : rows.shift()
  if (!header || header.querySelector('td')) return null
  return { header, rows }
}

function isDataFrameHead(table: HTMLTableElement, head: HTMLTableSectionElement): boolean {
  if (!table.classList.contains('data-frame') || head.rows.length !== 2) return false
  const types = Array.from(head.rows[1].cells).slice(1)
  return types.length > 0 && types.every((c) => c.title !== '' && c.title === c.textContent?.trim())
}

function hasSpans(table: HTMLTableElement): boolean {
  return table.querySelector('[colspan]:not([colspan="1"]), [rowspan]:not([rowspan="1"])') !== null
}

function enhance(table: HTMLTableElement, index: number): void {
  table.setAttribute(DONE, '')
  if (table.querySelector('table')) return
  const parts = tableParts(table)
  if (!parts || parts.rows.length <= MIN_ROWS) return
  const { header, rows } = parts
  // Moved into a <thead> so the header can stick.
  if (!table.tHead) table.createTHead().appendChild(header)
  const sortable = !hasSpans(table) && rows.length > 1
  const ncol = header.cells.length

  const numeric: boolean[] = []
  for (let c = 0; c < ncol; c++) {
    numeric.push(isNumericColumn(rows.map((r) => r.cells[c]?.textContent ?? '')))
    const cls = numeric[c] ? 'bvd-num' : 'bvd-text'
    const foot = table.tFoot ? Array.from(table.tFoot.rows) : []
    for (const row of [...Array.from(table.tHead!.rows), ...rows, ...foot]) {
      row.cells[c]?.classList.add(cls)
    }
  }

  const wrap = document.createElement('div')
  wrap.className = 'bvd-table'
  const tools = document.createElement('div')
  tools.className = 'bvd-table-tools'
  const input = document.createElement('input')
  input.type = 'search'
  input.className = 'bvd-table-filter'
  input.placeholder = 'Filter rows'
  input.setAttribute('aria-label', 'Filter table rows')
  input.autocomplete = 'off'
  const count = document.createElement('span')
  count.className = 'bvd-table-count'
  count.setAttribute('aria-live', 'polite')
  tools.append(input, count)
  const scroll = document.createElement('div')
  scroll.className = 'bvd-table-scroll'
  scroll.tabIndex = 0
  scroll.id = `bvd-table-${index}`
  input.setAttribute('aria-controls', scroll.id)

  table.parentNode!.insertBefore(wrap, table)
  wrap.append(tools, scroll)
  scroll.appendChild(table)
  // A DataFrame display's own wrapper always scrolls sideways; the box
  // scrolls instead, so that wrapper must not add a second scrollbar.
  if (wrap.parentElement?.matches('div.data-frame')) wrap.parentElement.style.overflowX = 'visible'
  table.classList.add('bvd-table-enhanced')
  // Table text, kept lower-cased once so filtering does not reread the DOM.
  const text = new Map(rows.map((r) => [r, (r.textContent ?? '').toLowerCase()]))

  const total = rows.length
  const showCount = (shown: number) => {
    count.textContent = shown === total ? `${total} rows` : `${shown} of ${total} rows`
  }
  showCount(total)

  input.addEventListener('input', () => {
    const terms = input.value.toLowerCase().split(/\s+/).filter(Boolean)
    let shown = 0
    for (const row of rows) {
      const hit = terms.every((t) => text.get(row)!.includes(t))
      row.hidden = !hit
      if (hit) shown++
    }
    showCount(shown)
  })

  if (!sortable) return
  const body = table.tBodies[0]
  // Every body row into the first tbody, in the original order.
  for (const row of rows) body.appendChild(row)
  const original = rows.slice()
  const keys = new Map<HTMLTableRowElement, string[]>(
    rows.map((r) => [r, Array.from(r.cells, (c) => (c.textContent ?? '').trim())])
  )
  let active = -1
  let dir: 'ascending' | 'descending' | 'none' = 'none'

  const cells = Array.from(header.cells)
  const buttons = cells.map((th, c) => {
    th.setAttribute('aria-sort', 'none')
    const button = document.createElement('button')
    button.type = 'button'
    button.className = 'bvd-sort'
    while (th.firstChild) button.appendChild(th.firstChild)
    th.appendChild(button)
    button.addEventListener('click', () => sortBy(c))
    return button
  })
  buttons.forEach((b) => b.setAttribute('title', 'Sort by this column'))

  function sortBy(c: number): void {
    if (active !== c) {
      dir = 'ascending'
    } else {
      dir = dir === 'ascending' ? 'descending' : dir === 'descending' ? 'none' : 'ascending'
    }
    active = dir === 'none' ? -1 : c
    cells.forEach((th, i) => th.setAttribute('aria-sort', i === active ? dir : 'none'))
    const order = original.slice()
    if (dir !== 'none') {
      const descending = dir === 'descending'
      order.sort((a, b) =>
        compareCells(keys.get(a)![c] ?? '', keys.get(b)![c] ?? '', numeric[c], descending)
      )
    }
    for (const row of order) body.appendChild(row)
  }
}

let counter = 0

export function enhanceTables(root: ParentNode = document): void {
  const tables = root.querySelectorAll<HTMLTableElement>(`.vp-doc table:not([${DONE}])`)
  tables.forEach((t) => enhance(t, counter++))
}

// Re-applies after client-side navigation and when content mounts late.
export function watchTables(): void {
  enhanceTables()
  let queued = false
  new MutationObserver(() => {
    if (queued) return
    queued = true
    requestAnimationFrame(() => {
      queued = false
      enhanceTables()
    })
  }).observe(document.body, { childList: true, subtree: true })
}
