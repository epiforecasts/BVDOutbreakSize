// Large report tables get a compact style, a scroll box with a sticky
// header, a row filter and sortable columns. Tables with no more than
// `MIN_ROWS` body rows are left alone, as is every table without
// JavaScript. The enhancement only adds a wrapper and a toolbar around the
// table and reorders or hides its rows, so the server-rendered markup is
// unchanged until the page has hydrated.

// About a screen of rows at the compact size. The province tables (eight
// provinces) and the ten-row summaries stay as they are.
const MIN_ROWS = 10

const DONE = 'data-bvd-table'

// A cell sorts by its first number: "1.2 (0.8–1.9)", "12%", "-2.21",
// "286,637" and "1e-3" all qualify. ISO dates sort as text.
const NUMBER = /^[^\d\p{L}+\-−.]*([+\-−]?(?:\d{1,3}(?:,\d{3})+|\d+)?(?:\.\d+)?(?:e[+\-−]?\d+)?)/iu
const DATE = /^\s*\d{4}-\d{2}-\d{2}/

function cellNumber(text: string): number {
  if (DATE.test(text)) return NaN
  const m = NUMBER.exec(text)
  if (!m || !/\d/.test(m[1])) return NaN
  return Number(m[1].replace(/,/g, '').replace(/−/g, '-'))
}

function isBlank(text: string): boolean {
  return /^\s*(|-|—|–|missing|NaN|n\/a)\s*$/i.test(text)
}

interface Parts {
  header: HTMLTableRowElement
  rows: HTMLTableRowElement[]
}

// The header row whose cells label the columns, and the body rows. A raw
// HTML table without a <thead> has its leading all-<th> row as the header.
function tableParts(table: HTMLTableElement): Parts | null {
  const rows: HTMLTableRowElement[] = []
  for (const body of Array.from(table.tBodies)) {
    rows.push(...Array.from(body.rows))
  }
  const header = table.tHead ? table.tHead.rows[0] : rows.shift()
  if (!header || header.querySelector('td')) return null
  return { header, rows }
}

function hasSpans(table: HTMLTableElement): boolean {
  return table.querySelector('[colspan]:not([colspan="1"]), [rowspan]:not([rowspan="1"])') !== null
}

function isNumericColumn(rows: HTMLTableRowElement[], col: number): boolean {
  let numeric = 0
  let filled = 0
  for (const row of rows) {
    const cell = row.cells[col]
    if (!cell) continue
    const text = cell.textContent ?? ''
    if (isBlank(text)) continue
    filled++
    if (!Number.isNaN(cellNumber(text))) numeric++
  }
  return filled > 0 && numeric / filled >= 0.8
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
    numeric.push(isNumericColumn(rows, c))
    const cls = numeric[c] ? 'bvd-num' : 'bvd-text'
    for (const row of [...Array.from(table.tHead!.rows), ...rows]) {
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
  const collator = new Intl.Collator(undefined, { numeric: true, sensitivity: 'base' })
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

  function compare(c: number, a: HTMLTableRowElement, b: HTMLTableRowElement): number {
    const x = keys.get(a)![c] ?? ''
    const y = keys.get(b)![c] ?? ''
    // Blank cells go last in both directions.
    const bx = isBlank(x)
    const by = isBlank(y)
    if (bx || by) return bx === by ? 0 : bx ? Infinity : -Infinity
    if (numeric[c]) {
      const nx = cellNumber(x)
      const ny = cellNumber(y)
      const ix = Number.isNaN(nx)
      const iy = Number.isNaN(ny)
      if (ix || iy) return ix === iy ? collator.compare(x, y) : ix ? Infinity : -Infinity
      return nx - ny
    }
    return collator.compare(x, y)
  }

  function sortBy(c: number): void {
    if (active !== c) {
      dir = 'ascending'
    } else {
      dir = dir === 'ascending' ? 'descending' : dir === 'descending' ? 'none' : 'ascending'
    }
    active = dir === 'none' ? -1 : c
    cells.forEach((th, i) => th.setAttribute('aria-sort', i === active ? dir : 'none'))
    let order = original.slice()
    if (dir !== 'none') {
      const sign = dir === 'ascending' ? 1 : -1
      order.sort((a, b) => {
        const v = compare(c, a, b)
        // Infinity marks a blank or non-numeric cell that stays last.
        if (!Number.isFinite(v)) return v > 0 ? 1 : -1
        return sign * v
      })
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
