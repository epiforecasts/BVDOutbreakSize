#!/usr/bin/env julia
#
# Scan the per-province isolation occupancy, bed capacity and 24h
# patient-flow figures out of the INSP situation-report PDFs and emit
# `data/province_care_scanned.csv`.
#
# Where these figures live has changed three times.
#
#   To SitRep 080: a per-province occupancy table ("Tableau IV/V/5/6/7"
#   depending on vintage, headed "Mouvement des patients dans les
#   établissements de soins" or "Occupation des structures de soins"), one
#   column per province plus (usually) an "Ensemble"/"Global" total column.
#   The row labels are stable ("Patients en isolement (Fin J)", "Nombre de
#   lits", "Total admissions (24h)", "Total sorties (24h)", "Taux
#   d'occupation") but the province columns present, their order, and
#   whether the Ensemble column exists at all, all change as provinces join
#   (Sud-Kivu from SitRep 080) and as the table format is re-typeset. A cell
#   the report marks "ND" (non disponible) is not printed and is recorded as
#   an empty field, not zero. This era does not give a usable per-province
#   breakdown of the 24h recovered/died-in-care counts: the "Sorties"
#   breakdown rows wrap unpredictably across lines and are not attempted, so
#   `recovered_24h`/`deaths_incare_24h` are always empty here.
#
#   SitReps 081-083: the table is dropped. A two-column page layout carries
#   "PRISE EN CHARGE HOLISTIQUE (PECH)" prose in its right-hand column,
#   headed "→ Province : ..." bullets, mostly Ituri and Nord-Kivu.
#   `pdftotext -layout` interleaves both columns onto one physical line, so
#   the right column is recovered by slicing every line at the character
#   offset the "PECH" heading itself is printed at.
#
#   SitReps 084 onward: single-column prose returns, in a numbered section
#   ("1.5" to "1.7" depending on vintage) titled "Continuité des soins" or,
#   from SitRep 109, "Prise en charge holistique". One bullet per province,
#   for example (SitRep 125): "En Ituri, 89 nouvelles admissions ... contre
#   107 sorties dont 18 guéries. Ainsi, 442 patients sont hospitalisés pour
#   1015 lits, soit un taux d'occupation de 43,5 %." Phrasing is far freer
#   than the later vintages' near-fixed template, especially 084-106
#   ("l'occupation atteint 548 patients pour 991 lits", "37 admissions, 4
#   guéris et 8 décès notifiés ce jour", "8 patients sont en isolement ...
#   pour 12 lits"), so several phrasings are matched explicitly per field.
#   From around SitRep 124, Nord-Kivu narrows its bed count to "structures
#   normées": "376 patients sont hospitalisés dont 221 dans les structures
#   normées avec une capacité d'accueil de 228 lits" records
#   `patients_isolated = 376`, `patients_isolated_normed = 221`,
#   `beds = 228` (the only capacity figure printed that day).
#
# A province bullet that carries no care keyword at all (isolement,
# hospitalisé, admission, sortie, lits, occupation, guéri, décès) is not a
# care entry and contributes no row. A bullet that does carry one of these
# but from which no field at all can be extracted is reported to stderr as
# UNPARSED and contributes no row; the script exits non-zero if any
# occurred, exactly like scripts/scan_province_lab.jl. A field simply not
# printed within an otherwise-parsed bullet or row is left empty, never 0.
#
# Requires `pdftotext` (poppler-utils) on the PATH, and the sitrep PDFs
# already downloaded (`task download-sitreps`).
#
# Usage:
#
#   julia --project=scripts scripts/scan_province_care.jl
#   julia --project=scripts scripts/scan_province_care.jl path/to/pdfdir
#
# The output CSV path can be overridden for a dry run with the
# `CARE_OUT_CSV` environment variable; it defaults to
# `data/province_care_scanned.csv`.

using TOML
using Dates
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const PDF_DIR = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(ROOT, "data", "sitrep_pdfs")
const MANIFEST = joinpath(ROOT, "data", "observations.toml")
const SITREP_CSV = joinpath(ROOT, "data", "insp_sitrep_scanned.csv")
const OUT_CSV = get(
    ENV, "CARE_OUT_CSV", joinpath(ROOT, "data", "province_care_scanned.csv")
)

## Provinces in patch order: the first is the primary (origin) patch.
const PROVINCES = [
    "ituri", "nord_kivu", "sud_kivu", "haut_uele", "tshopo",
    "bas_uele", "sud_ubangi",
]
const DISPLAY = Dict(
    "ituri" => "Ituri", "nord_kivu" => "Nord-Kivu", "sud_kivu" => "Sud-Kivu",
    "haut_uele" => "Haut-Uele", "tshopo" => "Tshopo", "bas_uele" => "Bas-Uele",
    "sud_ubangi" => "Sud Ubangi",
)
## Folded province-name alternation, shared by every matcher in this file.
const PROV_ALT =
    "(ituri|nord[ -]?kivu|sud[ -]?kivu|haut[ -]?uele|tshopo|bas[ -]?uele|" *
    "sud[ -]?ubangi)"
const KEYS = Dict(
    "ituri" => "ituri", "nordkivu" => "nord_kivu",
    "sudkivu" => "sud_kivu", "hautuele" => "haut_uele",
    "tshopo" => "tshopo", "basuele" => "bas_uele",
    "sudubangi" => "sud_ubangi",
)
## Folded spellings of each province, to spot a clause that has moved on to
## another province mentioned in passing within the same block.
const NAMES = Dict(
    "ituri" => "ituri", "nord-kivu" => "nord_kivu", "nord kivu" => "nord_kivu",
    "sud-kivu" => "sud_kivu", "sud kivu" => "sud_kivu",
    "haut-uele" => "haut_uele", "haut uele" => "haut_uele",
    "tshopo" => "tshopo", "bas-uele" => "bas_uele", "bas uele" => "bas_uele",
    "sud-ubangi" => "sud_ubangi", "sud ubangi" => "sud_ubangi",
)
## A bullet in either the colon-style "→ Province : ..." (081-083) or the
## prose-style "En/Au/À la Province, ..." (084 onward). The leading prefix
## word is optional so both styles are matched by one regex.
const ENTRY = Regex(
    "^[^0-9a-z]*(?:en|au|aux|a la|a l.?)?\\s*" * PROV_ALT * "\\s*(?::+|,)"
)

"""
Strip accents and lower-case, so the many spellings and cases in the
sitreps collapse onto one folded form.
"""
function fold(s::AbstractString)
    return lowercase(Base.Unicode.normalize(String(s); stripmark = true))
end

# ---------------------------------------------------------------------
# Era 3/4 (SitRep 084 onward): single-column "Continuité des soins" /
# "Prise en charge holistique" prose section.
# ---------------------------------------------------------------------

"""
The lines of the care section, from its heading to (not including) the
"Défis" subsection, found by title rather than number since it has been
numbered 1.5, 1.6 and 1.7 across vintages.
"""
## A section number's own components never reach double digits in these
## reports (the deepest seen is "1.5.2"), unlike a wrapped sentence that
## happens to start with a number followed by an upper-case acronym
## ("34 PPL ont bénéficié de trois repas chauds ..."), which otherwise
## satisfies the heading regex just as well as a real heading does and
## ends the section early (SitRep 126's Bas-Uélé/Tshopo/Sud-Kivu bullets,
## cut off right after this exact phrasing).
plausible_heading_number(s::AbstractString) =
    all(p -> (v = tryparse(Int, p); v !== nothing && v <= 20), split(s, '.'))

function care_lines(lines::Vector{String})
    head = r"^[^0-9\p{L}]*(\d+(?:\.\d+)*)\.?\s*[-—–]?\s*(\p{Lu}.*)$"
    start, num = nothing, ""
    for (i, l) in enumerate(lines)
        m = match(head, l)
        m === nothing && continue
        plausible_heading_number(m[1]) || continue
        t = fold(m[2])
        (
            startswith(t, "continuite des soins") ||
                startswith(t, "prise en charge holistique")
        ) || continue
        start, num = i, m[1]
        break
    end
    start === nothing && return String[]
    out = String[]
    for l in lines[(start + 1):min(start + 60, length(lines))]
        f = fold(l)
        m = match(head, l)
        if m !== nothing && plausible_heading_number(m[1])
            (startswith(m[1], num) && !occursin("defis", fold(m[2]))) ||
                break
        end
        push!(out, String(l))
    end
    return out
end

"""
The care section split into one block of text per province, the same
wrapped-bullet accumulation `scripts/scan_province_lab.jl` uses.
"""
function province_entries(lines::Vector{String})
    got = Dict{String, String}()
    order = String[]
    current = nothing
    for line in lines
        isempty(strip(line)) && continue
        f = fold(line)
        m = match(ENTRY, f)
        if m !== nothing
            name = KEYS[replace(m[1], r"[ -]" => "")]
            current = haskey(got, name) ? nothing : name
            if current !== nothing
                got[name] = f[(m.offset + ncodeunits(m.match)):end]
                push!(order, name)
            end
        elseif occursin(r"^\s*[^0-9a-z\s]", f)
            current = nothing
        elseif current !== nothing
            got[current] = got[current] * " " * f
        end
    end
    return got, order
end

## A number attached to a keyword this file cares about; used to decide
## whether a bullet from which nothing was extracted is a genuine parse
## failure or simply carries no care figures at all (e.g. an ambulance or
## staffing note).
const CARE_KEYWORD = Regex(
    "(isolement|hospitalis|admission|sorti|lits?\\b|" *
        "occupation|gueri|dec[ei]s)"
)

const CUMUL_CLAUSE = r"cumul|s.eleve a|s.etablit a"

"""
`(field => value)` pairs read out of one province's prose entry, or
`:unparsed` when the text plainly carries a care figure that none of the
patterns below can read.
"""
function parse_prose_entry(entry::AbstractString, name::AbstractString)
    t = replace(entry, r"\s+" => " ")
    others = [k for (k, v) in NAMES if v != name]
    keepc(c) = !occursin(CUMUL_CLAUSE, c) && !any(o -> occursin(o, c), others)
    clauses = [c for c in split(t, r"[;.] ") if keepc(c)]
    text = join(clauses, " ; ")

    fields = Dict{Symbol, Any}()

    ## --- admissions -----------------------------------------------------
    m = match(r"(\d+)\)?\s*nouvelles?\s*admissions?", text)
    if m !== nothing
        fields[:admissions_24h] = parse(Int, m[1])
    elseif occursin(r"aucune nouvelle admission", text)
        fields[:admissions_24h] = 0
    elseif occursin(r"\bune nouvelle admission\b", text)
        fields[:admissions_24h] = 1
    elseif occursin(r"\bdeux nouvelles admissions\b", text)
        fields[:admissions_24h] = 2
    else
        for c in clauses
            mm = match(r"(\d+) admissions?\b(?! cumul)", c)
            mm === nothing && continue
            fields[:admissions_24h] = parse(Int, mm[1])
            break
        end
    end

    ## --- discharges (sorties) --------------------------------------------
    m = match(r"(\d+)\)?\s*sorties?\b", text)
    if m !== nothing
        fields[:discharges_24h] = parse(Int, m[1])
    elseif occursin(r"aucune sortie", text)
        fields[:discharges_24h] = 0
    elseif occursin(r"\bune(?: \(1\))? sortie\b", text)
        fields[:discharges_24h] = 1
    end

    ## --- recovered (guéris) ----------------------------------------------
    m = match(r"dont (\d+) gueri[es]*", text)
    m === nothing && (m = match(r"\((\d+) gueri[es]*\)", text))
    m === nothing &&
        (m = match(r"(\d+) (?:nouveaux? )?gueri[es]*\b", text))
    if m !== nothing
        fields[:recovered_24h] = parse(Int, m[1])
    elseif occursin(r"aucun gueri", text)
        fields[:recovered_24h] = 0
    elseif occursin(r"\(gueri[es]*\)", text) && haskey(fields, :discharges_24h)
        ## Bare "(guéries)" with no count of its own: every discharge that
        ## day was a recovery.
        fields[:recovered_24h] = fields[:discharges_24h]
    end

    ## --- deaths in care ----------------------------------------------------
    m = match(r"(\d+) dec[ei]s\b", text)
    if m !== nothing
        fields[:deaths_incare_24h] = parse(Int, m[1])
    elseif occursin(r"aucun dec[ei]s", text)
        fields[:deaths_incare_24h] = 0
    end

    ## --- patients isolated / normed structures / beds --------------------
    PAIR = Regex.(
        [
            "(\\d+)\\)? ?patients? (?:sont )?(?:hospitalis[ei]s?|" *
                "pris en charge en hospitalisation)[^;]*?" *
                "dont (\\d+) dans les structures normees[^;]*?" *
                "capacite d.accueil\\s*\\d{0,2}\\s*de (\\d+) lits",
            "(\\d+)\\)? ?patients? (?:sont )?(?:hospitalis[ei]s?|" *
                "pris en charge en hospitalisation)[^;]*?pour (\\d+) lits",
            "atteint (\\d+) patients pour (\\d+) lits",
            "atteint (\\d+) lits sur (\\d+)",
            "(\\d+)\\)? ?malades? (?:sont |restent |demeurent )?" *
                "(?:hospitalis[ei]s?|en isolement)[^;]*?(?:sur|pour) (\\d+) lits",
            "(\\d+)\\)? ?patients? (?:sont )?en isolement[^;]*?" *
                "(?:sur|pour) (\\d+) lits",
            "(\\d+)\\)? ?cas suspects? (?:restent |sont |demeurent )?en isolement" *
                "[^;]*?(?:sur|pour) (\\d+) lits",
        ]
    )
    for re in PAIR
        m = match(re, text)
        m === nothing && continue
        if length(m.captures) == 3
            fields[:patients_isolated] = parse(Int, m[1])
            fields[:patients_isolated_normed] = parse(Int, m[2])
            fields[:beds] = parse(Int, m[3])
        else
            fields[:patients_isolated] = parse(Int, m[1])
            fields[:beds] = parse(Int, m[2])
        end
        break
    end
    if !haskey(fields, :patients_isolated)
        SOLO = Regex.(
            [
                "(\\d+)\\)? ?patients? (?:sont )?" *
                    "(?:hospitalis[ei]s?|pris en charge en hospitalisation)",
                "(\\d+)\\)? ?malades? (?:sont |restent |demeurent )?" *
                    "(?:hospitalis[ei]s?|en isolement|en hospitalisation)",
                "(\\d+)\\)? ?patients? (?:sont |restent |demeurent )?en isolement",
                "(\\d+)\\)? ?cas (?:suspects?|confirmes?) (?:restent |sont |demeurent |" *
                    "est )?(?:en isolement|hospitalis[ei]s?|suivis)",
                "(\\d+)\\)? ?patient (?:confirme )?est en (?:cours de soins|" *
                    "isolement)",
                "(\\d+)\\)? ?suspects? isoles?",
                "pour (\\d+)\\)? ?malades? en fin de journee",
            ]
        )
        for re in SOLO
            m = match(re, text)
            m === nothing && continue
            fields[:patients_isolated] = parse(Int, m[1])
            break
        end
    end
    if !haskey(fields, :beds)
        m = match(r"pour (\d+) lits", text)
        m === nothing && (m = match(r"(\d+) lits disponibles?", text))
        m === nothing || (fields[:beds] = parse(Int, m[1]))
    end

    ## --- occupancy rate ----------------------------------------------------
    m = match(r"occupation[^0-9]{0,60}(\d+(?:,\d+)?)\s*%", text)
    m === nothing &&
        (m = match(r"taux d.occupation[^0-9]{0,60}(\d+(?:,\d+)?)\s*%", text))
    m !== nothing &&
        (fields[:occupancy_rate_pct] = parse(Float64, replace(m[1], ',' => '.')))

    ## A vintage that states outright that this province's care figures were
    ## not documented that day ("portés en ND", "ne sont pas documentés")
    ## carries no usable numbers by design, not by a parse failure.
    no_data = occursin(r"portes? en nd", text) ||
        occursin(r"ne sont pas documentes", text)
    isempty(fields) && !no_data && occursin(CARE_KEYWORD, text) &&
        occursin(r"\d", text) && return :unparsed
    return fields
end

# ---------------------------------------------------------------------
# Era 2 (SitReps 081-083): two-column PECH prose.
# ---------------------------------------------------------------------

"""
The right-hand-column text of the "PRISE EN CHARGE HOLISTIQUE (PECH)" box,
recovered by slicing every following line at the character offset the
heading itself starts at, until the next two-column heading pair (detected
as an upper-cased right-hand slice) is reached.
"""
function pech_lines(text::AbstractString)
    raw_lines = collect(split(text, '\n'))
    lines = [collect(l) for l in raw_lines]  # Vector{Vector{Char}}
    offset = nothing
    start = nothing
    for (i, l) in enumerate(raw_lines)
        ## `fold` is pure ASCII once accents are stripped, so the byte
        ## index `findfirst` returns for it is also a valid character
        ## offset into the original line (which `fold` preserves the
        ## character count of for plain accented Latin text).
        r = findfirst("prise en charge holistique", fold(l))
        r === nothing && continue
        offset = first(r)
        start = i
        break
    end
    start === nothing && return String[]
    out = String[]
    for chars in lines[(start + 1):min(start + 80, length(lines))]
        s = strip(String(chars))
        (isempty(s) || occursin(r"^\d+$", s) || s == "\x0c") && continue
        ## A page break can re-typeset the two columns at a slightly
        ## different width, so a fresh right-column bullet (starting "→")
        ## re-anchors the offset rather than trusting the header's.
        arrows = findall(==('→'), chars)
        near = filter(p -> p >= offset - 15, arrows)
        !isempty(near) && (offset = last(near))
        length(chars) < offset && continue
        right = strip(String(chars[offset:end]))
        isempty(right) && continue
        letters = filter(isletter, right)
        if !isempty(letters) && letters == uppercase(letters) &&
                length(letters) > 3
            break
        end
        push!(out, right)
    end
    return out
end

# ---------------------------------------------------------------------
# Era 1 (to SitRep 080): per-province occupancy table.
# ---------------------------------------------------------------------

const TABLE_NOISE = r"\((?:fin j|24 ?h|j-1|cumul)\)"
## A complete number, allowing an embedded thousands-separator space
## ("1 050"), a decimal comma or point, and an optional percent sign; or an
## "ND" placeholder. Used both to find where each column starts and to read
## a value out of it, so the two always agree on what counts as one token.
const NUMTOK = r"\d{1,3}(?: \d{3})*(?:[.,]\d+)?%?|nd\b"

"""
`(names, has_ensemble, header_line)` for the occupancy table's column
header: the province names in left-to-right order, whether an
Ensemble/Global column follows them, and the line they were found on
(found by scanning backward from the anchor "Fin J" row for the nearest
line naming at least two provinces).
"""
function table_header(lines::Vector{String}, anchor::Int)
    prov_re = Regex(PROV_ALT)
    for i in anchor:-1:max(1, anchor - 30)
        f = fold(lines[i])
        ms = collect(eachmatch(prov_re, f))
        length(ms) >= 2 || continue
        names = [KEYS[replace(m.match, r"[ -]" => "")] for m in ms]
        has_ens = occursin(r"ensemble|global", f)
        return names, has_ens, i
    end
    return String[], false, nothing
end

"""
Column start `positions`, one per province plus one more for a trailing
Ensemble/Global column when `nprov_cols` says one is expected, read from
the first "complete" row at or after `from` (every column present, so
there is no ambiguity about which token is which column) rather than from
the header text itself.

A header word's own position is not a reliable column boundary: values are
right-aligned under a header label that can be much narrower than they
are, so slicing at the label's start clips a wide value's leading digit
(SitRep 056's occupancy row, "89,9" read as "9,9" under "Taux
d'occupation"), or, with a header naming several similarly-short provinces
close together, shifts a whole row over by one column (SitRep 079's
"Nombre de lits" row, where Nord-Kivu's "141" was clipped to "1" and the
remainder "41" bled into Haut-Uélé's own value). A row with no blank cell
prints every value in its true column, so its own token start positions
are exact wherever pdftotext places them.
"""
function reference_positions(
        lines::Vector{String}, from::Int, to::Int, nprov_cols::Int
    )
    for i in from:to
        f = fold(lines[i])
        occursin("cumul", f) && continue
        toks = collect(eachmatch(NUMTOK, f))
        length(toks) == nprov_cols && return [m.offset for m in toks]
    end
    return Int[]
end

"""
Column boundaries from a header's name-start `positions`: the midpoint
between each pair of consecutive names, rather than a name's own start,
because a numeric cell is right-aligned under its (usually short) header
word and so commonly starts to the left of it — slicing at the header's own
start position clips the leading digit of a wide value (SitRep 056's "89,9"
read as "9,9" under "Taux d'occupation", the value starting one character
left of where "Ituri" itself starts). The first boundary is the start of
the line and the last is one past its end, so the leftmost and rightmost
columns are never clipped either.
"""
function column_bounds(positions::Vector{Int})
    n = length(positions)
    bounds = Vector{Int}(undef, n + 1)
    bounds[1] = 1
    for i in 1:(n - 1)
        bounds[i + 1] = (positions[i] + positions[i + 1]) ÷ 2
    end
    bounds[n + 1] = typemax(Int)
    return bounds
end

"""
The value (a digit string, possibly with a decimal comma or point and a
percent sign, or `nothing` for a blank or "ND" cell) sliced out of `line`
between character positions `lo` (inclusive) and `hi` (exclusive). Label
text is stripped first so an hours annotation like "(24h)" glued to the row
label is never misread as a value, and a thousands separator ("1 049") is
merged before extracting so it is not split into two numbers.
"""
function column_value(line::AbstractString, lo::Int, hi::Int)
    chars = collect(line)
    lo = max(lo, 1)
    lo > length(chars) && return nothing
    seg = String(chars[lo:min(hi - 1, length(chars))])
    f = replace(fold(seg), TABLE_NOISE => "")
    m = match(NUMTOK, f)
    (m === nothing || m.match == "nd") && return nothing
    return replace(String(m.match), "%" => "", " " => "")
end

"""
`Dict(province => value_or_nothing)` for the first line at or after `from`
(and at or before `to`) matching `label_pred`, sliced at the header's
column positions; checks the following line too, for the vintages that
print the row label alone and its values on the next line. An empty `Dict`
means the label was never found this vintage.
"""
function find_row(
        lines::Vector{String}, from::Int, to::Int, label_pred,
        header::Vector{String}, positions::Vector{Int}
    )
    bounds = column_bounds(positions)
    slice(line) = Dict{String, Union{String, Nothing}}(
        p => column_value(line, bounds[i], bounds[i + 1])
            for (i, p) in enumerate(header)
    )
    for i in from:to
        f = fold(lines[i])
        occursin("cumul", f) && continue
        label_pred(f) || continue
        out = slice(lines[i])
        any(v -> v !== nothing, values(out)) && return out
        i < to && (out = slice(lines[i + 1]))
        any(v -> v !== nothing, values(out)) && return out
        return out
    end
    return Dict{String, Union{String, Nothing}}()
end

"""
`(province => (patients_isolated, beds, occupancy_pct, admissions,
discharges))` string-or-nothing tuples, read from the occupancy table
anchored on the "Patients en isolement (Fin J)" row. `nothing` for a table
that could not be found in `lines` at all (a different vintage's era).
"""
function table_entries(lines::Vector{String})
    anchor = nothing
    for (i, l) in enumerate(lines)
        f = fold(l)
        occursin("isolement", f) && occursin("(fin", f) && (anchor = i; break)
    end
    anchor === nothing && return nothing
    header, has_ens, header_line = table_header(lines, anchor)
    isempty(header) && return nothing

    to = min(anchor + 30, length(lines))
    from = max(1, anchor - 30)
    nprov_cols = length(header) + (has_ens ? 1 : 0)
    positions = reference_positions(lines, header_line + 1, to, nprov_cols)
    if isempty(positions)
        ## No row in the table happens to have every column filled; fall
        ## back to the header's own word-start positions, which are at
        ## least approximately right.
        prov_re = Regex(PROV_ALT)
        f = fold(lines[header_line])
        positions = [m.offset for m in eachmatch(prov_re, f)]
        me = match(r"ensemble|global", f)
        me !== nothing && push!(positions, me.offset)
    end

    ## The row label "(Fin J)" is itself sometimes wrapped mid-parenthesis
    ## ("... isolement (Fin" / values / "J)"), pushing the values this row's
    ## own label line down by one; the two-line fallback in `find_row`
    ## covers it once the search window includes that next line.
    r_isolated = find_row(
        lines, anchor, min(anchor + 1, length(lines)), f -> true,
        header, positions
    )
    ## SitRep 031 prints the values line ABOVE the label instead ("313 44 6
    ## 363" / "Patients en isolement (Fin J)", no numbers on the label's own
    ## line), so only try the line before once neither the label's own line
    ## nor the one after it produced a value.
    if isempty(r_isolated) || all(v -> v === nothing, values(r_isolated))
        anchor > 1 && (
            r_isolated = find_row(
                lines, anchor - 1, anchor - 1, f -> true, header, positions
            )
        )
    end

    to = min(anchor + 30, length(lines))
    from = max(1, anchor - 30)
    ## "Nombre de lits" is printed before "Patients en isolement (Fin J)" in
    ## most vintages but after it (following "dont confirmés"/"dont
    ## suspects") from around SitRep 062, so both sides of the anchor are
    ## searched.
    r_beds = find_row(
        lines, from, anchor - 1, f -> occursin("nombre de lits", f),
        header, positions
    )
    isempty(r_beds) && (
        r_beds = find_row(
            lines, anchor, min(anchor + 8, length(lines)),
            f -> occursin("nombre de lits", f), header, positions
        )
    )

    ## A genuine table row prints one percentage per column; a prose
    ## sentence that happens to mention "taux d'occupation" right after the
    ## table (SitRep 064 has no occupancy row at all, and its very next
    ## line is "Au total, 722 patients étaient en isolement, ... pour un
    ## taux d'occupation global de ...") carries at most one, and applying
    ## the table's column positions to it reads whatever number the prose
    ## happens to print first as if it were the first province's rate.
    ## Also kept to a tight window after the anchor, since every table that
    ## does carry this row prints it within a line or two of "Fin J".
    r_occ = find_row(
        lines, anchor, min(anchor + 6, length(lines)),
        f -> occursin(r"taux d.occupation", f) && count(==('%'), f) >= 2,
        header, positions
    )

    r_adm = find_row(
        lines, from, to,
        f -> occursin("total admissions", f) && occursin("24", f),
        header, positions
    )
    isempty(r_adm) && (
        r_adm = find_row(
            lines, from, to,
            f -> occursin("total admissions", f) && !occursin("24", f),
            header, positions
        )
    )

    r_sor = find_row(
        lines, anchor, to,
        f -> occursin("total sorties", f) && occursin("24", f),
        header, positions
    )
    isempty(r_sor) && (
        r_sor = find_row(
            lines, anchor, to,
            f -> occursin("total sorties", f) && !occursin("24", f),
            header, positions
        )
    )

    ## Beds may also come from a fraction printed under the occupancy row,
    ## e.g. "(237/315)", when no "Nombre de lits" row is printed this
    ## vintage.
    if isempty(r_beds) || all(v -> v === nothing, values(r_beds))
        frac_line = nothing
        for i in anchor:to
            f = fold(lines[i])
            occursin(r"\(\d+/\d+\)", f) && (frac_line = lines[i]; break)
        end
        if frac_line !== nothing
            fracs = [
                m.captures[2] for
                    m in eachmatch(r"\((\d+)/(\d+)\)", fold(frac_line))
            ]
            length(fracs) == length(header) &&
                (r_beds = Dict(zip(header, fracs)))
        end
    end

    return (
        isolated = r_isolated, beds = r_beds, occupancy = r_occ,
        admissions = r_adm, discharges = r_sor, header = header,
    )
end

# ---------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------

struct Row
    sitrep::String
    report_date::String
    province::String
    patients_isolated::Union{Int, Nothing}
    patients_isolated_normed::Union{Int, Nothing}
    beds::Union{Int, Nothing}
    occupancy_rate_pct::Union{Float64, Nothing}
    admissions_24h::Union{Int, Nothing}
    discharges_24h::Union{Int, Nothing}
    recovered_24h::Union{Int, Nothing}
    deaths_incare_24h::Union{Int, Nothing}
    source_section::String
    quote_text::String
end

csv_field(x::Nothing) = ""
csv_field(x::Integer) = string(x)
csv_field(x::AbstractFloat) = string(x)
function csv_field(x::AbstractString)
    q = replace(x, "\"" => "\"\"")
    return "\"$q\""
end

function trimmed_quote(s::AbstractString, n::Int = 300)
    q = replace(strip(s), r"\s+" => " ")
    return length(q) > n ? q[1:n] : q
end

"""
Map sitrep number to report date, from the hand-scanned headline CSV.
"""
function sitrep_dates()
    dates = Dict{String, String}()
    for (i, line) in enumerate(eachline(SITREP_CSV))
        i == 1 && continue
        f = split(line, ',')
        length(f) < 2 && continue
        dates[f[1]] = f[2]
    end
    return dates
end

function process_table_era!(
        rows::Vector{Row}, unparsed, sr::AbstractString, date::AbstractString, lines::Vector{String}
    )
    te = table_entries(lines)
    te === nothing && return false
    for p in te.header
        vals = Dict{Symbol, Any}()
        for (field, r) in (
                (:patients_isolated, te.isolated), (:beds, te.beds),
                (:occupancy_rate_pct, te.occupancy),
                (:admissions_24h, te.admissions), (:discharges_24h, te.discharges),
            )
            haskey(r, p) || continue
            v = r[p]
            v === nothing && continue
            vals[field] = field == :occupancy_rate_pct ?
                parse(Float64, replace(v, ',' => '.')) : parse(Int, v)
        end
        isempty(vals) && continue
        push!(
            rows,
            Row(
                sr, date, DISPLAY[p],
                get(vals, :patients_isolated, nothing), nothing,
                get(vals, :beds, nothing),
                get(vals, :occupancy_rate_pct, nothing),
                get(vals, :admissions_24h, nothing),
                get(vals, :discharges_24h, nothing),
                nothing, nothing,
                "occupancy table (Mouvement des patients / Occupation des " *
                    "structures de soins)",
                "",
            ),
        )
    end
    return true
end

function process_prose_era!(
        rows::Vector{Row}, unparsed, sr::AbstractString, date::AbstractString,
        lines::Vector{String}, section::String
    )
    entries, order = province_entries(lines)
    isempty(order) && return false
    for name in order
        entry = entries[name]
        r = parse_prose_entry(entry, name)
        if r === :unparsed
            push!(unparsed, (sr, DISPLAY[name], trimmed_quote(entry, 200)))
            continue
        end
        isempty(r) && continue
        push!(
            rows,
            Row(
                sr, date, DISPLAY[name],
                get(r, :patients_isolated, nothing),
                get(r, :patients_isolated_normed, nothing),
                get(r, :beds, nothing),
                get(r, :occupancy_rate_pct, nothing),
                get(r, :admissions_24h, nothing),
                get(r, :discharges_24h, nothing),
                get(r, :recovered_24h, nothing),
                get(r, :deaths_incare_24h, nothing),
                section, trimmed_quote(entry),
            ),
        )
    end
    return true
end

function main()
    isdir(PDF_DIR) || error(
        "no sitrep PDFs at $(PDF_DIR); run `task download-sitreps` first."
    )
    Sys.which("pdftotext") === nothing &&
        error("pdftotext not found; install poppler-utils.")

    dates = sitrep_dates()
    rows = Row[]
    unparsed = Tuple{String, String, String}[]
    processed = String[]
    skipped = String[]

    ## Several sitreps have more than one archived PDF (a hyphen- and an
    ## underscore-named file, or an original plus a "_v2" reissue that
    ## shares the original's digits). Only one file per sitrep number is
    ## scanned, the same last-one-sorted-wins rule
    ## `scripts/scan_province_lab.jl` gets for free from overwriting a Dict
    ## keyed by date; here the rows are accumulated in a plain vector, so
    ## the choice of file is made explicitly before scanning.
    by_sr = Dict{String, String}()
    for path in sort(
            filter(f -> endswith(f, ".pdf"), readdir(PDF_DIR; join = true))
        )
        m = match(r"(\d+)[_-]2026", basename(path))
        m === nothing && continue
        by_sr[m[1]] = path
    end

    for sr in sort(collect(keys(by_sr)); by = s -> parse(Int, s))
        haskey(dates, sr) || continue
        date = dates[sr]
        path = by_sr[sr]

        text = read(`pdftotext -layout $path -`, String)
        lines = String.(split(text, '\n'))

        ## The occupancy table, the PECH box and the prose section are
        ## alternative eras reporting the same figures, not complementary
        ## sources, and running every era unconditionally would double-count
        ## a vintage that carries more than one. The table is tried first
        ## and, when found, wins outright: its anchor row label ("Patients
        ## en isolement (Fin J)") is specific enough to the old table that
        ## it is never seen in the modern prose. This matters because a
        ## clutch of table-era vintages (050-079-ish) number their
        ## "Prise en charge holistique" section 4.x, the same title the
        ## post-084 prose section uses, and print the table directly under
        ## it with each province's name alone on its own line ("Ituri :")
        ## ahead of a prose paragraph and then the table itself; the prose
        ## bullet splitter reads that bare heading as an entry start and
        ## swallows the whole table into it, misreading the table's
        ## national total row as if it were that one province's figure.
        ## Trying the table anchor first, before that colon is ever read as
        ## a bullet, avoids the misattribution rather than detecting it
        ## after the fact.
        got_table = table_entries(lines) !== nothing &&
            process_table_era!(rows, unparsed, sr, date, lines)
        got_prose = false
        got_pech = false
        if !got_table
            got_prose = process_prose_era!(
                rows, unparsed, sr, date, care_lines(lines),
                "Continuité des soins / Prise en charge holistique section"
            )
            if !got_prose
                got_pech = process_prose_era!(
                    rows, unparsed, sr, date, pech_lines(text),
                    "PECH two-column box"
                )
            end
        end

        if got_prose || got_pech || got_table
            push!(processed, sr)
        else
            push!(skipped, sr)
        end
    end

    ## --- Write the CSV ---------------------------------------------------
    open(OUT_CSV, "w") do io
        println(
            io,
            "sitrep,report_date,province,patients_isolated," *
                "patients_isolated_normed,beds,occupancy_rate_pct," *
                "admissions_24h,discharges_24h,recovered_24h," *
                "deaths_incare_24h,source_section,quote"
        )
        for r in rows
            println(
                io,
                join(
                    [
                        csv_field(r.sitrep), csv_field(r.report_date),
                        csv_field(r.province),
                        r.patients_isolated === nothing ? "" :
                            csv_field(r.patients_isolated),
                        r.patients_isolated_normed === nothing ? "" :
                            csv_field(r.patients_isolated_normed),
                        r.beds === nothing ? "" : csv_field(r.beds),
                        r.occupancy_rate_pct === nothing ? "" :
                            csv_field(r.occupancy_rate_pct),
                        r.admissions_24h === nothing ? "" :
                            csv_field(r.admissions_24h),
                        r.discharges_24h === nothing ? "" :
                            csv_field(r.discharges_24h),
                        r.recovered_24h === nothing ? "" :
                            csv_field(r.recovered_24h),
                        r.deaths_incare_24h === nothing ? "" :
                            csv_field(r.deaths_incare_24h),
                        csv_field(r.source_section), csv_field(r.quote_text),
                    ],
                    ","
                )
            )
        end
    end
    println("Wrote $(length(rows)) rows to $OUT_CSV")
    println(
        "$(length(processed)) sitreps yielded at least one row; " *
            "$(length(skipped)) yielded none: ", join(skipped, ", ")
    )

    ## --- Cross-check against the national isolation_history --------------
    raw = TOML.parsefile(MANIFEST)
    nat = Dict(
        String(d) => v for
            (d, v) in zip(
                raw["isolation_history"]["dates"],
                raw["isolation_history"]["values"]
            )
    )
    by_date = Dict{String, Dict{String, Int}}()
    for r in rows
        r.patients_isolated === nothing && continue
        d = get!(by_date, r.report_date, Dict{String, Int}())
        d[r.province] = get(d, r.province, 0) + r.patients_isolated
    end
    println("\nCross-check against national isolation_history\n")
    @printf("%11s %5s %6s %6s  %s\n", "date", "nat", "prov", "diff", "provinces")
    println("-"^80)
    exceed = String[]
    for d in sort(collect(keys(nat)))
        haskey(by_date, d) || continue
        s = sum(values(by_date[d]))
        n = nat[d]
        note = s > n ? " <-- EXCEEDS NATIONAL" : ""
        !isempty(note) && push!(exceed, d)
        @printf(
            "%11s %5d %6d %6d  %s%s\n", d, n, s, s - n,
            join(sort(collect(keys(by_date[d]))), "+"), note
        )
    end

    ## --- Report unparsed entries -------------------------------------------
    if !isempty(unparsed)
        for (sr, prov, txt) in unparsed
            println(stderr, "UNPARSED sitrep=$sr province=$prov text=$txt")
        end
        println(
            stderr,
            "\n$(length(unparsed)) unparsed entr" *
                (length(unparsed) == 1 ? "y" : "ies") * "."
        )
    end
    if !isempty(exceed)
        println(
            stderr,
            "\n$(length(exceed)) date(s) where the province sum exceeds " *
                "the national isolation_history tile: ", join(exceed, ", ")
        )
    end
    isempty(unparsed) || exit(1)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
