# Table of change events for paper/main.qmd (@tbl-events), written to
# paper/generated/events_table.qmd from paper/data/change_events.csv.
#
#     julia --project=paper paper/scripts/events_table.jl
#
# Rows kept, at most MAX_ROWS ordered by date: the defects whose effect
# reached a released output, and the changes whose effect the record states
# as a move between two released medians with no other change sharing the
# release window. The evidence link and the human quote stay in the CSV.
# Each kept row's Event cell is the short label keyed by its evidence link
# below; its Effect cell is read from the CSV's stated move where there is
# one and otherwise from the label; its Releases live cell, for a defect,
# is the CSV's n_releases_live with the first release that carried it and
# the release that fixed it. The script stops if a kept row has no label,
# if a label carries a number the CSV row does not, if a defect row has no
# release span, or if more rows qualify than the cap.

using CSV
using DataFrames
using Dates

const MAX_ROWS = 20

## A defect whose effect never reached a released output.
const UNRELEASED = Regex(
    "not released|no released output|no committed row|fixed before|" *
        "fixed inside"
)
## A stated move between two released medians, `v1.5.0 4490 -> v1.6.0 2250`.
const MOVE = r"([vV]\d+\.\d+\.\d+) (\d+) -> ([vV]\d+\.\d+\.\d+) (\d+)"
## The move is shared with other changes merged between the same releases.
const SHARED = r"shared|not isolated"

released_defect(row) =
    row.kind == "defect" && !occursin(UNRELEASED, row.effect)
isolated_move(row) =
    row.kind != "defect" && occursin(MOVE, row.effect) &&
    !occursin(SHARED, row.effect)

## Short cell text by evidence link: (event, effect). The effect is used
## only where the CSV states no move.
const LABELS = Dict(
    "https://github.com/epiforecasts/BVDOutbreakSize/issues/69#issuecomment-4519775813" =>
        ("Data source moved to the INSP situation reports", ""),
    "https://github.com/epiforecasts/BVDOutbreakSize/issues/212#issuecomment-4647107920" =>
        (
        "Joint laboratory-queue fit sat above every single-stream fit",
        "Released with a warning banner, then replaced by the renewal model",
    ),
    "https://github.com/epiforecasts/BVDOutbreakSize/issues/144#issuecomment-4640589647" =>
        (
        "Renewal model with a weekly random walk on R replaced the " *
            "constant-growth model",
        "",
    ),
    "https://github.com/epiforecasts/BVDOutbreakSize/issues/235#issuecomment-4660405160" =>
        ("24h analysed laboratory volume fitted to anchor positivity", ""),
    "https://github.com/epiforecasts/BVDOutbreakSize/issues/406" =>
        (
        "Molecular-clock and growth priors updated to the " *
            "outbreak-specific genomic analysis",
        "",
    ),
    "https://github.com/epiforecasts/BVDOutbreakSize/issues/539" =>
        ("Growth-rate prior widened", ""),
    "https://github.com/epiforecasts/BVDOutbreakSize/issues/565#issuecomment-5314264922" =>
        (
        "Brief-format reports dropped the treatment flows and the " *
            "daily suspects, both frozen",
        "",
    ),
    "https://github.com/epiforecasts/BVDOutbreakSize/issues/606" =>
        (
        "Five defects in the one-week-ahead forecast, among them the " *
            "last innovation compounded over the horizon",
        "90% forecast coverage sat at one until the fix",
    ),
    "https://github.com/epiforecasts/BVDOutbreakSize/issues/623" =>
        (
        "Reporting-break corrections applied in the fit but not in " *
            "the scoring",
        "Frozen-fit forecasts scored an order of magnitude out",
    ),
    "https://github.com/epiforecasts/BVDOutbreakSize/pull/637#issuecomment-5531111868" =>
        (
        "Scoring fix carried confirmed but not occupancy break days " *
            "onto frozen fits",
        "Bed-stream scoring across breaks",
    ),
    "https://github.com/epiforecasts/BVDOutbreakSize/pull/637#issuecomment-5531338609" =>
        (
        "Bed shortfall could go negative under a positive offset",
        "Forecast bed-shortfall column",
    ),
    "https://github.com/epiforecasts/BVDOutbreakSize/pull/642#discussion_r3936181303" =>
        ("Onset triangle given a per-scan level error", ""),
    "https://github.com/epiforecasts/BVDOutbreakSize/pull/677#pullrequestreview-5200560195" =>
        (
        "Analysed laboratory volume capped below the modelled suspect " *
            "inflow, against the published methods",
        "Released fits before the provincial model carried the cap",
    ),
    "https://github.com/epiforecasts/BVDOutbreakSize/pull/412#issuecomment-4950535730" =>
        ("Four-province model became the headline", ""),
    "https://github.com/epiforecasts/BVDOutbreakSize/pull/745" =>
        ("Provincial release published with an unconverged headline fit", ""),
    "https://github.com/epiforecasts/BVDOutbreakSize/pull/742" =>
        (
        "Matched forecast scoring keyed without the made date",
        "Frozen confirmed-death skill 3.18 to 1.76",
    ),
    "https://github.com/epiforecasts/BVDOutbreakSize/issues/756" =>
        (
        "Provincial confirmed composition missed the report leg of " *
            "the delay",
        "Provincial split attributed to earlier days than the national " *
            "total",
    ),
    "https://github.com/epiforecasts/BVDOutbreakSize/issues/781" =>
        (
        "Sensitivity page still described the headline as three " *
            "provinces",
        "Released report text",
    ),
    "https://github.com/epiforecasts/BVDOutbreakSize/pull/875#pullrequestreview-5307701089" =>
        (
        "Onset digitiser read 10 to 23% short of the printed total " *
            "and mis-dated bars",
        "Onset curve data from v1.16.0 to v2.1.0",
    ),
    "https://github.com/epiforecasts/BVDOutbreakSize/issues/888#issuecomment-5822068785" =>
        (
        "Comparison table printed the current confirmed total as " *
            "observed at the calibration date of Chamla et al.",
        "Released sensitivity page, fix not yet released",
    ),
)

const FOUND_BY = Dict(
    "agent self-report" => "agent",
    "human review" => "human review",
    "human data check" => "human data check",
    "prospective evaluation" => "prospective evaluation",
    "automated review bot" => "automated reviewer",
    "CI gate" => "continuous integration",
    "outside contributor" => "outside contributor",
    "unrecorded" => "unrecorded",
)

## Thousands separators for the medians read from the CSV.
function fmt_count(n::Integer)
    s = string(n)
    parts = String[]
    while length(s) > 3
        pushfirst!(parts, s[(end - 2):end])
        s = s[1:(end - 3)]
    end
    pushfirst!(parts, s)
    return join(parts, ",")
end

## Every number in a label must appear in the CSV row it labels.
function check_numbers(label, row)
    text = replace(row.event * " " * row.effect, "," => "")
    for m in eachmatch(r"\d[\d.]*", label)
        occursin(m.match, text) ||
            error("`$(m.match)` in label `$label` is not in the CSV row")
    end
    return nothing
end

function effect_cell(row, label_effect)
    m = match(MOVE, row.effect)
    if m === nothing
        isempty(label_effect) &&
            error("no effect label for `$(row.evidence_url)`")
        return label_effect
    end
    from, a, to, b = m.captures
    return "Released median $(fmt_count(parse(Int, a))) to " *
        "$(fmt_count(parse(Int, b))) ($from to $to)"
end

## The tagged releases that carried a defect, from the first that did to
## the one before its fix; blank for the other rows.
function releases_live_cell(row)
    row.kind == "defect" || return ""
    n = row.n_releases_live
    ismissing(n) && error("no n_releases_live for `$(row.evidence_url)`")
    n == 0 && return "0, found before release"
    from, fix = row.live_from_release, row.fixed_in_release
    startswith(fix, "unreleased") && return "$n, from $from, fix unreleased"
    return "$n, from $from, fixed in $fix"
end

paper_dir() = normpath(joinpath(@__DIR__, ".."))

events = CSV.read(
    joinpath(paper_dir(), "data", "change_events.csv"), DataFrame
)
kept = filter(r -> released_defect(r) || isolated_move(r), events)
sort!(kept, :date)
nrow(kept) <= MAX_ROWS ||
    error("$(nrow(kept)) rows qualify, above the cap of $MAX_ROWS")

lines = String[
    "<!-- Written by paper/scripts/events_table.jl from " *
        "paper/data/change_events.csv; never edited by hand. -->",
    "",
    "| Date | Event | Effect | Found by | Decided by | Releases live |",
    "|---|---|---|---|---|---|",
]
for row in eachrow(kept)
    haskey(LABELS, row.evidence_url) ||
        error("no label for `$(row.evidence_url)`")
    event, effect = LABELS[row.evidence_url]
    check_numbers(event, row)
    check_numbers(effect, row)
    push!(
        lines,
        "| " * join(
            [
                Dates.format(row.date, "d U"), event,
                effect_cell(row, effect), FOUND_BY[row.detected_by],
                row.decided_by, releases_live_cell(row),
            ],
            " | "
        ) * " |"
    )
end
push!(lines, "")
push!(
    lines,
    ": Events in the development record whose effect reached a released " *
        "output, ordered by date (all 2026). " *
        "The $(nrow(kept)) rows are the subset of the $(nrow(events)) " *
        "recorded events that are defects in a released output or changes " *
        "whose move between two released medians the record states with " *
        "no other change sharing the release window; the full record is " *
        "in the repository (paper/data/change_events.csv). " *
        "Found by is the first record of the event on GitHub and " *
        "decided by is who set the response. " *
        "Releases live is the number of tagged releases that carried a " *
        "defect, from the first that did to the one before its fix, and " *
        "is blank for the other rows. " *
        "The evidence link for every row is in the record. {#tbl-events}"
)

out = joinpath(paper_dir(), "generated", "events_table.qmd")
mkpath(dirname(out))
write(out, join(lines, "\n") * "\n")
println("wrote ", nrow(kept), " rows to ", out)
