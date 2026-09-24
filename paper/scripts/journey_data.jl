# Development-journey data for the preprint.
#
# Run from the worktree root:
#
#     julia --project=. paper/scripts/journey_data.jl
#
# Uses only the Julia standard library. Every number is read from git,
# gh or a committed CSV; the command behind each table is quoted next
# to the code that runs it. Nothing is estimated.
#
# Outputs (all under paper/data/):
#   commits_weekly.csv  first-parent commits on main and merged PRs,
#                       per ISO week, split by human / agent / other
#   code_size.csv       lines of code and fitted-stream count per tag
#   releases.csv        headline C_T per released tag
#   data_events.csv     data events that changed what the model fits
#   fit_cost.csv        joint gradient and fit wall-clock benchmarks
#
# n_streams_fitted in code_size.csv counts the observation-stream slots
# of `bvd_joint` at that tag. An observation submodel is a
# `@model function NAME` defined in src/models/observations.jl at that
# tag. At v1.0.0 and v1.1.0 every model lives in docs/examples/analysis.jl,
# so the set is the one that v1.2.0 moved into src/models/observations.jl
# (exports_model, deaths_model, cases_model, exports_deaths_model,
# exports_detection_timing_model). The `bvd_joint` definition is the text
# from `@model function bvd_joint(` to the first line that is exactly
# `end`, with comment lines removed. A slot is a keyword argument whose
# default is an observation submodel name (`deaths = deaths_model`), and
# each slot counts once even when two slots share a default. An
# observation submodel called directly in the body
# (`exports_deaths_model(...)`) adds one unless that name is already a
# slot or a slot default. Where a slot is fitted only on one branch of
# an `if` (v1.3.0 fits either the cumulative or the daily exports slot)
# it still counts, so the count is the number of stream slots the joint
# can score, not the number scored in one fit.

using Dates

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(REPO, "paper", "data")
const GH_REPO = "epiforecasts/BVDOutbreakSize"
mkpath(OUT)

git(args...) = read(Cmd(["git", "-C", REPO, args...]), String)
gh(args...) = read(Cmd(["gh", "-R", GH_REPO, args...]), String)
splitlines(s) = filter(!isempty, split(s, '\n'))

csvcell(x) = occursin(r"[,\"\n]", string(x)) ?
    "\"" * replace(string(x), "\"" => "\"\"") * "\"" : string(x)

function write_csv(name, header, rows)
    open(joinpath(OUT, name), "w") do io
        println(io, join(header, ","))
        for r in rows
            println(io, join(csvcell.(r), ","))
        end
    end
    return println(name, ": ", length(rows), " rows")
end

# ---------------------------------------------------------------------
# Who is who. Git author names and emails vary; GitHub logins do not.
# ---------------------------------------------------------------------
const HUMAN_LOGINS = Set(["seabbs", "sbfnk", "kathsherratt", "SamuelBrand1"])
const AGENT_LOGINS = Set(["seabbs-bot", "sbfnk-bot"])
const HUMAN_EMAILS = Set(
    [
        "contact@samabbott.co.uk", "s.e.abbott12@gmail.com",
        "sam.abbott@lshtm.ac.uk",
        "sebfnk@gmail.com", "sebastian.funk@lshtm.ac.uk",
        "48288458+SamuelBrand1@users.noreply.github.com",
    ]
)
const AGENT_EMAILS = Set(
    [
        "signin@samabbott.co.uk", "noreply@anthropic.com",
        "242615673+sbfnk-bot@users.noreply.github.com",
    ]
)

# :human, :agent or :other (dependabot and anything unrecognised).
function classify_commit(name, email)
    email in AGENT_EMAILS && return :agent
    email in HUMAN_EMAILS && return :human
    occursin("dependabot", email) && return :other
    occursin(r"bot|claude"i, name) && return :agent
    name in (
        "Sam Abbott", "Sebastian Funk", "Samuel Brand",
        "Kath Sherratt", "Katharine Sherratt",
    ) && return :human
    return :other
end

function classify_login(login)
    login in HUMAN_LOGINS && return :human
    login in AGENT_LOGINS && return :agent
    return :other
end

# ---------------------------------------------------------------------
# 1. commits_weekly.csv
# ---------------------------------------------------------------------
# git -C <repo> log --first-parent --format=%H%x09%as%x09%an%x09%ae origin/main
# %as is the author date (short ISO). origin/main is used because the
# local main branch belongs to another checkout and may lag.
mainref = isempty(git("branch", "-r", "--list", "origin/main")) ?
    "main" : "origin/main"
commit_rows = map(
    splitlines(
        git(
            "log", "--first-parent",
            "--format=%H%x09%as%x09%an%x09%ae", mainref
        )
    )
) do l
    h, d, n, e = split(l, '\t')
    (sha = h, date = Date(d), kind = classify_commit(n, e))
end
# gh -R epiforecasts/BVDOutbreakSize pr list --state merged --limit 2000
#    --json number,mergedAt,author
#    --jq '.[] | [.number, .mergedAt, .author.login] | @tsv'
pr_rows = map(
    splitlines(
        gh(
            "pr", "list", "--state", "merged",
            "--limit", "2000", "--json", "number,mergedAt,author",
            "--jq", ".[] | [.number, .mergedAt, .author.login] | @tsv"
        )
    )
) do l
    n, t, a = split(l, '\t')
    (
        number = parse(Int, n), date = Date(first(t, 10)),
        kind = classify_login(a),
    )
end
week(d) = firstdayofweek(d)  # ISO week, Monday
first_week = week(minimum(r.date for r in commit_rows))
last_week = week(max(today(), maximum(r.date for r in commit_rows)))
weeks = first_week:Week(1):last_week
count_kind(rows, w, k) = count(r -> week(r.date) == w && r.kind == k, rows)
write_csv(
    "commits_weekly.csv",
    [
        "week_start", "n_commits_human", "n_commits_agent", "n_commits_other",
        "n_prs_merged_human", "n_prs_merged_agent", "n_prs_merged_other",
    ],
    [
        (
            w, count_kind(commit_rows, w, :human),
            count_kind(commit_rows, w, :agent), count_kind(commit_rows, w, :other),
            count_kind(pr_rows, w, :human), count_kind(pr_rows, w, :agent),
            count_kind(pr_rows, w, :other),
        ) for w in weeks
    ]
)
println(
    "  commits on ", mainref, ": ", length(commit_rows),
    "; merged PRs: ", length(pr_rows)
)

# ---------------------------------------------------------------------
# 2. code_size.csv
# ---------------------------------------------------------------------
# git -C <repo> tag --list 'v*' 'V*'   (release tags only, both cases)
# git -C <repo> log -1 --format=%cs <tag>   (committer date of the tagged commit)
# git -C <repo> ls-tree -r --name-only <tag>
# git -C <repo> show <tag>:<path>   (line count of the blob)
tags = filter(
    t -> occursin(r"^[vV]\d+\.\d+\.\d+$", t),
    splitlines(git("tag", "--list", "v*", "V*"))
)
vernum(t) = VersionNumber(lstrip(t, ['v', 'V']))
sort!(tags; by = vernum)
tagdate(t) = Date(strip(git("log", "-1", "--format=%cs", t)))
countlines_at(t, path) = countlines(IOBuffer(git("show", "$t:$path")))
loc(t, files, pat) = sum(
    countlines_at(t, f) for f in files if occursin(pat, f);
    init = 0
)

const V12_OBS = [
    "exports_model", "deaths_model", "cases_model",
    "exports_deaths_model", "exports_detection_timing_model",
]
model_names(src) = [
    m.captures[1] for m in
        eachmatch(r"^@model function ([A-Za-z_0-9]+)"m, src)
]
function joint_body(src)
    lines = split(src, '\n')
    i = findfirst(l -> startswith(l, "@model function bvd_joint("), lines)
    i === nothing && return ""
    j = findnext(l -> l == "end", lines, i)
    return join(filter(l -> !startswith(strip(l), "#"), lines[i:j]), "\n")
end
function n_streams(t, files)
    obs = "src/models/observations.jl" in files ?
        model_names(git("show", "$t:src/models/observations.jl")) : V12_OBS
    jointfile = "src/models/joint.jl" in files ? "src/models/joint.jl" :
        "docs/examples/analysis.jl"
    body = joint_body(git("show", "$t:$jointfile"))
    slots = [
        (m.captures[1], m.captures[2]) for m in
            eachmatch(r"^\s*([A-Za-z_0-9]+)\s*=\s*([A-Za-z_0-9]+_model),?\s*$"m, body)
            if m.captures[2] in obs
    ]
    bound = union(Set(first.(slots)), Set(last.(slots)))
    direct = Set(
        m.captures[1] for m in
            eachmatch(r"\b([A-Za-z_0-9]+_model)\(", body)
            if m.captures[1] in obs && !(m.captures[1] in bound)
    )
    return length(slots) + length(direct)
end
size_rows = map(tags) do t
    files = splitlines(git("ls-tree", "-r", "--name-only", t))
    docs_pat = any(startswith("docs/pages/"), files) ?
        r"^docs/pages/.*\.jl$" : r"^docs/examples/.*\.jl$"
    (
        tagdate(t), t,
        loc(t, files, r"^src/.*\.jl$"),
        loc(t, files, r"^src/models/.*\.jl$"),
        loc(t, files, r"^test/.*\.jl$"),
        loc(t, files, docs_pat),
        n_streams(t, files),
    )
end
write_csv(
    "code_size.csv",
    [
        "date", "tag", "loc_src", "loc_models", "loc_test", "loc_docs_pages",
        "n_streams_fitted",
    ], size_rows
)

# ---------------------------------------------------------------------
# 3. releases.csv
# ---------------------------------------------------------------------
# Read from data/released_estimates.csv (columns: tag, date, model,
# median, lo30, hi30, lo60, hi60, lo90, hi90). Its tag column carries
# the results-<tag> prefix, which is stripped to join on the release
# tag. model_version_label maps the file's `model` column: integral ->
# closed-form, renewal -> renewal, and patch/province -> provincial.
label(m) = m == "integral" ? "closed-form" :
    m == "renewal" ? "renewal" :
    m in ("patch", "province", "provincial") ? "provincial" :
    error("unknown model label $m")
rel_lines = splitlines(
    read(
        joinpath(REPO, "data", "released_estimates.csv"),
        String
    )
)
rel_header = split(rel_lines[1], ',')
col(name) = findfirst(==(name), rel_header)
rel_rows = map(rel_lines[2:end]) do l
    f = split(l, ',')
    tag = replace(f[col("tag")], r"^results-" => "")
    (
        tag, f[col("date")], label(f[col("model")]), f[col("median")],
        f[col("lo90")], f[col("hi90")],
    )
end
sort!(rel_rows; by = r -> vernum(r[1]))
write_csv(
    "releases.csv",
    [
        "tag", "date", "model_version_label", "headline_C_T_median",
        "lower90", "upper90",
    ], rel_rows
)

# ---------------------------------------------------------------------
# 4. data_events.csv
# ---------------------------------------------------------------------
# Transcribed from paper/.notes/steering-and-data.md section 4 and
# paper/.notes/model-timeline.md section 2g. Only events with a URL and
# a change to what the model fits are kept. human_decision is "yes"
# when those notes cite a human comment or approval on the event.
const GH = "https://github.com/epiforecasts/BVDOutbreakSize/"
data_events = [
    (
        "2026-05-27", "WHO AFRO 01 to INSP",
        "Data source moved from the WHO AFRO sitrep to INSP sitreps via the INRB-UMIE transcription",
        "New source for every stream (v1.2.0)",
        GH * "issues/69", "yes",
    ),
    (
        "2026-05-28", "014",
        "INSP reclassifies suspects (1077 to 906 to 349) and withdraws the suspected-death headline; SitRep 012 re-issued as v2",
        "Suspected cumulative streams frozen at 26 May; confirmed cases and deaths fitted (v1.3.0)",
        GH * "issues/168", "yes",
    ),
    (
        "2026-06-02", "",
        "GeneXpert cannot detect BDBV so early samples were false negative",
        "Test sensitivity prior widened and lowered to Beta(6, 2)",
        GH * "issues/131#issuecomment-4605308124", "yes",
    ),
    (
        "2026-06-06", "",
        "Uganda export stream lags the cut-off",
        "Exports frozen",
        GH * "issues/141#issuecomment-4639599194", "yes",
    ),
    (
        "2026-06-11", "",
        "24h analysed counts available only for some provinces on some days",
        "Analysed-volume stream fitted with partial days accepted",
        GH * "pull/257", "yes",
    ),
    (
        "2026-06-19", "",
        "DHIS2 reclassification steps the isolation occupancy down",
        "Opt-in occupancy break days with fitted level offsets",
        GH * "pull/372", "yes",
    ),
    (
        "2026-07-05", "",
        "Tableau 6 patient-movement table absent from 30 June, returns 5 to 8 July",
        "Treatment flows fitted as a vintage; confirmed and suspect occupancy split",
        GH * "issues/373#issuecomment-4956011911", "yes",
    ),
    (
        "2026-07-12", "059",
        "Analytique format: onset epidemic curve appears, page-1 suspected-death subtitle stops",
        "Onset reporting triangle fitted (v1.11.0); suspected daily deaths frozen at 11 July",
        GH * "issues/431", "yes",
    ),
    (
        "2026-07-22", "069",
        "Harmonisation step integrates a provincial base (+369 net against +97 gross cases)",
        "Opt-in confirmed break day",
        GH * "pull/485", "no",
    ),
    (
        "2026-08-06", "084",
        "Short MVEBDB brief format: no suspects-du-jour tile, no onset figure, Tableau 6 dropped at 081",
        "Treatment flows frozen at SitRep 080; suspected daily cases frozen at 083",
        GH * "issues/562", "yes",
    ),
    (
        "2026-08-12", "090",
        "National samples-analysed total dropped",
        "tests_analysed_daily_history frozen at 11 August, resumed when the total returned",
        GH * "issues/570", "yes",
    ),
    (
        "2026-08-13", "088-089",
        "Isolation census drops about 200 in a day",
        "Occupancy break day",
        GH * "pull/567#issuecomment-5291322716", "yes",
    ),
    (
        "2026-08-20", "098",
        "Onset figure embedded losslessly and reads about 7% high",
        "SitRep 098 excluded from the onset curve",
        GH * "issues/594", "yes",
    ),
    (
        "2026-09-07", "115-116",
        "Onset date printed past the report date creates a delay-0 cell; joint R-hat 2.644",
        "Loader bounds the pair window; SitRep 117 briefly excluded then restored",
        GH * "issues/671", "yes",
    ),
    (
        "2026-09-16", "",
        "Daily suspect series resumed from the alert-validation table and the headline joint failed to converge",
        "Centred walks and widened background innovation SD (v2.1.0)",
        GH * "pull/761", "no",
    ),
]
@assert length(data_events) <= 15
write_csv(
    "data_events.csv",
    ["date", "sitrep", "event", "effect", "evidence_url", "human_decision"],
    data_events
)

# ---------------------------------------------------------------------
# 5. fit_cost.csv
# ---------------------------------------------------------------------
# Numbers copied verbatim from the PR bodies and the github-actions
# "Benchmark comparison vs main" comments, read with
#   gh -R epiforecasts/BVDOutbreakSize pr view <n> --json body,comments,mergedAt
# `arm` says which side of the PR's own comparison the number is.
# `model` says what was timed, because the benchmark shapes differ:
# the production joint at the fit's data, a synthetic joint, or the CI
# benchmark suite's synthetic joint. A cell is empty where the source
# records no number; ratios stated without an absolute are in `note`.
const PR = GH * "pull/"
fit_cost = [
    (
        "2026-09-08", "#656", "before", "production joint, main", "21.0", "",
        PR * "656", "median of 20 gradient evaluations, idle machine",
    ),
    (
        "2026-09-08", "#656", "after", "production joint, branch", "10.2", "",
        PR * "656", "median of 20 gradient evaluations, idle machine",
    ),
    (
        "2026-09-16", "#716", "", "headline patch fit on CI", "", "320",
        PR * "716", "wall-clock of the last headline fit job before the PR",
    ),
    (
        "2026-09-16", "#717", "before", "production joint, 229 parameters", "16.05", "",
        PR * "717", "loaded machine; after arm given only as ratio 0.729 (median)",
    ),
    (
        "2026-09-22", "#810", "before", "synthetic joint bvd_joint(20, 2, 3, 5, 1, 4, 10; breakpoint = 14), 106 parameters", "5.655", "",
        PR * "810", "median of 400, quiet machine, rules off",
    ),
    (
        "2026-09-22", "#810", "after", "synthetic joint bvd_joint(20, 2, 3, 5, 1, 4, 10; breakpoint = 14), 106 parameters", "4.065", "",
        PR * "810", "median of 400, quiet machine, rules on",
    ),
    (
        "2026-09-22", "#810", "before", "CI benchmark suite joint, Mooncake", "1.83", "",
        PR * "810#issuecomment-5778816174", "AirspeedVelocity minimum time per call, main",
    ),
    (
        "2026-09-22", "#810", "after", "CI benchmark suite joint, Mooncake", "1.56", "",
        PR * "810#issuecomment-5778816174", "AirspeedVelocity minimum time per call, PR",
    ),
    (
        "2026-09-23", "#837", "before", "CI benchmark suite joint, Mooncake", "1.57", "",
        PR * "837#issuecomment-5794384416", "AirspeedVelocity minimum time per call, main",
    ),
    (
        "2026-09-23", "#837", "after", "CI benchmark suite joint, Mooncake", "1.38", "",
        PR * "837#issuecomment-5794384416", "AirspeedVelocity minimum time per call, PR",
    ),
    (
        "2026-09-24", "#856", "before", "CI benchmark suite joint, Mooncake", "1.35", "",
        PR * "856#issuecomment-5801344552", "AirspeedVelocity minimum time per call, main",
    ),
    (
        "2026-09-24", "#856", "after", "CI benchmark suite joint, Mooncake", "1.06", "",
        PR * "856#issuecomment-5801344552", "AirspeedVelocity minimum time per call, PR",
    ),
]
write_csv(
    "fit_cost.csv",
    [
        "date", "tag_or_pr", "arm", "model", "joint_gradient_ms",
        "joint_fit_minutes", "source_url", "note",
    ], fit_cost
)
