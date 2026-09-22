## The report's front matter (title, authors, dates, reading links, abstract,
## scope) is single-sourced in `README.md`, up to the `<!-- SHARED:END -->`
## marker. The home page shows the whole README and the offline
## `analysis.html` opens with the front matter, so both read it from here.

using BVDOutbreakSize: load_observations
import Dates

## The README with its live dates filled in: "Last updated" is the build date
## and "Data as of" is the loaded data cut-off, so a rebuild refreshes them
## without editing README.md.
function readme_with_dates()
    readme = read(joinpath(dirname(@__DIR__), "README.md"), String)
    built = Dates.format(Dates.today(), "d U yyyy")
    asof = Dates.format(load_observations().cutoff, "d U yyyy")
    return replace(
        readme,
        r"\*\*Last updated:\*\* [^.]*\." => "**Last updated:** $built.",
        r"\*\*Data as of:\*\* [^.]*\." => "**Data as of:** $asof."
    )
end

## The README up to the `<!-- SHARED:END -->` marker, dates filled in.
function front_matter()
    m = match(r"^(.*?)<!-- SHARED:END -->"s, readme_with_dates())
    m === nothing && error("README.md has no <!-- SHARED:END --> marker")
    return strip(m.captures[1])
end

## The abstract paragraph of the front matter, without its "**Abstract.**"
## label. It needs no dates, so it reads README.md directly rather than
## loading the data.
function readme_abstract()
    readme = read(joinpath(dirname(@__DIR__), "README.md"), String)
    m = match(r"^\*\*Abstract\.\*\* (.*?)\n\n"ms, readme)
    m === nothing && error("README.md has no **Abstract.** paragraph")
    return m.captures[1]
end
