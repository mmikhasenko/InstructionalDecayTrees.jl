using InstructionalDecayTrees
using Documenter

const DOCS = @__DIR__
const QMD = joinpath(DOCS, "wigner_su2_so3.qmd")
const GFM = joinpath(DOCS, "wigner_su2_so3.md")
const TUTORIAL = joinpath(DOCS, "src", "wigner_tutorial.md")

function render_wigner_tutorial!()
    cd(DOCS) do
        run(`quarto render $(basename(QMD)) --to gfm`)
    end
    isfile(GFM) || error("expected Quarto output at $(GFM)")
end

function documenter_tutorial_page(gfm_path::AbstractString)
    body = read(gfm_path, String)
    body = replace(body, r"^# (.+)$"m => s"# [\1](@id wigner)"; count = 1)
    meta = "```@meta\nCurrentModule = InstructionalDecayTrees\nEditURL = \"../wigner_su2_so3.qmd\"\n```\n\n"
    return meta * body
end

DocMeta.setdocmeta!(
    InstructionalDecayTrees,
    :DocTestSetup,
    :(using InstructionalDecayTrees);
    recursive = true,
)

render_wigner_tutorial!()
write(TUTORIAL, documenter_tutorial_page(GFM))

makedocs(;
    modules = [InstructionalDecayTrees],
    authors = "Mikhail Mikhasenko and contributors",
    repo = "https://github.com/mmikhasenko/InstructionalDecayTrees.jl/blob/{commit}{path}#{line}",
    sitename = "InstructionalDecayTrees.jl",
    doctest = false,
    checkdocs = :none,
    format = Documenter.HTML(;
        canonical = "https://mmikhasenko.github.io/InstructionalDecayTrees.jl",
        repolink = "https://github.com/mmikhasenko/InstructionalDecayTrees.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Wigner angles: SO(3) vs SU(2)" => "wigner_tutorial.md",
    ],
)
