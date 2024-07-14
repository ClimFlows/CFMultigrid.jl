using CFMultigrid
using Documenter

DocMeta.setdocmeta!(CFMultigrid, :DocTestSetup, :(using CFMultigrid); recursive=true)

makedocs(;
    modules=[CFMultigrid],
    authors="The ClimFlows contributors",
    sitename="CFMultigrid.jl",
    format=Documenter.HTML(;
        canonical="https://ClimFlows.github.io/CFMultigrid.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ClimFlows/CFMultigrid.jl",
    devbranch="main",
)
