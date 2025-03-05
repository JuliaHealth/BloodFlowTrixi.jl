using BloodFlowTrixi
using Documenter

DocMeta.setdocmeta!(BloodFlowTrixi, :DocTestSetup, :(using BloodFlowTrixi); recursive=true)

makedocs(;
    modules=[BloodFlowTrixi],
    authors="yolhan83 <yolhan@laposte.net>",
    sitename="BloodFlowTrixi.jl",
    format=Documenter.HTML(;
        canonical="https://github.com/JuliaHealth/BloodFlowTrixi.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tuto.md",
        "Mathematics" => "math.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaHealth/BloodFlowTrixi.jl",
    devbranch="master",
)
