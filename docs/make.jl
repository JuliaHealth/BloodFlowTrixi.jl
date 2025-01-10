using BloodFlowTrixi
using Documenter

DocMeta.setdocmeta!(BloodFlowTrixi, :DocTestSetup, :(using BloodFlowTrixi); recursive=true)

makedocs(;
    modules=[BloodFlowTrixi],
    authors="yolhan83 <yolhan@laposte.net>",
    sitename="BloodFlowTrixi.jl",
    format=Documenter.HTML(;
        canonical="https://yolhan83.github.io/BloodFlowTrixi.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tuto.md"
    ],
)

deploydocs(;
    repo="github.com/yolhan83/BloodFlowTrixi.jl",
    devbranch="master",
)
