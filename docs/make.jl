using Documenter
using ExampleJuggler
using CairoMakie
using ExtendableFEM

function make_all(; with_examples::Bool = true, modules = :all, run_examples::Bool = true, run_notebooks::Bool = false)

    module_examples = []

    if with_examples
        DocMeta.setdocmeta!(ExampleJuggler, :DocTestSetup, :(using ExampleJuggler); recursive = true)

        example_dir = joinpath(@__DIR__, "..", "examples")

        if modules === :all
            modules = readdir(example_dir)
        end

        #notebooks = ["PlutoTemplate.jl"
        #             "Example with Graphics" => "ExamplePluto.jl"]

        cleanexamples()

        module_examples = @docmodules(example_dir, modules, Plotter = CairoMakie)
        #html_examples = @docplutonotebooks(example_dir, notebooks, iframe=false)
        #pluto_examples = @docplutonotebooks(example_dir, notebooks, iframe=true)
    end

    makedocs(
        modules = [ExtendableFEM],
        sitename = "ExtendableFEM.jl",
        authors = "Christian Merdon, Jan Philipp Thiele",
        format = Documenter.HTML(; repolink = "https://github.com/WIAS-PDELib/ExtendableFEM.jl", mathengine = MathJax3()),
        clean = false,
        checkdocs = :none,
        warnonly = false,
        doctest = true,
        pages = [
            "Home" => "index.md",
            "Index" => "package_index.md",
            "Problem Description" => [
                "problemdescription.md",
                "tensordescription.md",
                "nonlinearoperator.md",
                "bilinearoperator.md",
                "linearoperator.md",
                "restrictions.md",
                "interpolateboundarydata.md",
                "homogeneousdata.md",
                "fixdofs.md",
                "combinedofs.md",
                "callbackoperator.md",
                "allindex.md",
            ],
            "Solving" => Any[
                "pdesolvers.md",
                "pdesolvers_dt.md",
                "parallel_assembly.md",
            ],
            "Postprocessing" => Any[
                "postprocessing.md",
                "itemintegrators.md",
                "faceinterpolator.md",
            ],
            #"Tutorial Notebooks" => notebooks,
            "Examples" => module_examples,
        ],
    )

    return cleanexamples()

end

make_all(; with_examples = true, run_examples = true, run_notebooks = false)

deploydocs(
    repo = "github.com/WIAS-PDELib/ExtendableFEM.jl",
)
