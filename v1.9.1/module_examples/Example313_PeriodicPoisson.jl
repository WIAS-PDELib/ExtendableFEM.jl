#=

# 312 : Periodic Poisson 3D
([source code](@__SOURCE_URL__))

This is a simple demonstration and validation of the new restriction based periodic boundary operator.

An unstructured cube grid is coupled along two axes periodically and restricted to non-zero Dirichlet values at the other to faces.

The result is a linear function and the error is measured directly.

Note that the result is independent of the periodic coupling! Therefore the correctness of the
new QR based column compression is covered by this example.

![](example313.png)
=#

module Example313_PeriodicPoisson

using ExtendableFEM
using ExtendableGrids
using SimplexGridFactory
using GridVisualize
using TetGen
using UnicodePlots
using StaticArrays
using LinearAlgebra
using Test #hide

const reg_left = 1
const reg_right = 2
const reg_front = 3
const reg_back = 4
const reg_bottom = 5
const reg_top = 6


function create_grid(; h)
    builder = SimplexGridBuilder(; Generator = TetGen)

    ## bottom points
    b00 = point!(builder, 0, 0, 0)
    b01 = point!(builder, 0, 1, 0)
    b10 = point!(builder, 1, 0, 0)
    b11 = point!(builder, 1, 1, 0)

    ## top points
    t00 = point!(builder, 0, 0, 1)
    t01 = point!(builder, 0, 1, 1)
    t10 = point!(builder, 1, 0, 1)
    t11 = point!(builder, 1, 1, 1)

    ## left face
    facetregion!(builder, reg_left)
    facet!(builder, b00, b01, t01, t00)

    ## right face
    facetregion!(builder, reg_right)
    facet!(builder, b10, b11, t11, t10)

    ## front face
    facetregion!(builder, reg_front)
    facet!(builder, b00, b10, t10, t00)

    ## back face
    facetregion!(builder, reg_back)
    facet!(builder, b01, b11, t11, t01)

    ## top face
    facetregion!(builder, reg_top)
    facet!(builder, t00, t10, t11, t01)

    ## bottom face
    facetregion!(builder, reg_bottom)
    facet!(builder, b00, b10, b11, b01)


    cellregion!(builder, 1)
    maxvolume!(builder, h)
    regionpoint!(builder, 0.5, 0.5, 0.5)

    return simplexgrid(builder)
end

function main(;
        Plotter = nothing,
        h = 1.0e-3,
        periodic = true,
        kwargs...
    )
    ## print options for better logs
    @info "selected options" periodic

    xgrid = create_grid(; h)

    FES = FESpace{H1P1{1}}(xgrid)

    ## problem description
    PD = ProblemDescription("Periodic Poisson Problem")
    u = Unknown("u"; name = "temperature")
    assign_unknown!(PD, u)

    assign_operator!(PD, BilinearOperator([grad(u)]; kwargs...))

    assign_restriction!(PD, BoundaryDataRestriction(u; value = -1, regions = [reg_bottom]))
    assign_restriction!(PD, BoundaryDataRestriction(u; value = +1, regions = [reg_top]))

    if periodic
        assign_restriction!(PD, CoupledDofsRestriction(u, reg_left, reg_right))
        assign_restriction!(PD, CoupledDofsRestriction(u, reg_front, reg_back))
    end

    ## solve
    sol, SC = solve(PD, FES; return_config = true, kwargs...)
    residual(SC) < 1.0e-10 || error("Residual is not zero!")

    @info "Lagrange residuals" SC.statistics[:restriction_residuals]

    function exact_error!(out, u, qpinfo)
        # exact solution here is u(x,y,z) = 2z - 1
        val = qpinfo.x[3] * 2 - 1
        out[1] = (val - u[1])^2
        return nothing
    end

    plt = plot([grid(u), id(u)], sol; Plotter, width = 1200, height = 800, scene3d = :LScene, slice = :y => 0.5)

    ErrorIntegrator = ItemIntegrator(exact_error!, [id(u)]; quadorder = 2)
    L2error = sqrt(sum(evaluate(ErrorIntegrator, sol)))

    @show L2error

    return L2error, plt

end

generateplots = ExtendableFEM.default_generateplots(Example313_PeriodicPoisson, "example313.png") #hide
function runtests()                                                                               #hide
    error1, _ = main(periodic = true)                                                             #hide
    @test error1 < 1.0e-12                                                                        #hide

    error2, _ = main(periodic = false)                                                            #hide
    @test error2 < 1.0e-12                                                                        #hide

    return nothing                                                                                #hide
end                                                                                               #hide

end # module
