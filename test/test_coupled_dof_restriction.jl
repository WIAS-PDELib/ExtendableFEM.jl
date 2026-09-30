using StaticArrays
using ExtendableGrids
using TetGen

function run_coupled_dofs_restriction_tests()

    # unit cube rotated by angle θ around its center (z-axis),
    # such that coupled faces are diagonal planes
    function create_diagonal_cube_grid(; h, θ)
        c, s = cos(θ), sin(θ)
        rot(x, y, z) = @SVector [
            c * x - s * y,
            s * x + c * y,
            z,
        ]

        builder = SimplexGridBuilder(; Generator = TetGen)

        ## bottom points
        b00 = point!(builder, rot(-1000, -1, -1)...)
        b01 = point!(builder, rot(-1000, 1, -1)...)
        b10 = point!(builder, rot(1000, -1, -1)...)
        b11 = point!(builder, rot(1000, 1, -1)...)

        ## top points
        t00 = point!(builder, rot(-1000, -1, 1)...)
        t01 = point!(builder, rot(-1000, 1, 1)...)
        t10 = point!(builder, rot(1000, -1, 1)...)
        t11 = point!(builder, rot(1000, 1, 1)...)

        ## left face
        facetregion!(builder, 1)
        facet!(builder, b00, b01, t01, t00)

        ## right face
        facetregion!(builder, 2)
        facet!(builder, b10, b11, t11, t10)

        ## front face
        facetregion!(builder, 3)
        facet!(builder, b00, b10, t10, t00)

        ## back face
        facetregion!(builder, 4)
        facet!(builder, b01, b11, t11, t01)

        ## top face
        facetregion!(builder, 6)
        facet!(builder, t00, t10, t11, t01)

        ## bottom face
        facetregion!(builder, 5)
        facet!(builder, b00, b10, b11, b01)

        cellregion!(builder, 1)
        maxvolume!(builder, h)
        regionpoint!(builder, 0.5, 0.5, 0.5)

        return simplexgrid(builder)
    end

    @testset "CoupledDofsRestriction on diagonal cube grids (H1P2)" begin
        for θ in (0.0, π / 8, π / 4)
            let
                xgrid = create_diagonal_cube_grid(; h = 1.0, θ)
                FES = FESpace{H1P2{1, 3}}(xgrid)

                PD = ProblemDescription("Periodic Poisson on diagonal cube with 𝜃 = $θ")
                u = Unknown("u"; name = "u")
                assign_unknown!(PD, u)

                assign_operator!(PD, BilinearOperator([grad(u)]))

                assign_restriction!(PD, BoundaryDataRestriction(u; regions = [5, 6]))
                assign_restriction!(PD, CoupledDofsRestriction(u, 1, 2))
                assign_restriction!(PD, CoupledDofsRestriction(u, 3, 4))

                # assembly of the coupling matrix and solve must not crash
                sol, SC = solve(PD, FES; return_config = true)

                # coupling must actually have been assembled
                @test haskey(SC.statistics, :restriction_residuals)
            end
        end
    end

    return nothing
end
