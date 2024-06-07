#=

# 203 : Poisson-Problem DG
([source code](@__SOURCE_URL__))

This example computes the solution ``u`` of the two-dimensional Poisson problem
```math
\begin{aligned}
-\Delta u & = f \quad \text{in } \Omega
\end{aligned}
```
with right-hand side ``f`` and inhomogeneous Dirichlet boundary conditions
chosen such that ``u(x,y) = x^3 - 3xy^2``.
This time the problem is solved on a given grid via the discontinuous Galerkin method.

The computed solution looks like this:

![](example203.svg)
=#

module Example203_PoissonProblemDG

using ExtendableFEM
using ExtendableGrids
using Symbolics
using Test #hide

## exact data for problem by Symbolics
function prepare_data(; μ = 1)

	@variables x y

	## exact solution
	u = x^3 - 3 * x * y^2

	## gradient
	∇u = Symbolics.gradient(u, [x, y])

	## Laplacian
	Δu = Symbolics.gradient(∇u[1], [x]) + Symbolics.gradient(∇u[2], [y])

	## right-hand side
	f = -μ * Δu[1]

	## build functions
	u_eval = build_function(u, x, y, expression = Val{false})
	∇u_eval = build_function(∇u, x, y, expression = Val{false})
	f_eval = build_function(f, x, y, expression = Val{false})

	return f_eval, u_eval, ∇u_eval[2]
end

function main(; dg = true, μ = 1.0, τ = 10.0, nrefs = 4, order = 2, bonus_quadorder = 2, Plotter = nothing, kwargs...)

	## prepare problem data
	f_eval, u_eval, ∇u_eval = prepare_data(; μ = μ)
	rhs!(result, qpinfo) = (result[1] = f_eval(qpinfo.x[1], qpinfo.x[2]))
	exact_u!(result, qpinfo) = (result[1] = u_eval(qpinfo.x[1], qpinfo.x[2]))
	exact_∇u!(result, qpinfo) = (∇u_eval(result, qpinfo.x[1], qpinfo.x[2]))

	## problem description
	PD = ProblemDescription()
	u = Unknown("u"; name = "potential")
	assign_unknown!(PD, u)
	assign_operator!(PD, BilinearOperator([grad(u)]; factor = μ, kwargs...))
	assign_operator!(PD, LinearOperator(rhs!, [id(u)]; bonus_quadorder = bonus_quadorder, kwargs...))
	assign_operator!(PD, InterpolateBoundaryData(u, exact_u!; bonus_quadorder = bonus_quadorder, regions = 1:4))

	## discretize
	xgrid = uniform_refine(grid_unitsquare(Triangle2D), nrefs)
	FES = FESpace{H1Pk{1, 2, order}}(xgrid; broken = dg)

	## add DG terms
	assign_operator!(PD, BilinearOperatorDG(dg_kernel(xgrid), [jump(id(u))], [average(grad(u))]; entities = ON_FACES, factor = -μ, kwargs...))
	assign_operator!(PD, BilinearOperatorDG(dg_kernelT(xgrid), [average(grad(u))], [jump(id(u))]; entities = ON_FACES, factor = -μ, kwargs...))
	assign_operator!(PD, LinearOperatorDG(dg_kernel_bnd(xgrid, exact_u!), [average(grad(u))]; entities = ON_BFACES, factor = -μ, bonus_quadorder = bonus_quadorder, kwargs...))
	assign_operator!(PD, BilinearOperatorDG(dg_kernel2(xgrid), [jump(id(u))]; entities = ON_FACES, factor = τ, kwargs...))
	assign_operator!(PD, LinearOperatorDG(dg_kernel2_bnd(xgrid, exact_u!), [id(u)]; entities = ON_BFACES, regions = 1:4, factor = τ, bonus_quadorder = bonus_quadorder, kwargs...))

	## solve
	sol = solve(PD, FES; kwargs...)

	## prepare error calculation
	function exact_error!(result, u, qpinfo)
		exact_u!(result, qpinfo)
		exact_∇u!(view(result, 2:3), qpinfo)
		result .-= u
		result .= result .^ 2
	end
	function dgjumps!(result, u, qpinfo)
		result .= u[1]^2/qpinfo.volume
	end
	ErrorIntegratorExact = ItemIntegrator(exact_error!, [id(u), grad(u)]; quadorder = 2 * order, params = [μ], kwargs...)
	DGJumpsIntegrator = ItemIntegratorDG(dgjumps!, [jump(id(u))]; entities = ON_IFACES, kwargs...)

	## calculate error
	error = evaluate(ErrorIntegratorExact, sol)
	dgjumps = sqrt(sum(evaluate(DGJumpsIntegrator, sol)))
	L2error = sqrt(sum(view(error, 1, :)))
	H1error = sqrt(sum(view(error, 2, :)) + sum(view(error, 3, :)))
	@info "L2 error = $L2error"
	@info "H1 error = $H1error"
	@info "dgjumps = $dgjumps"

	## plot
	plt = plot([id(u), grad(u)], sol; Plotter = Plotter)

	return L2error, plt
end

function dg_kernel(xgrid)
	facenormals = xgrid[FaceNormals]
	facecells = xgrid[FaceCells]
	facevolumes = xgrid[FaceVolumes]
	function closure(result, input, qpinfo)
		result[1] = (input[1] * facenormals[1, qpinfo.item] + input[2] * facenormals[2, qpinfo.item])
	end
end

function dg_kernelT(xgrid)
	facenormals = xgrid[FaceNormals]
	facecells = xgrid[FaceCells]
	facevolumes = xgrid[FaceVolumes]
	function closure(result, input, qpinfo)
		result[1:2] = input[1] .* view(facenormals, :, qpinfo.item)
	end
end

function dg_kernel_bnd(xgrid, uDb! = nothing)
	facenormals = xgrid[FaceNormals]
	facecells = xgrid[FaceCells]
	facevolumes = xgrid[FaceVolumes]
	function closure(result, qpinfo)
		uDb!(result, qpinfo)
		result[1:2] = result[1] .* view(facenormals, :, qpinfo.item)
	end
end
function dg_kernel2(xgrid)
	function closure(result, input, qpinfo)
		result .= input / qpinfo.volume
	end
end
function dg_kernel2_bnd(xgrid, uDb! = nothing)
	function closure(result, qpinfo)
		uDb!(result, qpinfo)
		result /= qpinfo.volume
	end
end

generateplots = default_generateplots(Example203_PoissonProblemDG, "example203.svg") #hide
function runtests() #hide
	L2error, ~ = main(; μ = 0.25, nrefs = 2, order = 2) #hide	
	@test L2error ≈ 0.00025771757957310844 #hide
end #hide
end # module
