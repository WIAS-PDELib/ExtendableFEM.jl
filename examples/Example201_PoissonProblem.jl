#= 

# 201 : Poisson-Problem
([source code](SOURCE_URL))

This example computes the solution ``u`` of the two-dimensional Poisson problem
```math
\begin{aligned}
-\Delta u & = f \quad \text{in } \Omega
\end{aligned}
```
with right-hand side ``f(x,y) \equiv xy`` and homogeneous Dirichlet boundary conditions
on the unit square domain ``\Omega`` on a given grid.

=#

module Example201_PoissonProblem

using ExtendableFEM
using ExtendableFEMBase
using ExtendableGrids
using GridVisualize

function f!(fval, qpinfo)
	fval[1] = 1#qpinfo.x[1] * qpinfo.x[2]
end

function main(; μ = 1.0, nrefs = 4, order = 2, Plotter = nothing, kwargs...)

	## problem description
	PD = ProblemDescription()
	u = Unknown("u"; name = "potential")
	assign_unknown!(PD, u)
	assign_operator!(PD, BilinearOperator([grad(u)]; factor = μ, kwargs...))
	assign_operator!(PD, LinearOperator(f!, [id(u)]; kwargs...))
	assign_operator!(PD, HomogeneousBoundaryData(u; regions = 1:4))

	## discretize
	xgrid = uniform_refine(grid_unitsquare(Triangle2D), nrefs)
	FES = FESpace{H1Pk{1, 2, order}}(xgrid)

	## solve
	sol = solve(PD, FES; kwargs...)

	## plot
	p = GridVisualizer(; Plotter = Plotter, layout = (1, 2), clear = true, size = (1000, 500))
	scalarplot!(p[1, 1], xgrid, nodevalues(sol[u])[:]; title = "u")
	scalarplot!(p[1, 2], xgrid, view(nodevalues(sol[u], Gradient; abs = true), 1, :), title = "∇u (abs + quiver)")
	vectorplot!(p[1, 2], xgrid, eval_func(PointEvaluator([grad(u)], sol)), vscale = 0.8, clear = false)
end

end # module