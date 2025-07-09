# CallbackOperator

`CallbackOperator` provides a flexible interface for injecting custom assembly logic into the finite element workflow. It delegates the assembly of the system matrix and right-hand side to a user-supplied callback function, allowing for advanced or nonstandard modifications that are difficult to express with standard operators.

## Constructor

```@autodocs
Modules = [ExtendableFEM]
Pages = ["common_operators/callback_operator.jl"]
Order   = [:type, :function]
```

## Example

```julia
function my_callback!(A, b, args; assemble_matrix=true, assemble_rhs=true, time=0, kwargs...)
    # Example: add a constant to the diagonal
    if assemble_matrix && A !== nothing
        for i in 1:min(size(A)...)
            A[i, i] += 1.0
        end
    end
    if assemble_rhs && b !== nothing
        b .+= 2.0
    end
end
op = CallbackOperator(my_callback!; u_args=[1], name="CustomOp")
```

## When to Use

- Custom or experimental PDE terms
- Coupling to external codes or data
- Advanced boundary or interface conditions
- Prototyping new assembly strategies

See also: [Example265](https://wias-pdelib.github.io/ExtendableFEM.jl/stable/examples/) for a practical use case.
