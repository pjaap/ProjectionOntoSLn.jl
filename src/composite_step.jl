struct CompositeStep <: AbstractProjectionMethod end
const composite_step = CompositeStep()


function project_onto_sln(::CompositeStep, a::AbstractVector; kwargs...)
    p0 = haskey(kwargs, :p_init) ? kwargs[:p_init] : find_inital_value(a)
    return compositeStepLoop(a, p0, kwargs[:maxIter], kwargs[:tolerance], kwargs[:debug])
end


function compositeStep(a::AbstractArray, pk::AbstractArray, tolerance)

    # this is the search direction
    pki = 1 ./ pk

    # objective function for root finding
    f = t -> (prod(a + t * pki) - 1)

    # bounds
    if prod(a) < 1
        x0, x1 = 0.0, 1.0 - prod(a)^(1 / length(a))
    else
        x0, x1 = maximum(-1 * (a .* pk)), 0.0
    end

    t, _ = regulaFalsi(f, x0, x1, 1000; xTol = tolerance)

    return a + t * pki
end

function compositeStepLoop(a::AbstractVector, p0::AbstractVector, maxIter::Integer, tol::AbstractFloat, debug::Bool)

    pk = p0

    if debug
        # remember all iterates
        pk_matrix = zeros(maxIter + 1, length(pk))
        pk_matrix[1, :] = pk
    end

    p_init = copy(pk)

    diff = 0

    for k in 1:maxIter

        pk1 = compositeStep(a, pk, tol)

        if debug
            pk_matrix[k + 1, :] = pk
            @show k norm(pk1 - pk) / norm(pk1)
        end

        diff = norm(pk1 - pk)

        if diff ≤ tol * norm(pk1) || diff < 1.0e-15
            return ProjectionResult(pk1, k, Dict(:distance_to_p0 => norm(p_init - pk1)))
        end

        pk = pk1
    end

    @show diff

    @show pk_matrix

    report_non_convergence("composite step", a)

    error("maximal number of iterations exceeded")

    return ProjectionResult(p0, maxIter, [])
end
