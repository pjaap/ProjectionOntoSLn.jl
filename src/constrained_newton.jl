#=
  This implements the constrained Newton Algorithm from the paper
=#

struct ConstrainedNewton <: AbstractProjectionMethod end
const constrained_newton = ConstrainedNewton()

function project_onto_sln(::ConstrainedNewton, a::AbstractVector; kwargs...)

    # compute initial value
    p_init = find_inital_value(a)
    p_log = log.(p_init)

    p0 = haskey(kwargs, :p_init) ? kwargs[:p_init] : p_log
    λ0 = haskey(kwargs, :λ_init) ? kwargs[:λ_init] : 0.0

    return constrainedNewtonLoop(a, p0, λ0, kwargs[:maxIter], kwargs[:tolerance], kwargs[:debug])
end


function constrainedNewtonUpdate(a, pk)

    # cache this result
    exp_pk = ℯ .^ pk

    rhs = vcat(-exp_pk .^ 2 + exp_pk .* a, 0.0)

    D = Diagonal(2exp_pk .^ 2 - exp_pk .* a)

    # this is lazy and fast
    D_inv = inv(D)

    ρ = -1 / tr(D_inv)
    u = vcat(diag(D_inv), -1.0)

    result = ρ * (u'rhs) * u

    M1 = Diagonal(vcat(diag(D_inv), 0))
    result += M1 * rhs

    # split p and λ
    return result[1:(end - 1)], result[end]
end

function constrainedNewtonLoop(a, p_k, λ_k, maximal_iterations, tolerance, debug)

    # keep that
    p_init = copy(p_k)

    if debug
        @show p_init exp.(p_init)
    end

    for i in 1:maximal_iterations
        dp_k, dλ_k = constrainedNewtonUpdate(a, p_k)

        p_k += dp_k
        λ_k += dλ_k

        if debug
            @show i norm(dp_k, Inf) p_k λ_k exp.(p_k)
        end

        norm_diff = norm(dp_k, Inf)
        norm_iter = norm(p_k, Inf)

        if norm_diff ≤ tolerance * norm_iter || norm_diff ≤ 1.0e-15
            return ProjectionResult(exp.(p_k), i)
        end


    end

    report_non_convergence("constrained newton", a)

    error("constrainedNewtonLoop did not converge")

    return ProjectionResult(exp.(p_k), maximal_iterations)
end
