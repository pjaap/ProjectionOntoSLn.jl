#=
  This implements the conforming descent algorithm
  for hyperbolic coordinates from the paper
=#


struct UnconstrainedNewton <: AbstractProjectionMethod end
const unconstrained_newton = UnconstrainedNewton()

function project_onto_sln(::UnconstrainedNewton, a::AbstractVector; kwargs...)
    return unconstrainedNewtonLoop(a; maxIter = kwargs[:maxIter], tolerance = kwargs[:tolerance], debug = kwargs[:debug])
end

# this implements the matrix B with size (n × n) with entries of type T
function B(n, ::Type{T}) where {T}

    # this can be dangerous, but we believe that the user uses a large type for 'n'
    if n > typemax(T)
        error("The matrix B cannot be constructed with")
    end

    if n < 1
        error("The matrix B must have positive size arguments")
    end

    # allocate a matrix
    M = zeros(T, n, n - 1)

    for i in 1:n
        for j in 1:(i - 1)
            M[i, j] = -j
        end
        for j in i:(n - 1)
            M[i, j] = n - j
        end
    end

    return M
end

# we assume that we do not need large values for 'n'
function B(n)
    return B(n, Int8)
end


# hyperbolic distance function to the point "a"
function dist_hyp(a, mu, Bn)
    return norm(a .- exp.(Bn * mu))^2
end


function unconstrainedNewtonLoop(a; maxIter, tolerance, debug)

    n = length(a)

    # initial guess is the transformed value of "a"
    p_init = find_inital_value(a)

    # transform to log coordinates
    p_log = log.(p_init)

    # apply B⁻¹
    mu_k = zeros(n - 1)
    for i in 1:(n - 1)
        mu_k[i] = (p_log[i] - p_log[i + 1]) / n
    end

    # keep that
    mu_init = copy(mu_k)

    Bn = B(n)

    # some convenient functions
    f(mu) = dist_hyp(a, mu, Bn)
    df(mu) = ForwardDiff.gradient(f, mu)
    ddf(mu) = ForwardDiff.hessian(f, mu)

    for i in 1:maxIter

        hess = ddf(mu_k)
        grad = df(mu_k)

        # stop if grad = 0
        if norm(grad, Inf) ≤ 1.0e-15
            return ProjectionResult(exp.(Bn * mu_k), i, Dict(:distance_to_p0 => norm(exp.(Bn * mu_k) - exp.(Bn * mu_init))))
        end

        # we do not put the minus sign here: we subtract c_k at all occurrences
        c_k = hess \ grad

        mu_k .-= c_k

        if debug
            @show i norm(c_k) mu_k exp.(Bn * mu_k)
        end

        norm_diff = norm(c_k, Inf)
        norm_iter = norm(mu_k, Inf)

        if norm_diff ≤ tolerance * norm_iter || norm_diff ≤ 1.0e-15
            return ProjectionResult(exp.(Bn * mu_k), i, Dict(:distance_to_p0 => norm(exp.(Bn * mu_k) - exp.(Bn * mu_init))))
        end

    end

    report_non_convergence("conforming descent", a)

    error("conformingDescentLoop did not converge")

    # assemble P = exp( B * mu )
    return ProjectionResult(exp.(Bn * mu_k), maxIter, [])
end
