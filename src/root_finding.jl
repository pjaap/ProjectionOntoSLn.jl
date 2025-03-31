struct RootFinding <: AbstractProjectionMethod end
const root_finding = RootFinding()

function project_onto_sln(::RootFinding, a::AbstractVector; kwargs...)
    return rootFinding(a; kwargs...)
end

function pPlus(a, λ)
    return a / 2 + sqrt.(a .^ 2 / 4 .- λ)
end

function pMinus(a, λ)
    p = zeros(typeof(λ), length(a))
    p[1:(end - 1)] = @views a[begin:(end - 1)] / 2 + sqrt.(a[begin:(end - 1)] .^ 2 / 4 .- λ)
    p[end] = a[end] / 2 - sqrt(a[end]^2 / 4 - λ)

    return p
end


# define the solution path and the search interval
function makePathAndRange(a)
    if prod(a) ≤ 1
        return λ -> pPlus(a, λ), -1.0, 0.0
    else
        function P(λ)
            if λ ≤ a[end]^2 / 4
                # P^- in the paper
                return pMinus(a, λ)
            else
                # P^+ in the paper
                return pPlus(a, a[end]^2 / 2 - λ)
            end
        end
        return P, 0.0, a[end]^2 / 2
    end
end


# this implements the root finding algorithm
function rootFinding(a; tolerance, debug, maxIter, kwargs...)

    # step 2: regula falsi
    P, x0, x1 = makePathAndRange(a)

    f(λ) = prod(P(λ)) - 1

    λ, iterations = regulaFalsi(f, x0, x1, maxIter; xTol = tolerance, debug)

    # step 3: construct p on the path
    return ProjectionResult(P(λ), iterations, Dict(:distance_to_p0 => norm(a - P(λ))))
end
