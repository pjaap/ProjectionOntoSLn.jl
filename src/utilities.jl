#=
    This file contains utilities used by multiple algorithms
=#


"""
 Find the root of "f" between "x0" and "x1" (the order does not matter)
 xTol: relative stopping tolerance for |x1-x0|
 yTol: stopping tolerance for |y=f(x)|
 maxIter: maximal number of iterations before error
"""
function regulaFalsi(
        f::Function,
        x0::Number,
        x1::Number,
        maxIter::Integer = 100;
        xTol::AbstractFloat = 1.0e-12,
        yTol::AbstractFloat = 1.0e-15,
        debug::Bool = false,
    )

    y0 = f(x0)
    y1 = f(x1)

    for step in 1:maxIter

        if debug
            # this macro is for easy debugging
            @show step x0 x1 y0 y1
        end

        # use mean value
        x = (x0 + x1) / 2

        # compute value at x and stop check
        y = f(x)

        if abs(x1 - x0) ≤ xTol || abs(y) ≤ yTol
            return x, step
        end

        # regula falsi
        if sign(y * y0) == 1
            x0 = x
            y0 = y
        else
            x1 = x
            y1 = y
        end
    end

    # if we got here, we ran out of iterations
    error("maximal number of iterations exceeded in regulaFalsi: $(x1 - x0)")

    return nothing, maxIter
end


"""
    find_inital_value(a)

    Find a robust and "good" initial p ∈ sl(d) for a given non-negative a ∈ C(n)
"""
function find_inital_value(a)

    sum_a = sum(a)
    n = length(a)
    if sum_a < n
        a = a .+ 1.0 .- sum_a / n # do not overwrite `a`
    end

    # safety against aᵢ=0
    a = a .+ 1.0e-15

    p = prod(a)^(-1 / n) * a

    if p[begin] > 1
        γ = log(a[begin]) / log(p[begin])
        p = p .^ γ
    end

    return p
end


# if P is closest to A then A-P is multiple of P^{-T}
function testStationarity(A::AbstractMatrix, P::AbstractMatrix)
    PiT = inv(P)'

    # this should be ≈1
    return abs((PiT ⋅ (A - P)) / (norm(PiT) * norm(A - P))) ≈ 1
end

# variant for the reduced vector case
function testStationarity(a::AbstractVector, p::AbstractVector)
    p_inv = 1 ./ p

    # this should be ≈1
    res = abs((p_inv' * (a - p)) / (norm(p_inv) * norm(a - p)))
    if res ≈ 1
        return true
    else
        @show abs((p_inv' * (a - p)) / (norm(p_inv) * norm(a - p)))
        return false

    end
end


# deviator of a matrix (overwrites A)
function deviator!(A::AbstractMatrix)
    n = size(A, 1)
    A .-= 1 / n * tr(A) * I(n)

    return nothing
end

# return deviator of a matrix
function deviator(A::AbstractMatrix{<:AbstractFloat})
    A_copy = copy(A)
    deviator!(A_copy)
    return A_copy
end

"""
    Generate a 'n X n' random matrix with determinant in the range determinant_range

    diag_only: return a diagonal matrix
"""
function generateRandomMatrix(n; det_range, diag_only, kwargs...)

    # Bates distribution with √n scaling
    rand_range = @. log(det_range) / √n

    # create a random matrix
    if diag_only
        B = Diagonal(rand(Uniform(rand_range...), n))
    else
        B = rand(Uniform(rand_range...), n, n)
    end

    return exp(B)
end


"""
    report_non_convergence(method, a, filename="non_convergence.log")

    Append to `filename` that a method did not converge for a given `a`
"""
function report_non_convergence(method, a, filename = "non_convergence.log")
    return open(filename, "a") do file
        write(file, "Method $method did not converge with dim = $(length(a)) and a = $a\n")
    end
end
