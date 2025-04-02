module ProjectionOntoSLn

using LinearAlgebra: exp, tr, I, norm, Diagonal, diag, svd
using Distributions: Uniform
using ForwardDiff: ForwardDiff
using StaticArrays: @SArray

"""
    ProjectionResult

    Data type holding all relevant information of the projection procedure.
"""
struct ProjectionResult

    # resulting vector 𝑝 ∈ sl(𝑛)
    projection::AbstractVector

    # number of iterations of the projection method
    iterations::Unsigned
end


"""
    AbstractProjectionMethod

    Abstract type for the algorithms computing the projection of a given matrix.
"""
abstract type AbstractProjectionMethod end

"""
    project_onto_sln(method::AbstractProjectionMethod, a::AbstractVector; kwargs...)

    Interface for the projection of a vector in 𝑎 ∈ ℝⁿ onto the set sl(𝑛) := { 𝑝 ∈ ℝⁿ : ∏ᵢ 𝑝ᵢ = 1 }.
    implement an overload for a custom method.

    Input:
        method: method of the vector projection
        a:      input vector 𝑎 ∈ ℝⁿ

    Expected output:
        result: ProjectionResult of the projection (or error)

"""
function project_onto_sln(method::AbstractProjectionMethod, a::AbstractVector; kwargs...)
    error("no method of project_onto_sln implemented for $(nameof(method))")
end

"""
    project_onto_sln(method::AbstractProjectionMethod, A::AbstractMatrix; svd_method = svd, kwargs...)

    General projection method of a matrix 𝐴 ∈ ℝⁿˣⁿ onto SL(𝑛).
    The matrix is decomposed by the singular value decomposition. The singular are then projected and the result
    is assembled back into matrix form.

    Input:
        method:     method of the vector projection
        A:          𝑛×𝑛 matrix to be projected onto SL(𝑛)
        svd_method: method to compute the SVD `U, σ, V = svd_method(A)`, defaults to `LinearAlgebra.svd`

    Output: 𝑛×𝑛 matrix 𝑃 ∈ SL(𝑛), the projection of 𝐴.
"""
function project_onto_sln(method::AbstractProjectionMethod, A::AbstractMatrix; svd_method = svd, kwargs...)

    # compute singular value decomposition
    U, σ, V = svd_method(A)

    # compute projection in sl(𝑛)
    result = project_onto_sln(method, σ; kwargs...)

    # transform the result back to ℝⁿˣⁿ
    return U * Diagonal(result.projection) * V'

end

export project_onto_sln

include("utilities.jl")
export generateRandomMatrix
export testStationarity

include("root_finding.jl")
export root_finding

"""
    project_onto_sln(a::AbstractVector)

    Default vector projection with no specified method.
    In this case `RootFinding` is used since it is the most robust and stable method.
"""
function project_onto_sln(a::AbstractVector)
    return project_onto_sln(root_finding, a; maxIter = 999, tolerance = 1.0e-12, debug = false)
end


"""
    project_onto_sln(a::AbstractVecOrMat)

    Default vector projection with no specified method.
    In this case `RootFinding` is used since it is the most robust and stable method.
"""
function project_onto_sln(a::AbstractVecOrMat)
    return project_onto_sln(root_finding, a; maxIter = 999, tolerance = 1.0e-12, debug = false)
end


include("composite_step.jl")
export composite_step

include("constrained_newton.jl")
export constrained_newton

include("unconstrained_newton.jl")
export unconstrained_newton

end # module ProjectionOntoSLn
