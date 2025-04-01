# ProjectionOntoSLn.jl
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

This project was created to solve the following optimization problem:

For a given matrix $A \in \mathbb R^{n \times n}$, compute $P \in \textrm{SL}(n)$ that minimizes the Frobenius distance
```math
\Vert A - P \Vert_F.
```
The Lie group $\textrm{SL}(n)$ is given by matrices with determinant one:

```math
\textrm{SL}(n) \coloneqq \{ P \in \mathbb R^{n\times n}:\ \det P = 1\}.

```

# Citation

The methods used in this project are explained and analyzed in our preprint paper [How to project onto SL(n)](https://arxiv.org/abs/2501.19310) with the DOI [
https://doi.org/10.48550/arXiv.2501.19310](
https://doi.org/10.48550/arXiv.2501.19310).

Please cite our work if you use the code:

```
@misc{jaap2025projectsln,
      title={How to project onto SL($n$)},
      author={Patrick Jaap and Oliver Sander},
      year={2025},
      eprint={2501.19310},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2501.19310},
}
```

# Installation

To include this project into our own Julia project, you simply run
```
pkg> add https://github.com/pjaap/ProjectionOntoSLn.jl
```

# Usage

The main function exported by this package is `project_onto_sln`.
This function can be called with a square matrix.

```
using ProjectionOntoSLn
using LinearAlgebra

A = rand(3,3)

P = project_onto_sln(A)

@show det(P) # will be close to 1
```

## Advanced Usage

Currently, there are four methods implemented to compute the projection:
- Direct root finding: `root_finding` (this is the default if no method is specified, as seen above)
- Iterative composite step: `composite_step`
- Newton with scalar constraint: `constrained_newton`
- An unconstrained Newton method: `unconstrained_newton`

The methods are explained in depth in the our [paper](https://arxiv.org/abs/2501.19310).

To use, e.g., the constrained Newton method, you can run

```
P = project_onto_sln(constrained_newton, A; maxIter=100, tolerance=1e-12, debug=false)
```


## Vector Variant

Internally, the matrix methods compute a singular value decomposition $U \Sigma V^T = A$ with a diagonal matrix $\Sigma = \textrm{diag}(\sigma)$.
In the [paper](https://arxiv.org/abs/2501) we show that the projection is done efficiently in the singular value space.

Therefore, you can call `project_onto_sln` also with a vector `a`

```
a = rand(3)

result = project_onto_sln(a)
p = result.projection

@show prod(p) # will be close to 1
```

Note, in the vector case, a solution object of type `ProjectionResult` is returned.
This also contains the number of iterations.

# Contributing

Feel free to open PRs in order to contribute to the project :)
