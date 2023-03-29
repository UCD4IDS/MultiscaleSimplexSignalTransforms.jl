# MultiscaleSimplexSignalTransforms.jl

## COPYRIGHT

Copyright 2023 The Regents of the University of California

Implemented by Eugene Shvarts

## SETUP

`MultiscaleSimplexSignalTransforms` is not yet in the `General` registry, so either add via the `Tetrapods` registry

```
] registry add https://github.com/UCD4IDS/TetrapodsRegistry
] add MultiscaleSimplexSignalTransforms
```

or add manually by URL

```
] add https://github.com/UCD4IDS/MultiscaleSimplexSignalTransforms.jl
```

## USAGE

- At top level are the graph basis dictionaries `kGHWT` and `kHGLET`, implementing the abstract type `MSST`, or `MultiscaleSimplicialTransform`.
- These dictionaries themselves rely on a method of partitioning simplicial complexes, which implements the abstract type `SCPartition`.
The implementations are `SubmatrixPartition` (the default), and `FullPartition`.
- Each possesses an array of configuration options, `Representation`, `SubRepresentation`, `Basis`, `PartitionInput`, `PartitionMethod`, `EigenMethod`, which are extensible and allow for configurable, stackable, and repeatable experiments.
- The fundamental adjacency data structures are `ZeroRegion` and `KRegion`, which implement the abstract type `Region`, and the fundamental spectral representation for these is `k_laplacian`.
- The structure of a simplicial complex is stored in a `SimplexTree`, and generally speaking, when a function takes a `SimplexTree`, it will happily accept some `g::AbstractGraph` instead by passing in `cliquecomplex(g, k)` for an appropriate `k`.

If you have some `g::AbstractGraph`, then the easiest way to get started with analyzing, say, signals on the triangles of `g` with all defaults set is to construct the basis dictionary `basis = kGHWT(KRegion(g, 2))`. Then you can investigate the dictionary vector at level `j`, location `k`, tag `l` with ordinary indexing (i.e., `basis[j,k,l]`), and you can obtain a dictionary of expansion coefficients for some triangle signal `s` with `analyze(basis, s)`.
