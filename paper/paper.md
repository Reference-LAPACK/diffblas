---
title: 'diffblas: algorithmically differentiated BLAS routines'
tags:
  - BLAS
  - automatic differentiation
  - numerical linear algebra
  - scientific computing
authors:
  - name: Sri Hari Krishna Narayanan^[corresponding author]
    orcid: 0000-0003-0388-5943
    affiliation: 1
  - name: Alexis Montoison
    orcid: 0000-0002-3403-5450
    affiliation: 1
  - name: Jean-Luc Bouchot
    orcid: 0000-0003-4523-3986
    affiliation: 2
affiliations:
 - name: Mathematics and Computer Science Division, Argonne National Laboratory, USA
   index: 1
 - name: Inria de Saclay, Palaiseau, France
   index: 2
date: 22 February 2026
bibliography: paper.bib

---

# Summary

`diffblas` is a library that provides BLAS routines algorithmically differentiated using algorithmic differentiation principles [@griewank2008] from their reference implementations in [LAPACK](https://github.com/Reference-LAPACK/lapack) on GitHub using the automatic differentiation tool Tapenade [@tapenade].
It supports four modes: forward (`_d`), vector forward (`_dv`), reverse (`_b`), and vector reverse (`_bv`).

In addition to differentiating the standard Fortran-style `BLAS` interface, `diffblas` also provides differentiated `CBLAS` routines, facilitating interoperability with C and other languages.
Its API mirrors BLAS / CBLAS, with additional arguments specifying differentiation variables, making integration into existing workflows straightforward.

`diffblas` calls the underlying standard `BLAS `implementation, and is agnostic to the backend (OpenBLAS, BLIS, MKL, Apple Accelerate), ensuring both performance and portability.
Precompiled artifacts are available on GitHub to simplify installation and usage across different programming environments.

By providing efficient and accurate derivatives of linear algebra operations, `diffblas` facilitates gradient-based optimization, sensitivity analysis, and derivative-based workflows in scientific computing.

# Statement of need

Linear algebra routines such as those in LAPACK are widely used in scientific computing, optimization, and machine learning. However, they do not provide derivatives, which are often required for gradient-based algorithms.

Existing approaches rely on hand-coded derivatives or generic automatic differentiation applied to high-level code, which can be inefficient or error-prone [@jonasson2020].

`diffblas` addresses this gap by providing algorithmically differentiated BLAS routines directly from reference LAPACK implementations and following relevant differntiation rules [@giles2008].

# State of the field

Automatic source-to-source differentiation tools, such as Tapenade [@tapenade], ADOL-C [@ADOLC], or TAF [@TAF], provide general mechanisms to compute derivatives of code.
However, naively applying AD to low-level BLAS or LAPACK routines can be inefficient due to loop structures, memory layout, and caching constraints [@jonasson2020].
`diffblas` addresses these challenges by generating differentiated routines directly from reference LAPACK implementations, combining the reliability and performance of LAPACK with the efficiency of algorithmic differentiation.

# Research impact statement

This work was inspired in part by a need to differentiate a Fortran code [@HFBTHO] that uses BLAS and LAPACK routines, and to use the differentiated application for gradient-based optimization.

Providing both the standard and CBLAS interfaces ensures that diffblas can be adopted across different programming environments, facilitating derivative computations in diverse scientific computing projects.
Precompiled artifacts on GitHub further simplify integration, enabling rapid deployment in multiple languages and scientific computing projects.

# Acknowledgements

This work was supported in part by the Applied Mathematics activity within the U.S. Department of Energy, Office of Science, Office
of Advanced Scientific Computing Research Applied Mathematics, and Office of Nuclear Physics SciDAC program under Contract No. DE-AC02-06CH11357. This work was supported in part by NSF CSSI grant 2104068.

# AI usage disclosure

Generative AI was used to ...
ChatGPT was used to check spelling, grammar, and clarity of the English text in this paper.

# References
