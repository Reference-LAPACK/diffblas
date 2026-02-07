---
title: 'diffblas: algorithmically differentiated BLAS routines'
tags:
  - automatic differentiation
  - linear algebra
authors:
  - name: Alexis Montoison^[corresponding author]
    orcid: 0000-0002-3403-5450
    affiliation: 1
  - name: Sri Hari Krishna Narayanan
    orcid: 0000-0003-0388-5943
    affiliation: 1
  - name: Jean-Luc Bouchot
    orcid: 0000-0003-4523-3986
    affiliation: 2
affiliations:
 - name: Argonne National Laboratory, Lemont, IL, USA.
   index: 1
 - name: Inria de Saclay, Palaiseau, France.
   index: 2
date: 6 February 2026
bibliography: paper.bib

---

# Summary

`diffblas`is a library that  provides (algorithmically) differentiated BLAS routines from their reference implementation in LAPACK using the automatic differentiation tool Tapenade [@tapenade] in four modes: forward (`_d`), vector forward (`_dv`), reverse (`_b`), and vector reverse (`_bv`).

# Statement of need

# Acknowledgements

This work was supported in part by the Applied Mathematics activity within the U.S. Department of Energy, Office of Science, Office
of Advanced Scientific Computing Research Applied Mathematics, and Office of Nuclear Physics SciDAC program under Contract No. DE-AC02-06CH11357. This work was supported in part by NSF CSSI grant 2104068.

# References
