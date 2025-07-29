# SeuratIntegrate (development version)

* Fix overestimation of number of zeros in sparse matrices: remaining 0s in the
matrix values are discarded before counting non-zero values
(`src/n_zeros_matrix.cpp`)

* Fix `choose_matrix_format()` (`utils.R`)

* Created a first set of unit tests covering `utils.R` (testthat)

* Add method to preserve a `Graph` object's attributes when using `as.Graph` on
a `Graph` object

* Fix deprecation warning in `GetNeighborsPerBatch` (`Graph` method)

* Minor fixes in docstrings


# SeuratIntegrate 0.4.1

* Revised score rescaling with a new option enabling min-max rescaling of ranks
(default) rather than scores directly (as in
[Luecken *et al.*, 2021](https://doi.org/10.1038/s41592-021-01336-8))

* Speed up Dijkstra's algorithm-like used in `ExpandNeighbours` with a new c++
implementation

* The most suited matrix format is automatically chosen for corrected counts
output by integration methods (should be dense matrix most of the time)

* Add support for scGraph metric ([Wang *et al.*, 2024](https://doi.org/10.1101/2024.04.02.587824))

* Improved speed of `CreateIntegrationGroups` for non-SCT assay in unambiguous cases


# SeuratIntegrate 0.4.0

* Initial public release
