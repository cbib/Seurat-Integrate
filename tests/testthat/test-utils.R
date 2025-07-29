test_that("matrix is converted to dgCMatrix", {
  m <- matrix(nrow = 0, ncol=0)
  expect_equal(as.dgcmatrix(m), as(m, "dgCMatrix"))

  m <- matrix(rnorm(500, 0, 20), nrow = 20, ncol = 25)
  expect_equal(as.dgcmatrix(m), as(m, "dgCMatrix"))

  m[sample.int(500, 400)] <- 0
  expect_equal(as.dgcmatrix(m), as(m, "dgCMatrix"))
})


test_that("row-wise matrix sorting works", {
  m <- matrix(nrow = 0, ncol=0)
  expect_equal(rowSort(m), m)

  m <- matrix(c(1,2,3,4,5,
                1,2,2,4,5), nrow=2, byrow = T)
  expect_equal(rowSort(m), m)
  expect_equal(rowSort(m, decreasing = TRUE), m[,5:1])

  expect_equal(rowSort(as.dgcmatrix(m)), m)
  expect_equal(rowSort(as.dgcmatrix(m), decreasing = TRUE), m[,5:1])

  m <- matrix(rnorm(500, 0, 20), nrow = 25, ncol = 20)
  expect_equal(rowSort(m), naive_row_sort(m))
  expect_equal(rowSort(m, decreasing = TRUE), naive_row_sort(m, decreasing = TRUE))

  m_ <- t(sapply(c(rep(F, 20), rep(T, 5)), sample.int, n = ncol(m), size = ncol(m)))
  expect_equal(rowSort(m_, m), naive_row_sort(m_, m))
  expect_equal(rowSort(m_, m, decreasing = TRUE), naive_row_sort(m_, m, decreasing = TRUE))

  # expected <- t(sapply(1:nrow(m), function(i) m_[i,order(m[i,])]))
  # expect_equal(rowSort(m_, m), expected)
  # expected <- t(sapply(1:nrow(m), function(i) m_[i,order(m[i,], decreasing = TRUE)]))
  # expect_equal(rowSort(m_, m, decreasing = TRUE), expected)
})

test_that("col-wise matrix sorting works", {
  m <- matrix(nrow = 0, ncol=0)
  expect_equal(colSort(m), m)

  m <- matrix(c(1,2,3,4,5,
                1,2,2,4,5), ncol=2, byrow = F)
  expect_equal(colSort(m), m)
  expect_equal(colSort(m, decreasing = TRUE), m[5:1,])

  expect_equal(colSort(as.dgcmatrix(m)), m)
  expect_equal(colSort(as.dgcmatrix(m), decreasing = TRUE), m[5:1,])

  m <- matrix(rnorm(500, 0, 20), nrow = 20, ncol = 25)
  expect_equal(colSort(m), naive_col_sort(m))
  expect_equal(colSort(m, decreasing = TRUE), naive_col_sort(m, decreasing = TRUE))

  m_ <- sapply(c(rep(F, 20), rep(T, 5)), sample.int, n = nrow(m), size = nrow(m))
  expect_equal(colSort(m_, m), naive_col_sort(m_, m))
  expect_equal(colSort(m_, m, decreasing = TRUE), naive_col_sort(m_, m, decreasing = TRUE))
  # expected <- sapply(1:ncol(m), function(i) m_[order(m[,i]), i])
  # expect_equal(colSort(m_, m), expected)
  # expected <- sapply(1:ncol(m), function(i) m_[order(m[,i], decreasing = TRUE), i])
  # expect_equal(colSort(m_, m, decreasing = TRUE), expected)
})

test_that("matrix sorting status testing works", {
  m <- matrix(nrow = 0, ncol=0)
  expect_error(rowSorted(m), "subscript out of bounds$")
  expect_error(colSorted(m), "subscript out of bounds$")

  m <- matrix(c(1,2,3,4,5,
                1,1,3,5,5,
                1,5,4,2,3), nrow=3, byrow = T)
  expect_equal(rowSorted(m), c(T, T, F))
  expect_equal(colSorted(m), c(T, F, T, F, F))


  m <- matrix(c( 6,   5,   4,   3,   0,
                 2,-.99,-.99,-.99,-.99,
                -8,  -5,   2,   3,  -9), nrow=3, byrow = T)
  expect_equal(rowSorted(m, decreasing = TRUE), c(T, T, F))
  expect_equal(colSorted(m, decreasing = TRUE), c(T, T, F, F, T))
})

test_that("matrix indices conversion works", {
  expect_equal(rowcol2idx(0,0,0), 0)
  expect_equal(rowcol2idx(0,3,0), 0)
  expect_equal(rowcol2idx(3,0,0), 3) # should be 0?

  expect_nan(idx2row(0, 0))
  expect_nan(idx2row(3, 0))
  expect_nan(idx2row(-3, 0))

  expect_inf(idx2col(0,0))
  expect_inf(idx2col(3,0))
  expect_inf(idx2col(-3,0))


  nrow <- sample(50:200, 1)
  ncol <- sample(100:400, 1)

  m <- matrix(1:(nrow*ncol), nrow = nrow, ncol = ncol)
  expect_equal(outer(1:nrow, 1:ncol, rowcol2idx, height = nrow), m)
  expect_equal(idx2col(m, nrow), matrix(rep(1:ncol, each = nrow), nrow = nrow))
  expect_equal(idx2row(m, nrow), matrix(rep(1:nrow, each = ncol), nrow = nrow, byrow = TRUE))

  m <- matrix(1:(nrow*ncol), nrow = ncol, ncol = nrow)
  expect_equal(outer(1:ncol, 1:nrow, rowcol2idx, height = ncol), m)
  expect_equal(idx2col(m, ncol), matrix(rep(1:nrow, each = ncol), nrow = ncol))
  expect_equal(idx2row(m, ncol), matrix(rep(1:ncol, each = nrow), nrow = ncol, byrow = TRUE))
})

test_that("matrix symmetry works", {
  m <- as.dgcmatrix(matrix(nrow = 0, ncol = 0))
  expect_squared(SymmetrizeKnn(m))

  m <- as.dgcmatrix(matrix(rnorm(500, 0, 20), nrow = 20, ncol = 25))
  expect_squared(SymmetrizeKnn(m))
  expect_error(suppressWarnings(SymmetrizeKnn(t(m))))

  # m <- as.dgcmatrix(
  #   matrix(c(1, -5, 3, 3, 4,
  #            5, -3, 2, 6, 7,
  #            2, -9, 7, -.5, -4,
  #            -1, 6, 6, 2, -5,
  #            -1, 6, 3, -7.2, .2),
  #          nrow = 5, byrow = TRUE)
  # )

  m <- as.dgcmatrix(
    matrix(c(1, .5, 3,  3,   4,
             5,  2, 2,  6,   7,
             2,  9, 7, .5, 4.8,
             1,  6, 6,  2,   5,
             7, 6, 3, 7.2,  .2),
           nrow = 5, byrow = TRUE)
  )

  # m_sym_max <- m
  # m_sym_min <- m
  # for (i in 1:5) {
  #   for(j in i:5) {
  #     if (m_sym_max[i,j] < m_sym_max[j,i]) {
  #       m_sym_max[i,j] <- m_sym_max[j,i]
  #     } else {
  #       m_sym_max[j,i] <- m_sym_max[i,j]
  #     }
  #     if (m_sym_min[i,j] > m_sym_min[j,i]) {
  #       m_sym_min[i,j] <- m_sym_min[j,i]
  #     } else {
  #       m_sym_min[j,i] <- m_sym_min[i,j]
  #     }
  #   }
  # }
  expect_equal(SymmetrizeKnn(m, use.max = TRUE),
               naive_positive_matrix_symmetry(m))
  expect_equal(SymmetrizeKnn(m, use.max = FALSE),
               naive_positive_matrix_symmetry(m, FALSE))

  m[1,2] <- 0
  m[2,1] <- 0
  m[5,3] <- 0 # < m[3,5]
  m[3,2] <- 0 # > m[2,3]
  # m_sym_max <- m
  # m_sym_min <- m
  # for (i in 1:5) {
  #   for(j in i:5) {
  #     if (m_sym_max[i,j] == 0) {
  #       m_sym_max[i,j] <- m_sym_max[j,i]
  #       m_sym_min[i,j] <- m_sym_min[j,i]
  #     }
  #     if (m_sym_max[j,i] == 0) {
  #       m_sym_max[j,i] <- m_sym_max[i,j]
  #       m_sym_min[j,i] <- m_sym_min[i,j]
  #     }
  #     if (m_sym_max[i,j] < m_sym_max[j,i]) {
  #       m_sym_max[i,j] <- m_sym_max[j,i]
  #     } else {
  #       m_sym_max[j,i] <- m_sym_max[i,j]
  #     }
  #     if (m_sym_min[i,j] > m_sym_min[j,i]) {
  #       m_sym_min[i,j] <- m_sym_min[j,i]
  #     } else {
  #       m_sym_min[j,i] <- m_sym_min[i,j]
  #     }
  #   }
  # }
  expect_equal(SymmetrizeKnn(m, use.max = TRUE),
               naive_positive_matrix_symmetry(m))
  expect_equal(SymmetrizeKnn(m, use.max = FALSE),
               naive_positive_matrix_symmetry(m, FALSE))
})

test_that("matrix symmetry works on  a Seurat object", {
  m <- matrix(c(0,1,7,0,0,2,0,2,0,0,0,0,1,1,0,4,0,0,5,2,0,0,0,2,0,0,0,1,1,1,0,0,0,0,7,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,
                4,0,0,5,2,3,1,0,1,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,8,0,0,5,8,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,6,7,
                5,0,0,0,0,0,0,0,1,3,4,1,1,1,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,4,1,0,0,4,1,0,0,0,0,2,0,4,0,0,0,7,0,0,1,0,
                0,0,0,0,0,0,1,0,0,2,0,0,0,0,1,0,2,0,1,0,1,1,0,1,1,6,2,1,0,0,0,0,4,0,4,0,8,0,0,0,0,0,0,1,0,1,0,1,1,0,
                4,0,1,0,1,0,0,0,0,0,3,0,1,2,0,5,0,0,0,0,0,0,0,2,0,0,0,0,0,8,0,4,1,0,2,0,7,0,0,0,0,0,0,2,6,3,0,0,1,2,
                0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,2,0,0,0,0,0,2,0,2,0,0,0,4,0,1,0,0,0,0,5,0,0,0,8,5,0,2,0,0,0,0,0,0,2,1,
                1,0,0,0,1,4,0,0,2,4,0,0,1,0,2,0,2,4,0,0,4,0,0,0,0,2,0,2,0,0,0,2,2,0,0,0,0,0,0,0,5,0,7,0,0,3,1,0,0,3,
                4,1,6,0,1,0,0,2,0,0,0,0,2,6,1,0,2,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,2,5,0,0,0,0,0,0,0,2,0,1,1,0,0,
                2,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,3,0,0,0,2,4,0,0,0,0,0,0,7,1,0,0,0,8,0,0,3,0,0,0,0,0,0,2,0,0,2,0,0,
                4,0,4,1,0,1,0,8,2,8,2,0,0,0,0,1,0,5,0,0,0,0,0,0,0,0,3,2,8,2,2,4,2,0,0,0,0,1,0,0,3,1,0,0,0,1,0,0,2,0,
                0,0,3,0,0,0,1,0,0,0,0,0,7,0,0,0,0,0,0,0,0,2,0,0,2,0,1,0,0,1,0,0,0,0,8,1,0,0,0,2,0,2,8,0,0,2,2,1,0,1,
                1,6,0,0,0,2,1,0,0,2,0,2,0,0,0,0,2,8,2,2,0,0,0,3,0,1,0,0,0,0,2,0,0,3,0,0,2,0,2,0,6,1,0,0,0,0,0,0,0,0,
                1,0,4,4,3,1,8,5,0,2,0,0,4,0,6,0,0,0,0,2,0,0,0,1,0,0,6,0,1,0,0,0,1,0,3,0,0,8,0,2,0,0,0,0,8,0,0,0,4,2,
                2,7,0,2,0,0,0,0,2,0,0,8,0,1,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,1,0,1,1,0,2,0,0,0,0,0,0,1,0,0,1,3,0,0,0,3,
                0,1,0,0,0,0,0,1,1,2,0,0,1,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,5,0,0,2,0,0,0,0,0,0,0,0,0,0,3,0,
                0,4,1,0,2,0,0,2,0,0,1,0,2,1,0,0,0,0,0,0,0,0,0,0,0,7,0,3,0,1,3,1,2,0,3,0,0,0,1,0,0,0,0,0,2,1,0,0,1,2,
                8,0,2,0,0,2,0,0,0,0,0,3,0,0,0,0,0,4,0,2,0,0,0,0,0,1,1,0,0,2,2,8,0,0,0,0,1,0,0,1,4,0,0,1,0,0,0,0,4,0,
                2,0,0,0,0,0,2,7,0,0,0,0,0,0,0,0,0,0,0,2,2,0,2,0,0,0,0,0,0,0,1,0,0,0,0,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,
                0,0,0,0,6,0,7,1,0,0,0,2,0,0,3,1,0,0,2,1,0,2,0,2,2,0,0,0,0,2,0,0,0,0,0,0,0,0,4,0,0,0,0,0,2,0,0,0,0,0,
                2,0,0,5,2,0,0,0,1,0,0,0,0,0,2,1,2,0,0,0,2,0,0,0,2,5,0,0,2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,0,0,0,4,2,
                1,0,3,2,2,0,0,0,0,4,0,5,0,0,0,0,2,0,0,0,0,2,5,0,1,0,1,1,0,0,3,4,0,4,0,0,0,0,0,0,0,0,1,2,0,0,0,0,1,0,
                1,0,0,0,2,7,4,0,4,0,0,0,0,4,0,0,5,0,2,4,2,7,0,0,0,0,1,0,5,3,1,3,0,6,0,0,5,0,5,0,0,0,3,0,0,0,1,1,0,0,
                0,0,6,6,0,0,0,2,0,3,0,0,4,0,0,0,0,0,5,2,3,3,0,4,0,0,2,4,0,1,1,0,0,0,0,0,0,0,0,1,2,0,0,0,0,0,0,0,0,0,
                0,0,0,0,1,0,0,0,0,2,0,0,2,0,0,1,1,1,0,5,1,2,0,0,0,0,1,0,0,1,2,0,0,0,0,2,4,0,8,2,0,2,0,1,0,2,0,0,0,0,
                0,1,0,0,0,2,1,0,2,1,7,0,2,6,1,2,3,0,4,4,1,0,0,1,0,0,2,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,2,1,4,0,0,
                0,0,0,7,2,0,0,1,0,0,0,0,0,3,0,0,0,0,3,0,3,0,5,0,4,1,0,2,0,0,0,0,0,0,6,1,1,0,8,2,0,0,0,0,0,0,0,0,1,1,
                5,0,1,0,2,1,0,0,0,5,0,1,1,0,0,0,1,0,0,0,0,3,0,0,1,1,2,0,0,2,0,0,0,0,0,0,1,0,0,2,0,0,0,1,0,0,0,0,0,0,
                2,0,0,2,0,0,0,0,0,5,0,1,0,0,0,0,0,2,0,0,0,0,0,0,4,4,0,3,0,0,0,0,0,1,0,0,0,2,0,0,0,1,0,0,1,0,0,0,0,4,
                2,0,0,3,0,0,4,0,0,0,0,0,2,2,1,6,0,0,1,6,0,5,0,2,0,0,0,6,5,0,1,0,4,0,2,0,0,2,2,0,0,0,0,0,0,1,0,1,0,0,
                1,0,0,0,0,0,0,1,0,0,0,0,2,0,0,0,1,0,0,2,0,0,0,0,0,0,0,0,0,0,1,0,2,0,0,1,0,0,0,0,0,2,0,0,1,1,0,0,0,2,
                0,0,1,0,0,6,0,0,0,2,0,0,1,0,0,0,1,1,0,0,1,0,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,2,2,0,0,0,1,3,0,0,0,1,
                0,0,0,0,0,0,6,0,0,2,0,0,0,0,0,0,0,0,0,0,3,0,0,6,0,2,0,1,0,0,1,2,0,2,1,0,4,0,0,0,0,0,4,0,1,0,2,1,0,5,
                1,0,0,5,0,0,0,0,0,0,0,0,0,0,2,0,2,0,1,0,0,2,0,0,0,1,0,1,0,0,2,0,0,0,0,3,0,0,0,0,0,1,1,0,0,0,0,1,0,0,
                0,0,0,2,0,5,0,5,0,0,2,2,0,0,0,0,0,2,5,2,0,2,0,2,1,2,0,0,0,0,0,0,1,0,0,0,3,0,0,8,0,0,0,0,0,0,0,1,0,0,
                5,0,2,0,8,2,7,5,0,0,0,0,0,0,0,0,0,0,0,1,0,2,0,0,0,2,1,0,0,0,0,0,0,1,1,0,2,0,0,0,0,0,1,0,0,1,0,0,6,1,
                1,0,0,0,0,0,0,5,0,0,0,0,0,1,0,0,1,0,1,0,0,0,2,1,0,0,0,1,0,0,0,0,0,2,0,5,0,4,0,0,0,0,3,0,1,0,0,0,3,0,
                0,0,0,2,1,1,0,0,0,0,0,1,2,1,1,0,2,0,0,0,0,1,0,4,4,0,2,0,2,0,0,0,2,0,0,0,0,6,0,0,0,6,0,0,0,1,0,0,0,0,
                1,0,0,0,7,0,0,2,0,0,2,0,2,0,1,0,0,0,1,0,1,0,0,2,0,0,1,0,0,1,0,0,1,0,1,0,0,4,1,0,0,0,0,1,0,1,1,0,0,4,
                0,0,0,0,0,8,5,0,3,0,0,1,1,0,0,0,0,5,0,0,0,1,0,2,0,0,0,0,0,0,0,1,0,3,0,0,0,5,5,1,0,0,0,2,0,0,0,0,1,0,
                1,2,0,0,0,1,0,0,6,0,0,0,2,0,0,0,0,0,0,0,2,0,1,0,1,0,0,0,0,0,0,8,0,0,0,0,1,0,0,5,0,0,1,0,0,0,0,0,2,0,
                0,0,0,0,0,0,0,0,0,5,0,0,0,2,0,0,0,0,0,0,0,0,2,0,0,0,5,0,0,0,0,2,0,0,1,0,2,0,1,1,0,0,0,0,0,0,0,0,1,0,
                1,0,2,0,1,2,1,0,4,0,0,0,0,0,0,0,1,7,0,0,4,0,7,1,0,1,1,2,1,0,0,0,0,0,0,2,7,1,0,0,0,0,0,0,1,2,1,1,0,0,
                0,0,0,0,0,0,1,0,0,0,0,0,0,0,2,1,3,0,6,2,5,0,0,1,7,0,1,5,0,0,0,0,5,0,0,1,0,2,0,0,0,0,4,0,0,0,0,1,0,5,
                0,0,0,0,0,0,0,0,0,1,0,1,0,1,7,0,3,0,0,0,2,0,2,2,2,0,0,0,0,0,0,0,2,0,0,1,2,0,0,0,0,0,0,0,7,1,0,0,5,0,
                0,0,2,0,0,0,0,0,2,0,1,1,0,0,3,0,2,0,0,0,1,0,1,0,0,0,0,0,4,2,0,4,0,0,0,5,2,5,0,0,0,0,0,3,6,0,0,0,2,0,
                0,0,0,1,2,0,0,0,0,2,0,2,5,0,0,1,8,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,4,0,0,0,0,3,0,0,0,5,0,0,1,0,0,0,5,0,
                3,0,0,3,0,0,0,0,0,0,2,0,0,0,1,2,0,8,0,0,0,0,0,0,1,0,0,5,0,6,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,5,4,5,0,
                0,3,0,2,0,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,8,1,0,1,2,0,4,0,0,0,2,2,0,0,5,0,0,0,0,0,0,5,2,0,2,1,0,3,5,
                0,0,0,2,4,0,0,3,5,0,0,0,0,1,0,2,1,0,1,0,0,0,0,0,0,0,0,2,0,4,0,0,0,0,2,0,0,2,0,0,0,0,0,0,2,0,0,1,5,0,
                5,0,4,0,2,0,0,1,2,0,0,0,5,0,1,0,2,2,2,0,2,0,0,0,1,0,0,8,1,0,5,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,3,0),
              nrow = 50, dimnames = list(paste0("G", 1:50),
                                         paste0("C", 1:50)))
  seu <- PreProcessSO(CreateSeuratObject(as.dgcmatrix(m)), use_sct = FALSE,
                      nfeatures = 30, npcs = 15, approx = FALSE)
  seu <- suppressMessages(FindNeighbors(seu, dims = 1:15, k.param = 5L, return.neighbor = T, graph.name = "n_g1", verbose = FALSE))
  seu[["n_g2"]] <- seu[["n_g1"]]
  seu[["n_g2"]]@nn.dist[2,5] <- 45.678
  seu[["n_g2"]]@nn.idx[2,5] <- 10L
  seu[["n_g2"]]@nn.dist[49,5] <- 123.45
  seu[["g_g1"]] <- as.Graph(seu[["n_g1"]])
  seu[["g_g2"]] <- as.Graph(seu[["n_g2"]])
  seu <- suppressMessages(FindNeighbors(seu, dims = 1:15, k.param = 5L, return.neighbor = FALSE, graph.name = "g_g3", verbose = FALSE))

  r1_max <- SymmetrizeKnn(seu, 'n_g1', use.max = TRUE)[['n_g1_symmetric']]
  r1_min <- SymmetrizeKnn(seu, 'n_g1', use.max = FALSE)[['n_g1_symmetric']]
  expect_true(isSymmetric(r1_max))
  expect_true(isSymmetric(r1_min))
  expect_identical(r1_max, r1_min)
  expect_identical(dimnames(r1_max), dimnames(seu[['g_g1']]))
  expect_identical(dimnames(r1_max), dimnames(seu[['g_g3']]))
  expect_false(isTRUE(all.equal(seu[['g_g1']], r1_max)))
  expect_false(is.kconstant(r1_max))
  expect_true(get.k(r1_max, "min") == 5L)

  r2_max <- SymmetrizeKnn(seu, 'n_g2', use.max = TRUE)[['n_g2_symmetric']]
  r2_min <- SymmetrizeKnn(seu, 'n_g2', use.max = FALSE)[['n_g2_symmetric']]
  expect_true(isSymmetric(r2_max))
  expect_true(isSymmetric(r2_min))
  expect_false(all(r2_max == r2_min))
  expect_identical(dimnames(r2_max), dimnames(r2_min))
  expect_identical(dimnames(r2_max), dimnames(seu[['g_g1']]))
  expect_identical(dimnames(r2_max), dimnames(seu[['g_g3']]))
  expect_false(isTRUE(all.equal(seu[['g_g2']], r2_max)))
  expect_false(is.kconstant(r2_max))
  expect_true(get.k(r2_max, "min") == 5L)

  r3_max <- SymmetrizeKnn(seu, 'g_g3', use.max = TRUE)[['g_g3_symmetric']]
  r3_min <- SymmetrizeKnn(seu, 'g_g3', use.max = FALSE)[['g_g3_symmetric']]
  expect_true(isSymmetric(r3_max))
  expect_true(isSymmetric(r3_min))
  expect_identical(r3_max, r3_min)
  expect_identical(dimnames(r3_max), dimnames(r3_min))
  expect_identical(dimnames(r3_max), dimnames(seu[['g_g1']]))
  expect_identical(dimnames(r3_max), dimnames(seu[['g_g3']]))
  expect_false(isTRUE(all.equal(seu[['g_g3']], r3_max)))
  expect_false(is.kconstant(r3_max))
  expect_true(get.k(r3_max, "min") == 5L)

  expected <- naive_positive_matrix_symmetry(seu[['g_g1']], use_max = TRUE,  assay.used = "RNA")
  expect_identical(r1_max, expected)
  expect_identical(SymmetrizeKnn(seu, graph.name = "g_g1", use.max = TRUE)[['g_g1_symmetric']], expected)
  expected <- naive_positive_matrix_symmetry(seu[['g_g1']], use_max = FALSE, assay.used = "RNA")
  expect_identical(r1_max, expected)
  expect_identical(SymmetrizeKnn(seu, graph.name = "g_g1", use.max = FALSE)[['g_g1_symmetric']], expected)

  expected <- naive_positive_matrix_symmetry(seu[['g_g2']], use_max = TRUE,  assay.used = "RNA")
  expect_identical(r2_max, expected)
  expect_identical(SymmetrizeKnn(seu, graph.name = "g_g2", use.max = TRUE)[['g_g2_symmetric']], expected)
  expected <- naive_positive_matrix_symmetry(seu[['g_g2']], use_max = FALSE, assay.used = "RNA")
  expect_identical(r2_min, expected)
  expect_identical(SymmetrizeKnn(seu, graph.name = "g_g2", use.max = FALSE)[['g_g2_symmetric']], expected)

  expected <- naive_positive_matrix_symmetry(seu[['g_g3']], use_max = TRUE,  assay.used = "RNA")
  expect_identical(r3_max, expected)
  expect_identical(SymmetrizeKnn(seu, graph.name = "g_g3", use.max = TRUE)[['g_g3_symmetric']], expected)
  expected <- naive_positive_matrix_symmetry(seu[['g_g3']], use_max = FALSE, assay.used = "RNA")
  expect_identical(r3_max, expected)
  expect_identical(SymmetrizeKnn(seu, graph.name = "g_g3", use.max = FALSE)[['g_g3_symmetric']], expected)
})

# test_that("L2 normalisation produce same results as Seurat", {
#   nrow <- sample(50:200, 1)
#   ncol <- sample(100:400, 1)
#
#   m <- matrix(rnorm(nrow*ncol, 3, sd = 2), nrow = nrow, ncol = ncol)
#   m[m < 0] <- 0
#
#   m <- as.matrix(dist(m))
#   expect_equal(NormaliseL2(m, MARGIN = 1), Seurat:::L2Norm(m, MARGIN = 1))
#   expect_equal(NormaliseL2(m, MARGIN = 2), Seurat:::L2Norm(m, MARGIN = 2))
#
#   m <- as.dgcmatrix(m)
#   expect_equal(NormaliseL2(m, MARGIN = 1), Seurat:::L2Norm(m, MARGIN = 1))
#   expect_equal(NormaliseL2(m, MARGIN = 2), Seurat:::L2Norm(m, MARGIN = 2))
# })

test_that("knn trimming works on a Neighbor Object", {
  knn_dist <- matrix(rep(0, 12), ncol = 1)
  knn_idx <- matrix(1:12, ncol = 1)
  knn_object <- SeuratObject:::Neighbor(nn.idx = knn_idx, nn.dist = knn_dist,
                                        cell.names = paste0("Cell_", 1:12))

  expect_equal(.cut.knn(knn_object, k.max = 1), naive_knn_sort(knn_object))
  expect_equal(.CutKnn(knn_object, k.max = 1, verbose = FALSE), naive_knn_sort(knn_object))
  expect_error(.cut.knn(knn_object, k.max = 2), "subscript out of bounds$")
  expect_equal(.CutKnn(knn_object, k.max = 2, verbose = FALSE), naive_knn_sort(knn_object))

  knn_idx <- cbind(knn_idx, matrix(c(2, 5, 6, 7, 8, 10, 12,
                                     1, 3, 4, 6, 8,  9, 10,
                                     2, 4, 5, 6, 7,  9, 12,
                                     2, 3, 6, 7, 9, 10, 11,
                                     1, 3, 6, 7, 8,  9, 11,
                                     1, 2, 3, 4, 5,  9, 11,
                                     1, 3, 4, 5, 6,  8, 11,
                                     1, 2, 5, 7, 9, 10, 12,
                                     2, 3, 4, 5, 6,  8, 10,
                                     1, 2, 4, 8, 9, 11, 12,
                                     2, 4, 5, 6, 7, 10, 12,
                                     1, 3, 5, 7, 8, 10, 11),
                                   nrow = 12, byrow = TRUE))
  knn_dist <- cbind(knn_dist, matrix(c(8.9, .1, .5, .7, .4, 2, 12,
                                       8.9, .0098, 4, 22, 6.3, 4, .5,
                                       .0098, 4.992, 10, 10, 2.987, 82, 47,
                                       4, 4.992, 7, 0.518, 2.091, 17.8, 9.2,
                                       .1, 10, 7, 5.0927, 40.17, 33, 7,
                                       .5, 22, 10, 7, 7, 3, 7,
                                       .7, 2.987, 0.518, 5.0927, 28.3, 2.987, 2.987,
                                       .4, 6.3, 40.17, 2.987, 2.4, 45.9, 2,
                                       4, 82, 2.091, 33, 3, 2.4, 26.1,
                                       2, .5, 17.8, 45.9, 26.1, 67, .8,
                                       72, 9.2, 7, 7, 2.987, 67, 58.789,
                                       12, 47, 86, 2, 2, .8, 58.789),
                                     nrow = 12, byrow = TRUE))

  knn_object <- SeuratObject:::Neighbor(nn.idx = knn_idx, nn.dist = knn_dist,
                                        cell.names = paste0("Cell_", 1:12))

  for (k_value in 1:8) {
    # expected <- knn_object_ordered
    # expected@nn.idx <- expected@nn.idx[,1:k_value, drop = F]
    # expected@nn.dist <- expected@nn.dist[,1:k_value, drop = F]
    expect_equal(.cut.knn(knn_object, k.max = k_value),
                 naive_knn_sort(knn_object, k_value))
    # Problem
    expect_equal(.CutKnn(knn_object, k.max = k_value, , verbose = FALSE),
                 if (k_value < 8) naive_knn_sort(knn_object, k_value) else knn_object)
    # .......
  }
})

test_that("knn trimming works on a Graph Object", {
  knn_dist <- matrix(rep(0, 12), ncol = 1)
  knn_idx <- matrix(1:12, ncol = 1)
  knn_object <- as.Graph(
    SeuratObject:::Neighbor(nn.idx = knn_idx, nn.dist = knn_dist,
                            cell.names = paste0("Cell_", 1:12))
  )
  knn_object@assay.used <- "RNA"

  expect_equal(.cut.knn(knn_object, k.max = 1), knn_object)
  expect_equal(.cut.knn(knn_object, k.max = 2), knn_object)

  # Problem
  expect_equal(as.dgcmatrix(.CutKnn(knn_object, k.max = 1, verbose = FALSE)), drop0(knn_object))
  expect_equal(as.dgcmatrix(.CutKnn(knn_object, k.max = 2, verbose = FALSE)), drop0(knn_object))
  # .......
  # expect_error(.cut.knn(knn_object, k.max = 2), "subscript out of bounds$")

  knn_idx <- cbind(knn_idx, matrix(c(2, 5, 6, 7, 8, 10, 12,
                                     1, 3, 4, 6, 8,  9, 10,
                                     2, 4, 5, 6, 7,  9, 12,
                                     2, 3, 6, 7, 9, 10, 11,
                                     1, 3, 6, 7, 8,  9, 11,
                                     1, 2, 3, 4, 5,  9, 11,
                                     1, 3, 4, 5, 6,  8, 11,
                                     1, 2, 5, 7, 9, 10, 12,
                                     2, 3, 4, 5, 6,  8, 10,
                                     1, 2, 4, 8, 9, 11, 12,
                                     2, 4, 5, 6, 7, 10, 12,
                                     1, 3, 5, 7, 8, 10, 11),
                                   nrow = 12, byrow = TRUE))
  knn_dist <- cbind(knn_dist, matrix(c(8.9, .1, .5, .7, .4, 2, 12,
                                       8.9, .0098, 4, 22, 6.3, 4, .5,
                                       .0098, 4.992, 10, 10, 2.987, 82, 47,
                                       4, 4.992, 7, 0.518, 2.091, 17.8, 9.2,
                                       .1, 10, 7, 5.0927, 40.17, 33, 7,
                                       .5, 22, 10, 7, 7, 3, 7,
                                       .7, 2.987, 0.518, 5.0927, 28.3, 2.987, 2.987,
                                       .4, 6.3, 40.17, 2.987, 2.4, 45.9, 2,
                                       4, 82, 2.091, 33, 3, 2.4, 26.1,
                                       2, .5, 17.8, 45.9, 26.1, 67, .8,
                                       72, 9.2, 7, 7, 2.987, 67, 58.789,
                                       12, 47, 86, 2, 2, .8, 58.789),
                                     nrow = 12, byrow = TRUE))

  knn_object <- as.Graph(
    SeuratObject:::Neighbor(nn.idx = knn_idx, nn.dist = knn_dist,
                            cell.names = paste0("Cell_", 1:12))
  )
  knn_object@assay.used <- "RNA"


  empty_knn <- sparseMatrix(i=integer(0), p = rep(0,12), x = numeric(0),
                            dims = c(12,12),
                            dimnames = list(paste0("Cell_", 1:12),
                                            paste0("Cell_", 1:12)))
  empty_knn <- as.Graph(empty_knn)
  empty_knn@assay.used <- 'RNA'
  expect_equal(.cut.knn(knn_object, k.max = 0), empty_knn)
  # Problem
  expect_error(.CutKnn(knn_object, k.max = 0, verbose = F))#, empty_knn)
  # .......
  for (k_value in 1:10) {
    expected <- naive_prune_graph(knn_object, k_value = k_value)
    expect_equal(.cut.knn(knn_object, k.max = k_value), expected)
  }
  knn_object <- as.Graph(SymmetrizeKnn(knn_object))
  knn_object@assay.used <- "RNA"
  expect_equal(.cut.knn(knn_object, k.max = 0), empty_knn)
  for (k_value in 1:10) {
    expected <- naive_prune_graph(knn_object, k_value = k_value)
    expect_equal(.cut.knn(knn_object, k.max = k_value), expected)
  }
})

test_that("knn trimming works on a Seurat Object", {
  seu <- PreProcessSO(liver_small, use_sct = TRUE, npcs = 10)
  seu <- suppressMessages(FindNeighbors(seu, dims = 1:10, return.neighbor = TRUE, verbose = FALSE))

  missing_err <- 'missing, with no default'
  null_err <- 'graph.name cannot be null'
  # method_err1 <- 'no applicable method for \'\\w+\' applied to an object of class'
  # method_err2 <- 'unable to find an inherited method for function .[[:alnum:]_\\.]+. for signature .object = .[[:alnum:]_\\.]+..'
  method_err <- paste(
    'no applicable method for .[[:alnum:]_\\.]+. applied to an object of class',
    'unable to find an inherited method for function .[[:alnum:]_\\.]+. for signature .object = .[[:alnum:]_\\.]+..',
    sep = "|"
  )
  graph_err <- 'not found in this Seurat object'

  expect_error(CutKnn(), missing_err)
  expect_error(CutKnn(graph.name = NULL), null_err)
  expect_error(CutKnn(graph.name = "whatever", k.max = 20), missing_err)
  expect_error(CutKnn(object = 3, graph.name = "whatever", k.max = 20), method_err)
  expect_error(CutKnn(object = seu, graph.name = 7, k.max = 20), method_err)
  expect_error(CutKnn(object = seu, graph.name = "whatever", k.max = 20), graph_err)
  expect_error(CutKnn(object = seu, graph.name = "whatever", k.max = 20), graph_err)

  expect_identical(
    CutKnn(object = seu, graph.name = "SCT.nn", k.max = 20, new.graph = 'pruned'),
    CutKnn(object = seu, graph.name = "SCT.nn", k.max = 50, 'pruned'))
  expect_identical(
    CutKnn(object = seu, graph.name = "SCT.nn", k.max = 20, new.graph = 'pruned')@neighbors$pruned,
    seu@neighbors$SCT.nn)

  expect_equal(CutKnn(object = seu, graph.name = "SCT.nn", k.max = 1,
                      new.graph = 'pruned', verbose = F)[['pruned']],
               naive_knn_sort(seu[['SCT.nn']], 2))

  for (k_value in 2:20) {
    expect_equal(CutKnn(object = seu, graph.name = "SCT.nn", k.max = k_value, new.graph = 'pruned', verbose = F)@neighbors$pruned,
                 naive_knn_sort(seu[['SCT.nn']], cut_k = k_value))
  }

  seu <- SymmetrizeKnn(seu, graph.name = 'SCT.nn')
  expect_equal(
    CutKnn(object = seu, graph.name = "SCT.nn_symmetric", k.max = 1, new.graph = 'pruned', verbose = F)@graphs$pruned,
    CutKnn(object = seu, graph.name = "SCT.nn_symmetric", k.max = 2, new.graph = 'pruned', verbose = F)@graphs$pruned
  )
  expect_equal(
    CutKnn(object = seu, graph.name = "SCT.nn_symmetric", k.max = 1, new.graph = 'pruned', verbose = F)@graphs$pruned,
    naive_prune_graph(seu[['SCT.nn_symmetric']], k_value = 1)
  )
  for (k_value in 3:20) {
    # Problem
    expect_equal(CutKnn(object = seu, graph.name = "SCT.nn_symmetric", k.max = k_value, new.graph = 'pruned', verbose = F)@graphs$pruned,
                 naive_prune_graph(seu[['SCT.nn_symmetric']], k_value - 1))
    # .......
  }
})

test_that("retrieval of best neighbours in a knn Graph is correct", {
  m <- as.dgcmatrix(matrix(rnorm(500, 50, 20), nrow = 20, ncol = 25))
  m[m < 0] <- rnorm(sum(m < 0), 20, 2)
  expect_error(get.k.best.neighbours(m, k.max = 30), "can't get best 30 neighbors for all cells")
  expect_equal(get.k.best.neighbours(m, k.max = 0), matrix(NA_real_, 0, 0))
  for (i in 1:25) {
    expect_equal(get.k.best.neighbours(m, k.max = i, graph.type = "distances"),
                 naive_get_closest(m, k.max = i, decreasing = FALSE))
  }

  for (i in 1:nrow(m)) {
    m[i, sample.int(25, 15)] <- 0
  }
  for (i in 1:8) {
    expect_equal(get.k.best.neighbours(m, k.max = i, graph.type = "distances"),
                 naive_get_closest(m, k.max = i, decreasing = FALSE))
  }
  for (i in 1:nrow(m)) {
    m[i, m[i,] != 0] <- rescale(- m[i, m[i,] != 0], to = 0:1)
  }
  for (i in 1:8) {
    expect_equal(get.k.best.neighbours(m, k.max = i, graph.type = "connectivities"),
                 naive_get_closest(m, k.max = i, decreasing = TRUE))
  }
})

test_that("computation of number of neighbours between batches works on a Neighbor object", {
  knn_dist <- matrix(rep(0, 12), ncol = 1)
  knn_idx <- matrix(1:12, ncol = 1)
  knn_object <- SeuratObject:::Neighbor(nn.idx = knn_idx, nn.dist = knn_dist,
                                        cell.names = paste0("Cell_", 1:12))
  expect_error(GetNeighborsPerBatch(knn_object), "missing, with no default")

  batch <- data.frame(batch = rep("A", 12), row.names = paste0("Cell_", 1:12))

  expect_error(GetNeighborsPerBatch(knn_object, batch, count.self = TRUE),
               "contrasts can be applied only to factors with 2 or more levels")
  batch[sample.int(12, 6),1] <- "B"
  expect_error(GetNeighborsPerBatch(knn_object, batch$batch))
  expect_error(GetNeighborsPerBatch(knn_object, tibble::remove_rownames(batch)),
               "subscript out of bounds")
  expected <- matrix(0, nrow = 2, ncol = 2, dimnames =
                       list(sort(paste0("batch", unique(batch$batch))),
                            sort(paste0("batch", unique(batch$batch)))))
  expect_equal(GetNeighborsPerBatch(knn_object, batch, count.self = FALSE),
               expected)
  diag(expected) <- 6
  expect_equal(GetNeighborsPerBatch(knn_object, batch, count.self = TRUE),
               expected)

  batch$batch <- sample(c(rep("A", 4), rep("B", 2), rep("C", 1), rep("D", 5)))
  expected <- matrix(0, nrow = 4, ncol = 4, dimnames =
                       list(sort(paste0("batch", unique(batch$batch))),
                            sort(paste0("batch", unique(batch$batch)))))
  expect_equal(GetNeighborsPerBatch(knn_object, batch, count.self = FALSE),
               expected)
  diag(expected) <- c(4, 2, 1, 5)
  expect_equal(GetNeighborsPerBatch(knn_object, batch, count.self = TRUE),
               expected)

  knn_dist <- matrix(rep(0, 12), ncol = 1)
  knn_idx <- matrix(1:12, ncol = 1)
  knn_idx <- cbind(knn_idx, matrix(c(2, 5, 6, 7, 8, 10, 12,
                                     1, 3, 4, 6, 8,  9, 10,
                                     2, 4, 5, 6, 7,  9, 12,
                                     2, 3, 6, 7, 9, 10, 11,
                                     1, 3, 6, 7, 8,  9, 11,
                                     1, 2, 3, 4, 5,  9, 11,
                                     1, 3, 4, 5, 6,  8, 11,
                                     1, 2, 5, 7, 9, 10, 12,
                                     2, 3, 4, 5, 6,  8, 10,
                                     1, 2, 4, 8, 9, 11, 12,
                                     2, 4, 5, 6, 7, 10, 12,
                                     1, 3, 5, 7, 8, 10, 11),
                                   nrow = 12, byrow = TRUE))
  knn_dist <- cbind(knn_dist, matrix(c(8.9, .1, .5, .7, .4, 2, 12,
                                       8.9, .0098, 4, 22, 6.3, 4, .5,
                                       .0098, 4.992, 10, 10, 2.987, 82, 47,
                                       4, 4.992, 7, 0.518, 2.091, 17.8, 9.2,
                                       .1, 10, 7, 5.0927, 40.17, 33, 7,
                                       .5, 22, 10, 7, 7, 3, 7,
                                       .7, 2.987, 0.518, 5.0927, 28.3, 2.987, 2.987,
                                       .4, 6.3, 40.17, 2.987, 2.4, 45.9, 2,
                                       4, 82, 2.091, 33, 3, 2.4, 26.1,
                                       2, .5, 17.8, 45.9, 26.1, 67, .8,
                                       72, 9.2, 7, 7, 2.987, 67, 58.789,
                                       12, 47, 86, 2, 2, .8, 58.789),
                                     nrow = 12, byrow = TRUE))

  knn_object <- SeuratObject:::Neighbor(nn.idx = knn_idx, nn.dist = knn_dist,
                                        cell.names = paste0("Cell_", 1:12))

  expected <- naive_knn_per_batch(knn_object, batch$batch)
  expect_equal(GetNeighborsPerBatch(knn_object, batch, count.self = TRUE),
               expected)

  diag(expected) <- diag(expected) - c(4, 2, 1, 5)
  expect_equal(GetNeighborsPerBatch(knn_object, batch, count.self = FALSE),
               expected)

  batch_f <- batch %>%
    mutate(batch = factor(batch, levels = LETTERS[1:10]))

  expected <- naive_knn_per_batch(knn_object, batch_f$batch)
  expect_equal(GetNeighborsPerBatch(knn_object, batch_f, count.self = TRUE),
               expected)
  batch_f <- batch %>%
    mutate(batch = factor(batch, levels = sample(LETTERS[1:10])))
  expect_equal(GetNeighborsPerBatch(knn_object, batch_f, count.self = TRUE),
               expected[paste0("batch", levels(batch_f$batch)),
                        paste0("batch", levels(batch_f$batch))])

  diag(expected)[1:4] <- diag(expected)[1:4] - c(4, 2, 1, 5)
  batch_f <- batch %>%
    mutate(batch = factor(batch, levels = LETTERS[1:10]))
  expect_equal(GetNeighborsPerBatch(knn_object, batch_f, count.self = FALSE),
               expected)
  batch_f <- batch %>%
    mutate(batch = factor(batch, levels = sample(LETTERS[1:10])))
  expect_equal(GetNeighborsPerBatch(knn_object, batch_f, count.self = FALSE),
               expected[paste0("batch", levels(batch_f$batch)),
                        paste0("batch", levels(batch_f$batch))])

})

test_that("computation of number of neighbours between batches works on a Graph object", {
  knn_dist <- matrix(rep(0, 12), ncol = 1)
  knn_idx <- matrix(1:12, ncol = 1)
  knn_object <- SeuratObject:::Neighbor(nn.idx = knn_idx, nn.dist = knn_dist,
                                        cell.names = paste0("Cell_", 1:12))

  expect_error(GetNeighborsPerBatch(as.Graph(knn_object)), "missing, with no default")

  batch <- data.frame(batch = rep("A", 12), row.names = paste0("Cell_", 1:12))

  expect_error(GetNeighborsPerBatch(as.Graph(knn_object), batch, count.self = TRUE),
               "contrasts can be applied only to factors with 2 or more levels")
  batch[sample.int(12, 6),1] <- "B"
  expect_error(GetNeighborsPerBatch(as.Graph(knn_object), batch$batch))
  expect_error(GetNeighborsPerBatch(as.Graph(knn_object), tibble::remove_rownames(batch)),
               "subscript out of bounds")

  expect_equal(GetNeighborsPerBatch(as.Graph(knn_object), batch, count.self = FALSE),
               GetNeighborsPerBatch(knn_object, batch, count.self = FALSE))

  expect_equal(GetNeighborsPerBatch(as.Graph(knn_object), batch, count.self = TRUE),
               GetNeighborsPerBatch(knn_object, batch, count.self = TRUE))

  batch$batch <- sample(c(rep("A", 4), rep("B", 2), rep("C", 1), rep("D", 5)))

  expect_equal(GetNeighborsPerBatch(as.Graph(knn_object), batch, count.self = FALSE),
               GetNeighborsPerBatch(knn_object, batch, count.self = FALSE))

  expect_equal(GetNeighborsPerBatch(as.Graph(knn_object), batch, count.self = TRUE),
               GetNeighborsPerBatch(knn_object, batch, count.self = TRUE))

  knn_dist <- matrix(rep(0, 12), ncol = 1)
  knn_idx <- matrix(1:12, ncol = 1)
  knn_idx <- cbind(knn_idx, matrix(c(2, 5, 6, 7, 8, 10, 12,
                                     1, 3, 4, 6, 8,  9, 10,
                                     2, 4, 5, 6, 7,  9, 12,
                                     2, 3, 6, 7, 9, 10, 11,
                                     1, 3, 6, 7, 8,  9, 11,
                                     1, 2, 3, 4, 5,  9, 11,
                                     1, 3, 4, 5, 6,  8, 11,
                                     1, 2, 5, 7, 9, 10, 12,
                                     2, 3, 4, 5, 6,  8, 10,
                                     1, 2, 4, 8, 9, 11, 12,
                                     2, 4, 5, 6, 7, 10, 12,
                                     1, 3, 5, 7, 8, 10, 11),
                                   nrow = 12, byrow = TRUE))
  knn_dist <- cbind(knn_dist, matrix(c(8.9, .1, .5, .7, .4, 2, 12,
                                       8.9, .0098, 4, 22, 6.3, 4, .5,
                                       .0098, 4.992, 10, 10, 2.987, 82, 47,
                                       4, 4.992, 7, 0.518, 2.091, 17.8, 9.2,
                                       .1, 10, 7, 5.0927, 40.17, 33, 7,
                                       .5, 22, 10, 7, 7, 3, 7,
                                       .7, 2.987, 0.518, 5.0927, 28.3, 2.987, 2.987,
                                       .4, 6.3, 40.17, 2.987, 2.4, 45.9, 2,
                                       4, 82, 2.091, 33, 3, 2.4, 26.1,
                                       2, .5, 17.8, 45.9, 26.1, 67, .8,
                                       72, 9.2, 7, 7, 2.987, 67, 58.789,
                                       12, 47, 86, 2, 2, .8, 58.789),
                                     nrow = 12, byrow = TRUE))

  knn_object <- SeuratObject:::Neighbor(nn.idx = knn_idx, nn.dist = knn_dist,
                                        cell.names = paste0("Cell_", 1:12))

  expect_equal(GetNeighborsPerBatch(as.Graph(knn_object), batch, count.self = TRUE),
               GetNeighborsPerBatch(knn_object, batch, count.self = TRUE))
  expect_equal(GetNeighborsPerBatch(as.Graph(knn_object), batch, count.self = FALSE),
               GetNeighborsPerBatch(knn_object, batch, count.self = FALSE))

  batch_f <- batch %>%
    mutate(batch = factor(batch, levels = LETTERS[1:10]))

  expect_equal(GetNeighborsPerBatch(as.Graph(knn_object), batch_f, count.self = TRUE),
               GetNeighborsPerBatch(knn_object, batch_f, count.self = TRUE))
  expect_equal(GetNeighborsPerBatch(as.Graph(knn_object), batch_f, count.self = FALSE),
               GetNeighborsPerBatch(knn_object, batch_f, count.self = FALSE))
  batch_f <- batch %>%
    mutate(batch = factor(batch, levels = sample(LETTERS[1:10])))
  expect_equal(GetNeighborsPerBatch(as.Graph(knn_object), batch_f, count.self = TRUE),
               GetNeighborsPerBatch(knn_object, batch_f, count.self = TRUE))

})

test_that("computation of number of neighbours between batches works on a Seurat Object", {
  seu <- PreProcessSO(liver_small, use_sct = TRUE, npcs = 10)
  seu <- suppressMessages(FindNeighbors(seu, dims = 1:10, return.neighbor = TRUE, verbose = FALSE))
  seu <- SymmetrizeKnn(seu, "SCT.nn")
  seu[['SCT.nn_dist']] <- as.Graph(seu[['SCT.nn']])
  seu <- suppressMessages(FindNeighbors(seu, dims = 1:10, return.neighbor = FALSE, verbose = FALSE))

  expect_error(GetNeighborsPerBatch(seu), "batch.var = \"missing\"")
  expect_error(GetNeighborsPerBatch(seu, "First_author"), "not found in this Seurat object")
  expect_error(GetNeighborsPerBatch(seu, "boulgiboulga"), "undefined columns selected")
  expect_error(GetNeighborsPerBatch(seu, seu[[]][, "First_author"]), "undefined columns selected")

  expected <- naive_knn_per_batch(seu[['SCT.nn']], seu$First_author, 'First_author')

  expect_equal(GetNeighborsPerBatch(seu, 'First_author', graph.name = 'SCT.nn'),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, seu[[]][,'First_author', drop = F], graph.name = 'SCT.nn'),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, seu[[]][, c('First_author', 'ID_sample')], graph.name = 'SCT.nn'),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, 'First_author', graph.name = 'SCT.nn_dist'),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, 'First_author', graph.name = 'SCT_nn'),
               expected)
  expect_equal(GetPropInterBatch(seu, 'First_author', graph.name = 'SCT_nn'),
               naive_prop_inter(expected))
  expect_equal(GetPropInterBatch(seu, 'First_author', graph.name = 'SCT_nn', per.batch = FALSE),
               naive_prop_inter(expected, per.batch = FALSE))
  expect_equal(GetPropIntraBatch(seu, 'First_author', graph.name = 'SCT_nn'),
               naive_prop_intra(expected))
  expect_equal(GetPropIntraBatch(seu, 'First_author', graph.name = 'SCT_nn', per.batch = FALSE),
               naive_prop_intra(expected, per.batch = FALSE))

  diag(expected) <- diag(expected) - table(seu$First_author)[c('Saleh', 'Song')]

  expect_equal(GetNeighborsPerBatch(seu, 'First_author', graph.name = 'SCT.nn', count.self = FALSE),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, seu[[]][,'First_author', drop = F], graph.name = 'SCT.nn', count.self = FALSE),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, seu[[]][, c('First_author', 'ID_sample')], graph.name = 'SCT.nn', count.self = FALSE),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, 'First_author', graph.name = 'SCT.nn_dist', count.self = FALSE),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, 'First_author', graph.name = 'SCT_nn', count.self = FALSE),
               expected)
  expect_equal(GetPropInterBatch(seu, 'First_author', graph.name = 'SCT_nn', count.self = FALSE),
               naive_prop_inter(expected))
  expect_equal(GetPropInterBatch(seu, 'First_author', graph.name = 'SCT_nn', count.self = FALSE, per.batch = FALSE),
               naive_prop_inter(expected, per.batch = FALSE))
  expect_equal(GetPropIntraBatch(seu, 'First_author', graph.name = 'SCT_nn', count.self = FALSE),
               naive_prop_intra(expected))
  expect_equal(GetPropIntraBatch(seu, 'First_author', graph.name = 'SCT_nn', count.self = FALSE, per.batch = FALSE),
               naive_prop_intra(expected, per.batch = FALSE))


  expected[] <- c(1651, 262, 262, 2931)
  expect_equal(GetNeighborsPerBatch(seu, 'First_author', graph.name = 'SCT.nn_symmetric'),
               expected)
  expect_equal(GetPropInterBatch(seu, 'First_author', graph.name = 'SCT.nn_symmetric'),
               naive_prop_inter(expected))
  expect_equal(GetPropInterBatch(seu, 'First_author', graph.name = 'SCT.nn_symmetric', per.batch = FALSE),
               naive_prop_inter(expected, per.batch = FALSE))
  expect_equal(GetPropIntraBatch(seu, 'First_author', graph.name = 'SCT.nn_symmetric'),
               naive_prop_intra(expected))
  expect_equal(GetPropIntraBatch(seu, 'First_author', graph.name = 'SCT.nn_symmetric', per.batch = FALSE),
               naive_prop_intra(expected, per.batch = FALSE))
  diag(expected) <- diag(expected) - table(seu$First_author)[c('Saleh', 'Song')]
  expect_equal(GetNeighborsPerBatch(seu, 'First_author', graph.name = 'SCT.nn_symmetric', count.self = FALSE),
               expected)
  expect_equal(GetPropInterBatch(seu, 'First_author', graph.name = 'SCT.nn_symmetric', count.self = FALSE),
               naive_prop_inter(expected))
  expect_equal(GetPropInterBatch(seu, 'First_author', graph.name = 'SCT.nn_symmetric', count.self = FALSE, per.batch = FALSE),
               naive_prop_inter(expected, per.batch = FALSE))
  expect_equal(GetPropIntraBatch(seu, 'First_author', graph.name = 'SCT.nn_symmetric', count.self = FALSE),
               naive_prop_intra(expected))
  expect_equal(GetPropIntraBatch(seu, 'First_author', graph.name = 'SCT.nn_symmetric', count.self = FALSE, per.batch = FALSE),
               naive_prop_intra(expected, per.batch = FALSE))
  expected[] <- c(2513, 710, 710, 5527)
  expect_equal(GetNeighborsPerBatch(seu, 'First_author', graph.name = 'SCT_snn'),
               expected)
  expect_equal(GetPropInterBatch(seu, 'First_author', graph.name = 'SCT_snn'),
               naive_prop_inter(expected))
  expect_equal(GetPropInterBatch(seu, 'First_author', graph.name = 'SCT_snn', per.batch = FALSE),
               naive_prop_inter(expected, per.batch = FALSE))
  expect_equal(GetPropIntraBatch(seu, 'First_author', graph.name = 'SCT_snn'),
               naive_prop_intra(expected))
  expect_equal(GetPropIntraBatch(seu, 'First_author', graph.name = 'SCT_snn', per.batch = FALSE),
               naive_prop_intra(expected, per.batch = FALSE))
  diag(expected) <- diag(expected) - table(seu$First_author)[c('Saleh', 'Song')]
  expect_equal(GetNeighborsPerBatch(seu, 'First_author', graph.name = 'SCT_snn', count.self = FALSE),
               expected)
  expect_equal(GetPropInterBatch(seu, 'First_author', graph.name = 'SCT_snn', count.self = FALSE),
               naive_prop_inter(expected))
  expect_equal(GetPropInterBatch(seu, 'First_author', graph.name = 'SCT_snn', count.self = FALSE, per.batch = FALSE),
               naive_prop_inter(expected, per.batch = FALSE))
  expect_equal(GetPropIntraBatch(seu, 'First_author', graph.name = 'SCT_snn', count.self = FALSE),
               naive_prop_intra(expected))
  expect_equal(GetPropIntraBatch(seu, 'First_author', graph.name = 'SCT_snn', count.self = FALSE, per.batch = FALSE),
               naive_prop_intra(expected, per.batch = FALSE))

  expected <- naive_knn_per_batch(seu[['SCT.nn']], seu$ID_sample, 'ID_sample')

  expect_equal(GetNeighborsPerBatch(seu, 'ID_sample', graph.name = 'SCT.nn'),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, seu[[]][,'ID_sample', drop = F], graph.name = 'SCT.nn'),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, seu[[]][, c('ID_sample', 'First_author')], graph.name = 'SCT.nn'),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, 'ID_sample', graph.name = 'SCT.nn_dist'),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, 'ID_sample', graph.name = 'SCT_nn'),
               expected)
  expect_equal(GetPropInterBatch(seu, 'ID_sample', graph.name = 'SCT_nn'),
               naive_prop_inter(expected))
  expect_equal(GetPropInterBatch(seu, 'ID_sample', graph.name = 'SCT_nn', per.batch = FALSE),
               naive_prop_inter(expected, per.batch = FALSE))
  expect_equal(GetPropIntraBatch(seu, 'ID_sample', graph.name = 'SCT_nn'),
               naive_prop_intra(expected))
  expect_equal(GetPropIntraBatch(seu, 'ID_sample', graph.name = 'SCT_nn', per.batch = FALSE),
               naive_prop_intra(expected, per.batch = FALSE))

  diag(expected) <- diag(expected) - table(seu$ID_sample)[c('B_13_TL', 'B_18_TL', 'N_01_TL', 'N_02_TL')]

  expect_equal(GetNeighborsPerBatch(seu, 'ID_sample', graph.name = 'SCT.nn', count.self = FALSE),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, seu[[]][,'ID_sample', drop = F], graph.name = 'SCT.nn', count.self = FALSE),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, seu[[]][, c('ID_sample', 'First_author')], graph.name = 'SCT.nn', count.self = FALSE),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, 'ID_sample', graph.name = 'SCT.nn_dist', count.self = FALSE),
               expected)
  expect_equal(GetNeighborsPerBatch(seu, 'ID_sample', graph.name = 'SCT_nn', count.self = FALSE),
               expected)
  expect_equal(GetPropInterBatch(seu, 'ID_sample', graph.name = 'SCT_nn', count.self = FALSE),
               naive_prop_inter(expected))
  expect_equal(GetPropInterBatch(seu, 'ID_sample', graph.name = 'SCT_nn', count.self = FALSE, per.batch = FALSE),
               naive_prop_inter(expected, per.batch = FALSE))
  expect_equal(GetPropIntraBatch(seu, 'ID_sample', graph.name = 'SCT_nn', count.self = FALSE),
               naive_prop_intra(expected))
  expect_equal(GetPropIntraBatch(seu, 'ID_sample', graph.name = 'SCT_nn', count.self = FALSE, per.batch = FALSE),
               naive_prop_intra(expected, per.batch = FALSE))


  expected[] <- c(843, 295, 46, 53, 295, 218, 67, 96, 46, 67, 575, 399, 53, 96, 399, 1558)
  expect_equal(GetNeighborsPerBatch(seu, 'ID_sample', graph.name = 'SCT.nn_symmetric'),
               expected)
  expect_equal(GetPropInterBatch(seu, 'ID_sample', graph.name = 'SCT.nn_symmetric'),
               naive_prop_inter(expected))
  expect_equal(GetPropInterBatch(seu, 'ID_sample', graph.name = 'SCT.nn_symmetric', per.batch = FALSE),
               naive_prop_inter(expected, per.batch = FALSE))
  expect_equal(GetPropIntraBatch(seu, 'ID_sample', graph.name = 'SCT.nn_symmetric'),
               naive_prop_intra(expected))
  expect_equal(GetPropIntraBatch(seu, 'ID_sample', graph.name = 'SCT.nn_symmetric', per.batch = FALSE),
               naive_prop_intra(expected, per.batch = FALSE))
  diag(expected) <- diag(expected) - table(seu$ID_sample)[c('B_13_TL', 'B_18_TL', 'N_01_TL', 'N_02_TL')]
  expect_equal(GetNeighborsPerBatch(seu, 'ID_sample', graph.name = 'SCT.nn_symmetric', count.self = FALSE),
               expected)
  expect_equal(GetPropInterBatch(seu, 'ID_sample', graph.name = 'SCT.nn_symmetric', count.self = FALSE),
               naive_prop_inter(expected))
  expect_equal(GetPropInterBatch(seu, 'ID_sample', graph.name = 'SCT.nn_symmetric', count.self = FALSE, per.batch = FALSE),
               naive_prop_inter(expected, per.batch = FALSE))
  expect_equal(GetPropIntraBatch(seu, 'ID_sample', graph.name = 'SCT.nn_symmetric', count.self = FALSE),
               naive_prop_intra(expected))
  expect_equal(GetPropIntraBatch(seu, 'ID_sample', graph.name = 'SCT.nn_symmetric', count.self = FALSE, per.batch = FALSE),
               naive_prop_intra(expected, per.batch = FALSE))
  expected[] <- c(1103, 523, 121, 140, 523, 364, 179, 270, 121, 179, 919, 969, 140, 270, 969, 2670)
  expect_equal(GetNeighborsPerBatch(seu, 'ID_sample', graph.name = 'SCT_snn'),
               expected)
  expect_equal(GetPropInterBatch(seu, 'ID_sample', graph.name = 'SCT_snn'),
               naive_prop_inter(expected))
  expect_equal(GetPropInterBatch(seu, 'ID_sample', graph.name = 'SCT_snn', per.batch = FALSE),
               naive_prop_inter(expected, per.batch = FALSE))
  expect_equal(GetPropIntraBatch(seu, 'ID_sample', graph.name = 'SCT_snn'),
               naive_prop_intra(expected))
  expect_equal(GetPropIntraBatch(seu, 'ID_sample', graph.name = 'SCT_snn', per.batch = FALSE),
               naive_prop_intra(expected, per.batch = FALSE))
  diag(expected) <- diag(expected) - table(seu$ID_sample)[c('B_13_TL', 'B_18_TL', 'N_01_TL', 'N_02_TL')]
  expect_equal(GetNeighborsPerBatch(seu, 'ID_sample', graph.name = 'SCT_snn', count.self = FALSE),
               expected)
  expect_equal(GetPropInterBatch(seu, 'ID_sample', graph.name = 'SCT_snn', count.self = FALSE),
               naive_prop_inter(expected))
  expect_equal(GetPropInterBatch(seu, 'ID_sample', graph.name = 'SCT_snn', count.self = FALSE, per.batch = FALSE),
               naive_prop_inter(expected, per.batch = FALSE))
  expect_equal(GetPropIntraBatch(seu, 'ID_sample', graph.name = 'SCT_snn', count.self = FALSE),
               naive_prop_intra(expected))
  expect_equal(GetPropIntraBatch(seu, 'ID_sample', graph.name = 'SCT_snn', count.self = FALSE, per.batch = FALSE),
               naive_prop_intra(expected, per.batch = FALSE))

})

test_that("assessement of whether k is constant is correct", {
  knn_dist <- matrix(rep(0, 12), ncol = 1)
  knn_idx <- matrix(1:12, ncol = 1)
  knn_object <- SeuratObject:::Neighbor(nn.idx = knn_idx, nn.dist = knn_dist,
                                        cell.names = paste0("Cell_", 1:12))
  expect_true(is.kconstant(knn_object))
  expect_true(is.kconstant(as.Graph(knn_object)))

  knn_idx <- cbind(knn_idx, matrix(c(2, 5, 6, 7, 8, 10, 12,
                                     1, 3, 4, 6, 8,  9, 10,
                                     2, 4, 5, 6, 7,  9, 12,
                                     2, 3, 6, 7, 9, 10, 11,
                                     1, 3, 6, 7, 8,  9, 11,
                                     1, 2, 3, 4, 5,  9, 11,
                                     1, 3, 4, 5, 6,  8, 11,
                                     1, 2, 5, 7, 9, 10, 12,
                                     2, 3, 4, 5, 6,  8, 10,
                                     1, 2, 4, 8, 9, 11, 12,
                                     2, 4, 5, 6, 7, 10, 12,
                                     1, 3, 5, 7, 8, 10, 11),
                                   nrow = 12, byrow = TRUE))
  knn_dist <- cbind(knn_dist, matrix(c(8.9, .1, .5, .7, .4, 2, 12,
                                       8.9, .0098, 4, 22, 6.3, 4, .5,
                                       .0098, 4.992, 10, 10, 2.987, 82, 47,
                                       4, 4.992, 7, 0.518, 2.091, 17.8, 9.2,
                                       .1, 10, 7, 5.0927, 40.17, 33, 7,
                                       .5, 22, 10, 7, 7, 3, 7,
                                       .7, 2.987, 0.518, 5.0927, 28.3, 2.987, 2.987,
                                       .4, 6.3, 40.17, 2.987, 2.4, 45.9, 2,
                                       4, 82, 2.091, 33, 3, 2.4, 26.1,
                                       2, .5, 17.8, 45.9, 26.1, 67, .8,
                                       72, 9.2, 7, 7, 2.987, 67, 58.789,
                                       12, 47, 86, 2, 2, .8, 58.789),
                                     nrow = 12, byrow = TRUE))

  knn_object <- SeuratObject:::Neighbor(nn.idx = knn_idx, nn.dist = knn_dist,
                                        cell.names = paste0("Cell_", 1:12))
  expect_true(is.kconstant(knn_object))
  g <- as.Graph(knn_object)
  expect_true(is.kconstant(g))
  g <- SymmetrizeKnn(g)
  expect_false(is.kconstant(g))
})

test_that("fetched values of k are accurate", {
  knn_dist <- matrix(rep(0, 12), ncol = 1)
  knn_idx <- matrix(1:12, ncol = 1)
  knn_object <- SeuratObject:::Neighbor(nn.idx = knn_idx, nn.dist = knn_dist,
                                        cell.names = paste0("Cell_", 1:12))
  expect_equal(get.k(knn_object, FUN = "min"), get.k(knn_object, FUN = "max"))
  expect_equal(get.k(knn_object, FUN = "min"), 1)
  expect_equal(get.k(knn_object, FUN = "range"), c(1, 1))
  expect_equal(get.k(knn_object, FUN = "all"), rep(1, 12))
  g <- as.Graph(knn_object)
  expect_equal(get.k(g, FUN = "min"), get.k(g, FUN = "max"))
  expect_equal(get.k(g, FUN = "min"), 1)
  expect_equal(get.k(g, FUN = "range"), c(1, 1))
  expect_equal(get.k(g, FUN = "all"), setNames(rep(1, 12), paste0("Cell_", 1:12)))

  seu <- PreProcessSO(liver_small, use_sct = TRUE, npcs = 10)
  for (k.param in c(1, 2, 20, 50, 100, 200)) {
    seu <- suppressMessages(FindNeighbors(seu, dims = 1:10, return.neighbor = TRUE, verbose = FALSE, k.param = k.param))
    expect_equal(get.k(seu[["SCT.nn"]], FUN = "min"), get.k(seu[["SCT.nn"]], FUN = "max"))
    expect_equal(get.k(seu[["SCT.nn"]], FUN = "min"), k.param)
    expect_equal(get.k(seu[["SCT.nn"]], FUN = "range"), c(k.param, k.param))
    expect_equal(get.k(seu[["SCT.nn"]], FUN = "all"), rep(k.param, 200))

    seu[['SCT_nn']] <- as.Graph(seu[['SCT.nn']])
    expect_equal(get.k(seu[["SCT_nn"]], FUN = "min"), get.k(seu[["SCT_nn"]], FUN = "max"))
    expect_equal(get.k(seu[["SCT_nn"]], FUN = "min"), k.param)
    expect_equal(get.k(seu[["SCT_nn"]], FUN = "range"), c(k.param, k.param))
    expect_equal(get.k(seu[["SCT_nn"]], FUN = "all"), setNames(rep(k.param, 200), Cells(seu)))
  }

})

test_that("the determination of whether this could be a connectivity graph works", {
  for (v in c(-Inf, -2, -1, 0, .5, 1, 12, Inf)) {
    knn_dist <- matrix(rep(v, 12), ncol = 1)
    knn_idx <- matrix(1:12, ncol = 1)
    knn_object <- SeuratObject:::Neighbor(nn.idx = knn_idx, nn.dist = knn_dist,
                                          cell.names = paste0("Cell_", 1:12))
    if (v <= 1 & v >= 0) {
      expect_true(could.be.connectivity(knn_object))
      expect_true(could.be.connectivity(as.Graph(knn_object)))
    } else {
      expect_false(could.be.connectivity(knn_object))
      expect_false(could.be.connectivity(as.Graph(knn_object)))
    }
  }

  knn_dist <- cbind(matrix(rep(0, 12), ncol = 1),
                    matrix(c(2,.7, 3.5,.982, 5, 2.3, .009, .2, 1.057,2.858,3.7,
                             2,1.780,6,5.16,.285,8.406,1.3,4,4.156,6.086,2.32,9,
                             3.918,2.344,1.067,4.453,1.893,1.059,1.653,.605,.37,
                             1.01,5,8.741, 10.43), ncol = 3))
  knn_idx <- cbind(matrix(1:12, ncol = 1),
                   matrix(c(3,10,4,11,10,8,10,6,2,11,10,6,3,1,8,10,11,12,3,6,5,
                            3,4,1,6,10,12,5,12,3,5,4,6,6,3,9), ncol = 3, byrow = TRUE))
  knn_object <- SeuratObject:::Neighbor(nn.idx = knn_idx, nn.dist = knn_dist,
                                        cell.names = paste0("Cell_", 1:12))

  expect_false(could.be.connectivity(knn_object))
  expect_false(could.be.connectivity(as.Graph(knn_object)))

  knn_dist <- rescale(knn_dist)
  knn_object_bounded <- SeuratObject:::Neighbor(
    nn.idx = knn_idx, nn.dist = knn_dist, cell.names = paste0("Cell_", 1:12))

  expect_false(could.be.connectivity(knn_object_bounded))
  expect_false(could.be.connectivity(as.Graph(knn_object_bounded)))
  expect_true(could.be.connectivity(knn_object_bounded, check.symmetry = FALSE))
  expect_true(could.be.connectivity(as.Graph(knn_object_bounded), check.symmetry = FALSE))

  expect_true(could.be.connectivity(SymmetrizeKnn(knn_object_bounded)))
  expect_true(could.be.connectivity(SymmetrizeKnn(as.Graph(knn_object_bounded))))

  expect_true(could.be.connectivity(compute.umap.connectivities(knn_object@nn.dist, knn_object@nn.idx, sorted.dist = F, verbose = FALSE)))
})


test_that("Seurat-inspired functions' outputs are identical", {
  m <- matrix(nrow=0, ncol=0)
  expect_equal(n_zeros_mat(m), 0)
  expect_equal(n_zeros_mat(as.dgcmatrix(m)), 0)

  m <- matrix(0, ncol = 20, nrow = 100)
  expect_equal(n_zeros_mat(m), 2000)
  expect_equal(n_zeros_mat(as.dgcmatrix(m)), 2000)

  m <- matrix(runif(2000, min = -10, max = 10), ncol = 20, nrow = 100)
  nz <- sum(m == 0)
  expect_equal(n_zeros_mat(m), nz)
  expect_equal(n_zeros_mat(as.dgcmatrix(m)), nz)
  m <- round(m, 0)
  nz <- sum(m == 0)
  expect_equal(n_zeros_mat(m), nz)
  expect_equal(n_zeros_mat(as.dgcmatrix(m)), nz)

  m[sample(seq_len(length(m)), 20)] <- NA
  nz <- sum(m == 0, na.rm = TRUE)
  expect_equal(n_zeros_mat(m), nz)
  expect_equal(n_zeros_mat(as.dgcmatrix(m)), nz)

  possibilities <- expand.grid(seq_len(nrow(m)), seq_len(ncol(m)))
  idx <- possibilities[sample(seq_len(nrow(possibilities)), 200),]
  m <- sparseMatrix(i = idx$Var1, j = idx$Var2,
                    x = sample(0:5, nrow(idx), replace = TRUE,
                               prob = c(10,1,1,1,1,1)),
                    dims = dim(m), dimnames = dimnames(m))
  nz <- sum(m == 0, na.rm = TRUE)
  expect_equal(n_zeros_mat(m), nz)
  m <- sparseMatrix(i = idx$Var1, j = idx$Var2,
                    x = sample(c(0:5, NA_integer_), nrow(idx), replace = TRUE,
                               prob = c(10,1,1,1,1,1,3)),
                    dims = dim(m), dimnames = dimnames(m))
  nz <- sum(m == 0, na.rm = TRUE)
  expect_equal(n_zeros_mat(m), nz)
  expect_equal(n_zeros_mat(as.data.frame(as.matrix(m))), nz)
  expect_nan(n_zeros_mat(liver_small))
#
#   m <- as.dgcmatrix(m)
#   i <- slot(object = m, name = "i") + 1
#   x <- slot(object = m, name = "x")
#   j <- findInterval(seq(x)-1, slot(object = m, name = "p")[-1]) + 1
#   i <- c(i, sample(seq_len(nrow(m)), 20))
#   j <- c(j, sample(seq_len(ncol(m)), 20))
#   x <- c(x, rep(0, 20))
})

test_that("matrix format selection is correct", {
  m <- matrix(nrow=0, ncol=0)
  expect_identical(choose_matrix_format(m), m)
  expect_identical(choose_matrix_format(as.dgcmatrix(m)), as.dgcmatrix(m))

  m <- matrix(nrow=10, ncol=10, dimnames = list(paste0("A", 1:10),
                                                paste0("A", 1:10)))
  expect_identical(choose_matrix_format(m), m)
  expect_true(is.matrix(choose_matrix_format(as.dgcmatrix(m))))
  expect_true(all(is.na(choose_matrix_format(as.dgcmatrix(m)))))
  expect_identical(dimnames(choose_matrix_format(m)), dimnames(m))


  possibilities <- expand.grid(seq_len(5e2), seq_len(5e2))
  idx <- possibilities[sample(seq_len(nrow(possibilities)), ceiling(nrow(possibilities) * .345)),]
  m <- sparseMatrix(i = idx$Var1, j = idx$Var2,
                    x = rnorm(nrow(idx), mean = 3, sd = 12),
                    dims = c(5e2, 5e2), dimnames = list(paste0("A", 1:5e2),
                                                        paste0("A", 1:5e2)))
  expect_identical(choose_matrix_format(m), m)
  expect_identical(choose_matrix_format(as.dgcmatrix(m)), m)

  idx <- possibilities[sample(seq_len(nrow(possibilities)), floor(nrow(possibilities) * .351)),]
  m <- sparseMatrix(i = idx$Var1, j = idx$Var2,
                    x = rnorm(nrow(idx), mean = 3, sd = 12),
                    dims = c(5e2, 5e2), dimnames = list(paste0("A", 1:5e2),
                                                        paste0("A", 1:5e2)))
  expect_identical(choose_matrix_format(m), m)
  expect_identical(choose_matrix_format(as.matrix(m)), m)


  m <- sparseMatrix(i = sample(seq_len(5e4), 1e5, T),
                    j = sample(seq_len(5e4), 1e5, T),
                    x = 1,
                    dims = c(5e4, 5e4), dimnames = list(paste0("A", 1:5e4),
                                                        paste0("A", 1:5e4)))
  expect_identical(choose_matrix_format(m), m)

  max_int_32bit <- 2^31 - 1
  l <- 1e10
  n_0s <- 4e9
  n_not0s <- l - n_0s
  expect_true(n_not0s > max_int_32bit)
  expect_false(n_0s/l < .35)
  expect_true((n_not0s > max_int_32bit) | (n_0s/l < .35))

  l <- 1e12
  n_0s <- 5e9
  n_not0s <- l - n_0s
  expect_true(n_not0s > max_int_32bit)
  expect_true(n_0s/l < .35)
  expect_true((n_not0s > max_int_32bit) | (n_0s/l < .35))

})

test_that("Seurat-inspired functions' outputs are identical", {
  seu <- PreProcessSO(liver_small, split = 'First_author', use_sct = FALSE)
  expect_equal(CreateIntegrationGroups(seu[['RNA']], Layers(seu, search = "counts"), "scale.data"),
               Seurat:::CreateIntegrationGroups(seu[['RNA']], Layers(seu, search = "counts"), "scale.data"))
  seu <-  PreProcessSO(seu, use_sct = TRUE)
  expect_equal(CreateIntegrationGroups(seu[['SCT']], Layers(seu, search = "counts"), "scale.data"),
               Seurat:::CreateIntegrationGroups(seu[['SCT']], Layers(seu, search = "counts"), "scale.data"))

  seu <- PreProcessSO(liver_small, split = 'ID_sample', use_sct = FALSE)
  expect_equal(CreateIntegrationGroups(seu[['RNA']], Layers(seu[['RNA']], search = "counts"), "scale.data"),
               Seurat:::CreateIntegrationGroups(seu[['RNA']], Layers(seu[['RNA']], search = "counts"), "scale.data"))
  expect_equal(CreateIntegrationGroups(seu[['RNA']], Layers(seu[['RNA']], search = "counts")[-2], "scale.data"),
               Seurat:::CreateIntegrationGroups(seu[['RNA']], Layers(seu[['RNA']], search = "counts")[-2], "scale.data"))
  cmap <- cmap_backup <- unclass(seu@assays$RNA@cells)
  cmap[sample(sample(seq_len(nrow(cmap)), 30)), 1:(ncol(cmap) - 1)] <-
    sample(c(TRUE, FALSE), 30*(ncol(cmap) - 1), replace = TRUE, prob = c(1,.3))
  seu@assays$RNA@cells <- as(cmap, "LogMap")
  expect_equal(CreateIntegrationGroups(seu[['RNA']], Layers(seu[['RNA']], search = "counts"), "scale.data"),
               Seurat:::CreateIntegrationGroups(seu[['RNA']], Layers(seu[['RNA']], search = "counts"), "scale.data"))
  seu@assays$RNA@cells <- as(cmap_backup, "LogMap")
  seu <-  PreProcessSO(seu, use_sct = TRUE)
  expect_equal(CreateIntegrationGroups(seu[['SCT']], Layers(seu[['SCT']], search = "counts"), "scale.data"),
               Seurat:::CreateIntegrationGroups(seu[['SCT']], Layers(seu[['SCT']], search = "counts"), "scale.data"))

  seu <- PreProcessSO(liver_small, split = 'orig.ident', use_sct = FALSE)
  expect_error(CreateIntegrationGroups(seu[['RNA']], Layers(seu[['RNA']], search = "counts"), "scale.data"),
               'attempt to set an attribute on NULL')
  seu <-  PreProcessSO(seu, use_sct = TRUE)
  expect_equal(CreateIntegrationGroups(seu[['SCT']], Layers(seu[['SCT']], search = "counts"), "scale.data"),
               Seurat:::CreateIntegrationGroups(seu[['SCT']], Layers(seu[['SCT']], search = "counts"), "scale.data"))

  nrow <- sample(50:200, 1)
  ncol <- sample(100:400, 1)

  m <- matrix(rnorm(nrow*ncol, 3, sd = 2), nrow = nrow, ncol = ncol)
  m[m < 0] <- 0

  m <- as.matrix(dist(m))
  expect_equal(NormaliseL2(m, MARGIN = 1), Seurat:::L2Norm(m, MARGIN = 1))
  expect_equal(NormaliseL2(m, MARGIN = 2), Seurat:::L2Norm(m, MARGIN = 2))

  m <- as.dgcmatrix(m)
  expect_equal(NormaliseL2(m, MARGIN = 1), Seurat:::L2Norm(m, MARGIN = 1))
  expect_equal(NormaliseL2(m, MARGIN = 2), Seurat:::L2Norm(m, MARGIN = 2))

  expect_equal(NormaliseL2(seu[["RNA"]]@layers$scale.data, 1),
               Seurat:::L2Norm(seu[["RNA"]]@layers$scale.data, 1))
  expect_equal(NormaliseL2(seu[["RNA"]]@layers$scale.data, 2),
               Seurat:::L2Norm(seu[["RNA"]]@layers$scale.data, 2))
})

test_that("converting a Graph to a Graph preserves everything", {
  knn_dist <- matrix(rep(0, 12), ncol = 1)
  knn_idx <- matrix(1:12, ncol = 1)
  knn_object <- as.Graph(SeuratObject:::Neighbor(
    nn.idx = knn_idx, nn.dist = knn_dist, cell.names = paste0("Cell_", 1:12)))
  expect_equal(knn_object, as.Graph(knn_object))
  slot(knn_object, "assay.used") <- "RNA"
  expect_equal(knn_object, as.Graph(knn_object))

  knn_idx <- cbind(knn_idx, matrix(c(2, 5, 6, 7, 8, 10, 12,
                                     1, 3, 4, 6, 8,  9, 10,
                                     2, 4, 5, 6, 7,  9, 12,
                                     2, 3, 6, 7, 9, 10, 11,
                                     1, 3, 6, 7, 8,  9, 11,
                                     1, 2, 3, 4, 5,  9, 11,
                                     1, 3, 4, 5, 6,  8, 11,
                                     1, 2, 5, 7, 9, 10, 12,
                                     2, 3, 4, 5, 6,  8, 10,
                                     1, 2, 4, 8, 9, 11, 12,
                                     2, 4, 5, 6, 7, 10, 12,
                                     1, 3, 5, 7, 8, 10, 11),
                                   nrow = 12, byrow = TRUE))
  knn_dist <- cbind(knn_dist, matrix(c(8.9, .1, .5, .7, .4, 2, 12,
                                       8.9, .0098, 4, 22, 6.3, 4, .5,
                                       .0098, 4.992, 10, 10, 2.987, 82, 47,
                                       4, 4.992, 7, 0.518, 2.091, 17.8, 9.2,
                                       .1, 10, 7, 5.0927, 40.17, 33, 7,
                                       .5, 22, 10, 7, 7, 3, 7,
                                       .7, 2.987, 0.518, 5.0927, 28.3, 2.987, 2.987,
                                       .4, 6.3, 40.17, 2.987, 2.4, 45.9, 2,
                                       4, 82, 2.091, 33, 3, 2.4, 26.1,
                                       2, .5, 17.8, 45.9, 26.1, 67, .8,
                                       72, 9.2, 7, 7, 2.987, 67, 58.789,
                                       12, 47, 86, 2, 2, .8, 58.789),
                                     nrow = 12, byrow = TRUE))
  knn_object <- as.Graph(SeuratObject:::Neighbor(
    nn.idx = knn_idx, nn.dist = knn_dist, cell.names = paste0("Cell_", 1:12)))
  expect_equal(knn_object, as.Graph(knn_object))
  slot(knn_object, "assay.used") <- "RNA"
  expect_equal(knn_object, as.Graph(knn_object))
})

