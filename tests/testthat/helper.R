PreProcessSO <- function(obj, split = NULL, use_sct = FALSE, nfeatures = NULL,
                          npcs = NULL, approx = TRUE) {
  split <- split %||% FALSE
  use_sct <- use_sct %||% FALSE
  npcs <- npcs %||% 0
  nfeatures <- nfeatures %||% ifelse(use_sct, 3e3, 2e3)
  withr::local_options(list(lifecycle_verbosity = 'quiet'))

  if (! isFALSE(split)) {
    obj <- split(obj, obj[[]][[split]])
  }

  if (use_sct) {
    obj <- SCTransform(obj, variable.features.n = nfeatures, verbose = F)
  } else {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = nfeatures, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
  }

  if (npcs > 0) {
    obj <- RunPCA(obj, npcs = npcs, approx = approx, verbose = FALSE)
  }
  return(obj)
}

naive_row_sort <- function(mat, by = NULL, decreasing = FALSE) {
  if (is.null(by)) {
    t(apply(mat, 1, sort, decreasing = decreasing))
  } else {
    t(sapply(1:nrow(mat), function(i) mat[i, order(by[i,], decreasing = decreasing)]))
  }
}
naive_col_sort <- function(mat, by = NULL, decreasing = FALSE) {
  if (is.null(by)) {
    apply(mat, 2, sort, decreasing = decreasing)
  } else {
    sapply(1:ncol(mat), function(i) mat[order(by[,i], decreasing = decreasing), i])
  }
}

naive_knn_sort <- function(obj, cut_k = NULL) {
  n <- nrow(obj@nn.idx)
  obj_k <- ncol(obj@nn.idx)
  cut_k <- min(cut_k %||% obj_k, obj_k)
  if (obj_k > 1) {
    oobj <- SeuratObject:::Neighbor(
      nn.idx = t(sapply(1:n, function(i) {
        obj@nn.idx[i,order(obj@nn.dist[i,])]
      }))[,seq_len(cut_k), drop = FALSE],
      nn.dist = t(sapply(1:n, function(i) {
        obj@nn.dist[i,order(obj@nn.dist[i,])]
      }))[,seq_len(cut_k), drop = FALSE],
      alg.idx = obj@alg.idx,
      alg.info = obj@alg.info,
      cell.names = obj@cell.names)
  } else {
    oobj <- obj
  }
  return(oobj)
}

naive_positive_matrix_symmetry <- function(mat, use_max = TRUE, assay.used = NULL) {
  stopifnot("The matrix is not square" = ncol(mat) == nrow(mat),
            "The matrix is not finite" = all(is.finite(mat)))
  n <- ncol(mat)
  mat[mat == 0] <- c(Inf, -Inf)[use_max + 1]
  `%<|>%` <- c(`>`, `<`)[[use_max + 1]]
  for (i in seq_len(n - 1)) {
    for (j in seq(from = i + 1, to = n)) {
      if (mat[i,j] %<|>% mat[j,i]) {
        mat[i,j] <- mat[j,i]
      } else {
        mat[j,i] <- mat[i,j]
      }
    }
  }
  mat[is.infinite(mat)] <- 0
  mat <- assay.used %iff% as.Graph(mat) %||% mat
  assay.used %iff% {slot(mat, 'assay.used') <- assay.used}
  return(mat)
}

naive_get_closest <- function(mat, k.max, decreasing = FALSE) {
  res <- apply(mat, 1, function(x) {
    o <- order(x, decreasing = decreasing)
    if (decreasing) {
      o <- o[1:(sum(x != 0))]
    } else {
      o <- o[(sum(x == 0)+1):length(x)]
    }

    as.numeric(o[1:k.max])
  })
  matrix(res, ncol = k.max, byrow = TRUE)
}

naive_prune_graph <- function(obj, k_value) {
  i <- slot(object = obj, name = "i") + 1
  x <- slot(object = obj, name = "x")
  p <- slot(object = obj, name = "p")
  j <- findInterval(seq(x)-1,p[-1]) + 1

  new_x <- new_j <- new_i <- numeric(0)
  for (r in seq_len(nrow(obj))) {
    x_ <- x[i == r]
    j_ <- j[i == r]
    i_ <- i[i == r]
    o <- order(x_)[1:min(k_value, length(x_))]
    new_x <- c(new_x, x_[o])
    new_j <- c(new_j, j_[o])
    new_i <- c(new_i, i_[o])
  }

  o <- order(new_j,new_i)
  obj@i <- as.integer(new_i[o] - 1)
  obj@x <- new_x[o]
  obj@p <- as.integer(cumsum(c(0, table(factor(new_j[o], 1:nrow(obj))))))
  obj
}

naive_knn_per_batch <- function(knn_obj, batch_var, batch_name = "batch") {
  batches <- levels(batch_var) %||% unique(batch_var)
  n_batch <- length(batches)
  res <- matrix(0, nrow = n_batch, ncol = n_batch, dimnames =
                list(sort(paste0(batch_name, batches)),
                     sort(paste0(batch_name, batches))))
  i_max <- nrow(knn_obj@nn.idx)
  j_max <- ncol(knn_obj@nn.idx)
  if (j_max < 2) {
    return(res)
  }
  for (i in seq_len(i_max)) {
    i_ <- paste0(batch_name, batch_var[i])
    for (j in seq_len(j_max)) {
      j_ <- paste0(batch_name, batch_var[knn_obj@nn.idx[i, j]])
      res[i_,j_] <- res[i_,j_] + 1
    }
  }
  return(res)
}

naive_prop_intra <- function(mat, per.batch = TRUE) {
  if (per.batch) {
    diag(mat) / rowSums(mat)
  } else {
    sum(diag(mat))/sum(mat)
  }
}
naive_prop_inter <- function(mat, per.batch = TRUE) {
  n_tot <- rowSums(mat)
  n_inter <- rowSums(`diag<-`(mat, 0))
  if (per.batch) {
    n_inter / n_tot
  } else {
    sum(n_inter) / sum(n_tot)
  }
}
