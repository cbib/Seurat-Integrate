#' Run bbknn on Seurat's \link[SeuratObject]{Assay5} object through \code{\link[Seurat]{IntegrateLayers}}
#'
#' @description
#' A wrapper to run \code{bbknn} on multi-layered Seurat V5 object.
#' Requires a conda environment with \code{bbknn} and necessary dependencies
#'
#' @inheritParams integration-method
#' @param conda_env Path to conda environment to run bbknn (should also
#' contain the scipy python module).  By default, uses the conda environment
#' registered for bbknn in the conda environment manager
#' @param new.graph Name of the \link[SeuratObject:Graph-class]{Graph object}
#' @param reconstructed.assay Name for the \code{assay} containing the corrected
#' expression matrix
#' @param ndims Number of dimensions for the new PCA computed on first output of
#' bbknn. 50 by default. Ignored when \code{ridge_regression = FALSE}
#' @param ndims.use Number of dimensions from \code{orig} to use for bbknn, and
#' from newly computed PCA when \code{ridge_regression = TRUE}.
#' @param ridge_regression When set to \code{TRUE} (default), new clusters are
#' computed on the output of bbknn, then a ridge regression is performed to
#' remove technical variables while preserving biological variables. Then, a new
#' bbknn run is performed.
#' @param graph.use Which graph(s) of bbknn to output. At least one of
#' "\code{connectivities}" or "\code{distances}". If both are provided (default)
#' and \code{ridge_regression = TRUE}, the first one ("\code{connectivities}" by
#' default, recommended) is used as input for computing clusters.
#' @param ... Additional arguments to be passed to \code{bbknn.bbknn()}. When
#' \code{ridge_regression = TRUE}, also accepts arguments to pass to
#' \code{Seurat::FindClusters()}, \code{Seurat::RunPCA()} and
#' \code{bbknn.ridge_regression()}. See \strong{Details} section
#'
#' @return A list containing at least one of:
#' \itemize{
#'   \item 1 or 2 new Graph(s) of name
#'   [\code{new_graph}]_scale.data_[\code{graph.use}]
#'   corresponding to the output(s) of the first run of \pkg{bbknn}
#'   \item a new Assay of name \code{reconstructed.assay} with corrected counts
#'   for each feature from \code{scale.layer}.
#'   \item a new DimReduc (PCA) of name \code{new.reduction} (key set to
#'   \code{reduction.key})
#'   \item 1 or 2 new Graph(s) of name
#'   [\code{new_graph}]_ridge.residuals_[\code{graph.use}]
#'   corresponding to the output(s) of the second run of \pkg{bbknn}
#' }
#'
#' [\code{graph.use}] can take two values (either "\code{connectivities}" or
#' "\code{distances}"), depending on the \code{graph.use} parameter.
#'
#' When called via \code{\link[Seurat]{IntegrateLayers}}, a Seurat object with
#' the new reduction and/or assay is returned
#'
#' @details
#' This wrappers calls three python functions through \pkg{reticulate}.
#' Find the \pkg{bbknn}-specific arguments there:
#' \itemize{
#'   \item bbknn function:
#'   \href{https://bbknn.readthedocs.io/en/latest/bbknn.bbknn.html}{bbknn.bbknn}, which relies on
#'   \href{https://bbknn.readthedocs.io/en/latest/bbknn.matrix.bbknn.html}{bbknn.matrix.bbknn}
#'   \item ridge regression:
#'   \href{https://bbknn.readthedocs.io/en/latest/bbknn.ridge_regression.html}{bbknn.ridge_regression}, which relies on
#'   \href{https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Ridge.html}{sklearn.linear_model.Ridge}
#' }
#'
#' @importFrom reticulate use_condaenv import r_to_py py_to_r
#' @importFrom Matrix t
#' @importFrom Seurat CreateDimReducObject FindClusters
#' @importFrom SeuratObject JoinLayers GetAssayData
#'
#' @export
#'@note This function requires the
#' \href{https://bbknn.readthedocs.io/en/latest/index.html}{\pkg{bbknn}} package
#' to be installed (along with \pkg{scipy})
#'
#' @examples
#' \dontrun{
#' # Preprocessing
#' obj <- SeuratData::LoadData("pbmcsca")
#' obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj)
#'
#' # After preprocessing, we integrate layers:
#' obj <- IntegrateLayers(object = obj, method = bbknnIntegration,
#'                        conda_env = 'bbknn', groups = obj[[]],
#'                        groups.name = 'Method')
#'
#' # To disable the ridge regression and subsequent steps:
#' obj <- IntegrateLayers(object = obj, method = bbknnIntegration,
#'                        conda_env = 'bbknn', groups = obj[[]],
#'                        groups.name = 'Method', ridge_regression = FALSE)
#' }
#'
#' @references Polański, K., Young, M. D., Miao, Z., Meyer, K. B., Teichmann,
#' S. A. & Park, J.-E. BBKNN: fast batch alignment of single cell transcriptomes.
#' Bioinformatics 36, 964–965 (2019).
#' \href{https://doi.org/10.1093/bioinformatics/btz625}{DOI}
#'
#' @seealso \code{\link[Seurat]{IntegrateLayers}}, \code{\link[Seurat]{writing-integration}}

bbknnIntegration <- function(
    object,
    orig,
    groups = NULL,
    groups.name = NULL,
    layers = 'data',
    scale.layer = 'scale.data',
    conda_env = NULL,
    new.graph = 'bbknn',
    new.reduction = "pca.bbknn",
    reduction.key = "bbknnPCA_",
    reconstructed.assay = 'bbknn.ridge',
    ndims = 50L,
    ndims.use = 30L,
    ridge_regression = T,
    graph.use = c("connectivities", "distances"),
    verbose = TRUE,
    seed.use = 42L,
    ...
) {
  graph.use <- match.arg(graph.use, several.ok = TRUE)
  seed.use <- seed.use %iff% as.integer(seed.use)

  args.bbknn <- c("trim", "annoy_n_trees", "pynndescent_n_neighbors",
                  "metric", "set_op_mix_ratio",
                  # "pynndescent_random_state", "metric", "set_op_mix_ratio",
                  "local_connectivity")
  args.ridge <- c("chunksize")
  args.findCluster <- c("modularity.fxn", "initial.membership", "node.sizes",
                        "resolution", "method", "algoritm", "n.start", "n.iter",
                        "group.singletons", "temp.file.location",
                        "edge.file.name")
  args.runPCA <- c("weight.by.var", "approx")

  varargs <- list(...)

  conda_bin <- "auto"
  if (is.null(conda_env) || is.na(conda_env) || isFALSE(conda_env)) {
    if (! isValid(conda_status$current[["bbknn"]], do.check = TRUE)) {
      abort(message = paste("bbknn conda environment is not valid. Either",
                            "set", sQuote("conda_env"), "argument or create",
                            "the environment via the conda manager"))
    }
    message("Using conda from conda environment manager\n"[verbose], appendLF = FALSE)
    conda_env <- conda_status$current[["bbknn"]][["conda.env.path"]]$value
    conda_bin <- conda_status$current[["bbknn"]][["conda.bin"]]$value
  }

  use_condaenv(conda_env, conda = conda_bin, required = TRUE)
  anndata <- reticulate::import("anndata",convert=FALSE)
  bbknn <- reticulate::import("bbknn", convert=FALSE)
  sc <- reticulate::import("scanpy",convert=FALSE)
  scipy <-  reticulate::import('scipy', convert = FALSE)

  if (!inherits(x = orig, what = 'DimReduc')) {
    abort(message = "'orig' must be a dimensional reduction")
  }
  new.graph <- new.graph %||% 'bbknn'
  new.reduction <- new.reduction %||% FALSE
  return.ridge.reduc <- ridge_regression & (! isFALSE(new.reduction))
  reconstructed.assay <- reconstructed.assay %||% FALSE
  return.new.assay <- ridge_regression & (! isFALSE(reconstructed.assay))

  X_pca <- Embeddings(object = orig)
  ndims <- as.integer(ndims %||% 50L)
  ndims.use <- as.integer(ndims.use %||% ncol(x = X_pca))
  if (ndims.use > ncol(x = X_pca)) {
    warning(sprintf(fmt =
                      "'ndims.use' must not be greater than number of dimensions in 'orig'. Setting to %d",
                    ncol(x = X_pca)), call. = T)
    ndims.use <- ncol(x = X_pca)
  }
  if (ndims < ndims.use) {
    warning(sprintf(fmt =
                      "'ndims' must not be greater than 'ndims.use'. Setting to %d",
                    ndims.use), call. = T)
    ndims <- ndims.use
  }
  layers <- Layers(object = object, search = layers %||% "data")
  scale.layer <- scale.layer %||% "scale.data"
  if (! scale.layer %in% Layers(object = object, search = scale.layer)) {
    abort(message = paste(sQuote(x = scale.layer), "not in object layers"))
  }

  groups <- groups %||% CreateIntegrationGroups(object = object,
                                                layers = layers,
                                                scale.layer = scale.layer)
  groups.name <- groups.name %||% colnames(groups)[1]
  groups.name <- intersect(colnames(groups), groups.name)
  if (! length(x = groups.name)) {
    abort(message = "'groups.name' not in 'groups' data frame")
  }
  if (length(x = groups.name) > 1) {
    groups.name <- groups.name
    warning(paste("more 'groups.name' that expected. Using the first one",
                  sQuote(x = groups.name)), call. = FALSE, immediate. = TRUE)
  }
  message("Preparing adata object..."[verbose], appendLF = FALSE)
  scaled.mat <- GetAssayData(object = object, layer = scale.layer)
  features <- rownames(x = scaled.mat)
  adata <- sc$AnnData(
    X   = scipy$sparse$csr_matrix(
      Matrix::t(x = scaled.mat)
    ),
    obs = reticulate::r_to_py(x = groups[, groups.name, drop = FALSE]),
    var = r_to_py(data.frame(features = features))
  )
  adata$X = adata$X$toarray()
  adata$obsm["X_pca"] = reticulate::r_to_py(x = X_pca)
  message("done.\n"[verbose], appendLF = FALSE)

  message("Running 1st instance of bbknn..."[verbose], appendLF = FALSE)
  args <- c(list(adata = adata, batch_key = groups.name, use_rep='X_pca',
                 key_added = NULL, copy = FALSE, n_pcs = ndims.use,
                 pynndescent_random_state = r_to_py(
                   varargs[["pynndescent_random_state"]] %||% seed.use)),
            varargs[intersect(names(varargs), args.bbknn)])
  do.call(bbknn$bbknn, args)
  message("done.\n"[verbose], appendLF = FALSE)

  bbknn_graphs <- sapply(graph.use, function(graph.key) {
    bbknn_graph <- reticulate::py_to_r(x = adata$obsp[graph.key])
    dimnames(bbknn_graph) <- list(
      as.character(reticulate::py_to_r(x = adata$obs_names$values)),
      as.character(reticulate::py_to_r(x = adata$obs_names$values))
    )
    bbknn_graph <- as.Graph(x = bbknn_graph)
    bbknn_graph@assay.used <- "RNA"
    bbknn_graph
  }, simplify = FALSE)

  graph.name <- sprintf(fmt = "%s_scale.data_%s", new.graph, graph.use)
  output.list <- setNames(bbknn_graphs, graph.name)

  if (ridge_regression) {
    message("Computing clusters..."[verbose], appendLF = FALSE)
    args <- c(list(object = bbknn_graphs[[1]], verbose = verbose,
                   random.seed = seed.use),
              varargs[intersect(names(varargs), args.findCluster)])
    clusters4ridge <- do.call(FindClusters, args)
    message("done.\n"[verbose], appendLF = FALSE)

    message("Computing ridge regression..."[verbose], appendLF = FALSE)
    adata$obs['clusters4ridge'] <- reticulate::r_to_py(x = clusters4ridge)
    args <- c(list(adata = adata,  batch_key = groups.name,
                   confounder_key = 'clusters4ridge',
                   random_state = r_to_py(seed.use)),
              varargs[intersect(names(varargs), args.ridge)])
    do.call(bbknn$ridge_regression, args)
    message("done.\n"[verbose], appendLF = FALSE)

    bbknn_assay <- object
    bbknn_assay@key <- "bbknnridge_"
    new.data <- t(x = reticulate::py_to_r(x = adata$X))
    rownames(x = new.data) <- features
    bbknn_assay <- SetAssayData(object = bbknn_assay, layer = scale.layer,
                                new.data = new.data)
    message("Computing PCA..."[verbose], appendLF = FALSE)
    args <- c(list(object = bbknn_assay, assay = reconstructed.assay,
                   reduction.name = new.reduction, reduction.key = reduction.key,
                   npcs = ndims, verbose = FALSE, seed.use = seed.use),
              varargs[intersect(names(varargs), args.runPCA)])
    bbknn_ridgePCA <- do.call(RunPCA, args)
    message("done.\n"[verbose], appendLF = FALSE)

    adata$obsm['X_pca'] = reticulate::r_to_py(
      x = Embeddings(object = bbknn_ridgePCA))

    message("Running 2nd instance of bbknn..."[verbose], appendLF = FALSE)
    args <- c(list(adata = adata, batch_key = groups.name, use_rep='X_pca',
                   key_added = NULL, copy = FALSE, n_pcs = ndims.use,
                   pynndescent_random_state = r_to_py(
                     varargs[["pynndescent_random_state"]] %||% seed.use)),
              varargs[intersect(names(varargs), args.bbknn)])
    do.call(bbknn$bbknn, args)
    message("done.\n"[verbose], appendLF = FALSE)

    bbknn_graphs <- sapply(graph.use, function(graph.key) {
      bbknn_graph <- reticulate::py_to_r(x = adata$obsp[graph.key])
      dimnames(bbknn_graph) <- list(
        as.character(reticulate::py_to_r(x = adata$obs_names$values)),
        as.character(reticulate::py_to_r(x = adata$obs_names$values))
      )
      bbknn_graph <- as.Graph(x = bbknn_graph)
      bbknn_graph@assay.used <- "RNA"
      bbknn_graph
    }, simplify = FALSE)
    if (return.new.assay) {
      output.list[[reconstructed.assay]] <- bbknn_assay
    }
    if (return.ridge.reduc) {
      output.list[[new.reduction]] <- bbknn_ridgePCA
    }
    graph.name <- sprintf(fmt = "%s_ridge.residuals_%s", new.graph, graph.use)
    for (i in 1:length(graph.name)) {
      output.list[[graph.name[i]]] <- bbknn_graphs[[i]]
    }
  }

  return(output.list)
}

attr(x = bbknnIntegration, which = 'Seurat.method') <- 'integration'
