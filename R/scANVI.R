#' @include kwargs.R
#'
NULL

#' Run scANVI on Seurat's \link[SeuratObject]{Assay5} object through \code{\link[Seurat]{IntegrateLayers}}
#'
#' @description
#' A wrapper to run \code{scANVI} on multi-layered Seurat V5 object.
#' Requires a conda environment with \code{scvi-tools} and necessary dependencies
#'
#' \strong{Recommendations}: use raw counts and all features
#' (\code{features = Features(object), layers = "counts"})
#'
#' @inheritParams integration-method
#' @param groups A \bold{named} data frame with grouping information. Can also
#' contain cell labels to guide scANVI.
#' @param labels.name Column name from \code{groups} data frame that stores
#' cell label information. If \code{labels.name = NULL}, all cells are assigned
#' the same label.
#' @param labels.null One value of \code{groups$labels.name} that indicates unlabeled
#' observations. \code{labels.null = NULL} means all labels are valid. Only applies
#' when \code{labels.name != NULL}.
#' @param layers Name of the layers to use in the integration.
#' \bold{'counts'} is highly recommended
#' @param conda_env Path to conda environment to run scANVI (should also
#' contain the scipy python module).  By default, uses the conda environment
#' registered for scANVI in the conda environment manager
#' @param torch.intraop.threads Number of intra-op threads available to torch
#' when training on CPU instead of GPU. Set via \code{torch.set_num_threads()}.
#' @param torch.interop.threads Number of intra-op threads available to torch
#' when training on CPU instead of GPU. Set via \code{torch.set_num_interop_threads()}.
#' Can only be changed once, on first call.
#' @param model.save.dir Path to a directory to save the model to. Uses
#' \code{SCANVI.save()}. Does not save anndata. Note that neither the trainer
#' optimizer state nor the trainer history are saved.
#' \code{model.save.dir = NULL} (default) disables saving the model.
#' @param ndims.out Number of dimensions for \code{new.reduction} output.
#' Corresponds to \code{n_latent} argument in the original API of SCANVI
#' @param n_hidden Number of nodes per hidden layer.
#' @param n_layers Number of hidden layers used for encoder and decoder NNs.
#' @param dropout_rate Dropout rate for neural networks.
#' @param dispersion One of the following:
#' \itemize{
#'  \item \code{gene}: dispersion parameter of NB is constant per gene across cells (default)
#'  \item \code{gene-batch}: dispersion can differ between different batches
#'  \item \code{gene-label}: dispersion can differ between different labels
#'  \item \code{gene-cell}: dispersion can differ for every gene in every cell
#' }
#' @param gene_likelihood One of the following:
#' \itemize{
#'  \item \code{zinb}: Zero-inflated negative binomial distribution (default)
#'  \item \code{nb}: Negative binomial distribution
#'  \item \code{poisson}: Poisson distribution
#' }
#' @param linear_classifier When switched to \code{TRUE}, uses a single linear layer for
#' classification instead of a multi-layer perceptron.
#' @param max_epochs Number of passes through the dataset for semisupervised training.
#' @param train_size  Size of training set in the range \code{[0.0, 1.0]}
#' @param batch_size Minibatch size to use during training.
#' @param verbose.scvi Verbosity level of scANVI. From quietest to talkiest:
#' CRITICAL, ERROR, WARNING, INFO (default), DEBUG, NOTSET
#' @param ... Additional arguments to be passed to \code{scvi.model.SCANVI},
#' \code{SCANVI.setup_anndata} or \code{SCANVI.train} (see \strong{Details} section)
#'
#' @return A list containing:
#' \itemize{
#'   \item a new DimReduc of name \code{new.reduction} (key set to
#'   \code{reduction.key}) consisting of the latent space of the model with
#'   \code{ndims.out} dimensions.
#' }
#' When called via \code{\link[Seurat]{IntegrateLayers}}, a Seurat object with
#' the new reduction and/or assay is returned
#'
#' @details
#' This wrappers calls three python functions through \pkg{reticulate}.
#' Find the \pkg{scVANVI}-specific arguments there:
#' \itemize{
#'   \item model initiation:
#'   \href{https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCANVI.html#scvi-model-scanvi}{scvi.model.SCANVI}, which relies on
#'   \href{https://docs.scvi-tools.org/en/stable/api/reference/scvi.module.SCANVAE.html#scvi-module-scanvae}{scvi.module.SCANVAE} which in turn relies on
#'   \href{https://docs.scvi-tools.org/en/stable/api/reference/scvi.module.VAE.html#scvi-module-vae}{scvi.module.VAE}
#'   \item anndata setup:
#'   \href{https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCANVI.html#scvi.model.SCANVI.setup_anndata}{SCANVI.setup_anndata}
#'   \item training:
#'   \href{https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCANVI.html#scvi.model.SCANVI.train}{SCANVI.train}
#' }
#'
#' @importFrom reticulate use_condaenv import r_to_py py_to_r
#' @importFrom Matrix t
#' @importFrom Seurat CreateDimReducObject
#' @importFrom SeuratObject JoinLayers GetAssayData
#'
#' @export
#' @note This function requires the
#' \href{https://scvi-tools.org/}{\pkg{scvi-tools}} package
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
#' obj <- IntegrateLayers(object = obj, method = scANVIIntegration,
#'                        features = Features(obj), conda_env = 'scvi-tools',
#'                        layers = 'counts', groups = obj[[]], groups.name = 'Method',
#'                        labels.name = 'CellType', labels.null = 'Unassigned')
#'
#' # To enable saving the model, add other 'nuisance' factors and increase number of threads used:
#' obj <- IntegrateLayers(object = obj, method = scANVIIntegration,
#'                        features = Features(obj), conda_env = 'scvi-tools',
#'                        layers = 'counts', groups = obj[[]], groups.name = "Method",
#'                        labels.name = "CellType", labels.null = "Unassigned",
#'                        categorical_covariate_keys = "Experiment",
#'                        continuous_covariate_keys = "percent.mito",
#'                        ncores = 8, model.save.dir = '~/Documents/scANVI.model')
#' }
#'
#' @references Kingma, D. P., Rezende, D. J., Mohamed, S. & Welling, M. Semi-
#' Supervised Learning with Deep Generative Models. Preprint at arXiv (2014).
#' \href{https://doi.org/10.48550/arXiv.1406.5298}{DOI}
#' @references Xu, C., Lopez, R., Mehlman, E., Regier, J., Jordan, M. I. & Yosef, N.
#' Probabilistic harmonization and annotation of single‐cell transcriptomics data
#' with deep generative models. Molecular Systems Biology 17, (2021).
#' \href{https://doi.org/10.15252/msb.20209620}{DOI}
#'
#' @seealso \code{\link[Seurat]{IntegrateLayers}}, \code{\link[Seurat]{writing-integration}}

scANVIIntegration <- function(
    object,
    groups = NULL,
    groups.name = NULL,
    labels.name = NULL,
    labels.null = NULL,
    features = NULL,
    layers = 'counts',
    scale.layer = 'scale.data',
    conda_env = NULL,
    new.reduction = 'integrated.scANVI',
    reduction.key = "scANVIlatent_",
    torch.intraop.threads = 4L,
    torch.interop.threads = NULL,
    # cuda.cores = NULL,
    model.save.dir = NULL,
    # scvi.model.SCANVI
    ndims.out = 10,
    n_hidden = 128L,
    n_layers = 1L,
    dropout_rate = 0.1,
    dispersion = c('gene', 'gene-batch', 'gene-label', 'gene-cell'),
    gene_likelihood = c("zinb", "nb", "poisson"),
    linear_classifier = FALSE,
    # SCANVI.train
    max_epochs = NULL,
    train_size = 0.9,
    batch_size = 128L,
    seed.use = 42L,
    verbose = TRUE,
    verbose.scvi = c("INFO", "NOTSET", "DEBUG", "WARNING", "ERROR", "CRITICAL"),
    ...){
  dispersion <- match.arg(arg = dispersion)
  gene_likelihood <- match.arg(arg = gene_likelihood)
  verbose.scvi <- toupper(verbose.scvi)
  verbose.scvi <- match.arg(arg = verbose.scvi)
  varargs <- list(...)

  conda_bin <- "auto"
  if (is.null(conda_env) || is.na(conda_env) || isFALSE(conda_env)) {
    if (! isValid(conda_status$current[["scanvi"]], do.check = TRUE)) {
      abort(message = paste("scANVI conda environment is not valid. Either",
                            "set", sQuote("conda_env"), "argument or create",
                            "the environment via the conda manager"))
    }
    message("Using conda from conda environment manager\n"[verbose], appendLF = FALSE)
    conda_env <- conda_status$current[["scanvi"]][["conda.env.path"]]$value
    conda_bin <- conda_status$current[["scanvi"]][["conda.bin"]]$value
  }

  use_condaenv(conda_env, conda = conda_bin, required = TRUE)
  sc <-  import('scanpy', convert = FALSE)
  torch <- import("torch", convert=FALSE)
  scvi <-  import('scvi', convert = FALSE)
  seed.use %iff% { scvi$settings$seed = as.integer(x = seed.use) }
  # cuda.cores %iff% { scvi$settings$num_threads = as.integer(x = cuda.cores) }
  scvi$settings$verbosity = args$scANVI$verbose[verbose.scvi]
  scipy <-  import('scipy', convert = FALSE)

  ncores.blas.old <- blas_get_num_procs()
  ncores.omp.old <- omp_get_num_procs()
  if ((torch.intraop.threads %iff% !is.na(as.integer(torch.intraop.threads))) %||% FALSE) {
    blas_set_num_threads(1L)
    omp_set_num_threads(1L)
    torch$set_num_threads(as.integer(torch.intraop.threads))
  }
  if ((torch.interop.threads %iff% !is.na(as.integer(torch.interop.threads))) %||% FALSE) {
    blas_set_num_threads(1L)
    omp_set_num_threads(1L)
    tryCatch({torch$set_num_interop_threads(as.integer(torch.interop.threads))},
             error = function(e) {
               warning("Number of inter-op threads was already set to ",
                       torch$get_num_interop_threads(),
                       "or parallel work has started. Cannot be changed, passing",
                       call. = FALSE, immediate. = TRUE)
             })
  }
  message(sprintf("%d intra-op and %d inter-op threads available to torch\n",
                  py_to_r(torch$get_num_threads()),
                  py_to_r(torch$get_num_interop_threads()))[verbose],
          appendLF = FALSE)

  layers <- Layers(object = object, search = layers %||% "counts")
  scale.layer <- scale.layer %||% "scale.data"
  groups <- groups %||% abort('A metadata table with cell type annotations is required')

  if (! inherits(x = groups, what = "data.frame")) {
    # groups is supposedly a vector, a matrix or a list
    groups <- as.data.frame(groups)
  }
  groups.name %||% {
    groups <- cbind(groups, CreateIntegrationGroups(object = object,
                                                    layers = layers,
                                                    scale.layer = scale.layer))
    groups.name <- colnames(groups)[ncol(groups)]
  }
  groups.name <- intersect(colnames(groups), groups.name)
  labels.name <- intersect(colnames(groups), labels.name)
  if(! length(groups.name)) {
    abort(message="'groups.name' not in 'groups' data frame")
  }
  if(! length(labels.name)) {
    abort(message="'labels.name' not in 'groups' data frame")
  }
  if(length(groups.name) > 1) {
    groups.name <- groups.name[1]
    warning(paste("more 'groups.name' that expected. Using the first one",
                  sQuote(x = groups.name)), call. = FALSE, immediate. = TRUE)
  }
  if(length(labels.name) > 1) {
    labels.name <- labels.name[1]
    warning(paste("more 'labels.name' that expected. Using the first one",
                  sQuote(x = labels.name)), call. = FALSE, immediate. = TRUE)
  }

  layer <- unique(sub("\\..*", "", layers %||% "counts"))
  if(length(layer) > 1) {
    abort(message="cannot find a consensus layer")
  }

  object <- JoinLayers(object = object, layers = layer)

  adata <- sc$AnnData(
    X   = scipy$sparse$csr_matrix(
      t( GetAssayData(object, layer = layer)[features ,] )
    ),
    obs = r_to_py(groups),
    var = r_to_py(data.frame(features = features))
  )
  # ValueError raised if unlabeled_category = NULL/None => use unlikely cell label)
  labels.null <- labels.null %||% "Gloubi-boulga"
  args.call <- c(list(adata = adata, labels_key = r_to_py(labels.name),
                      unlabeled_category = r_to_py(labels.null),
                      layer = r_to_py(NULL),
                      batch_key = r_to_py(groups.name)),
                 varargs[intersect(names(varargs), args$scANVI$setup_anndata)])
  do.call(scvi$model$SCANVI$setup_anndata, args.call)

  args.call <- c(list(adata = adata, n_hidden = r_to_py(as.integer(n_hidden)),
                      n_latent = r_to_py(as.integer(ndims.out)),
                      n_layers = r_to_py(as.integer(n_layers)),
                      dropout_rate = r_to_py(dropout_rate),
                      dispersion = r_to_py(dispersion),
                      gene_likelihood = r_to_py(gene_likelihood),
                      linear_classifier = r_to_py(linear_classifier)),
                 varargs[intersect(names(varargs), args$scANVI$scANVI)])
  model = do.call(scvi$model$SCANVI, args.call)

  max_epochs <- max_epochs %iff% as.integer(x = max_epochs)
  args.call <- c(list(max_epochs = r_to_py(max_epochs),
                      train_size = r_to_py(train_size),
                      batch_size = r_to_py(as.integer(batch_size))),
                 varargs[intersect(names(varargs), args$scANVI$train)])

  do.call(model$train, args.call)


  model.save.dir %iff% model$save(dir_path = path.expand(model.save.dir),
                                  overwrite = TRUE,
                                  save_anndata = F) # buggy, str conversion not working

  latent = model$get_latent_representation()
  latent <- as.matrix(latent)
  rownames(latent) <- py_to_r(adata$obs$index$values)
  colnames(latent) <- paste0(new.reduction, "_", 1:ncol(latent))
  suppressWarnings(latent.dr <- CreateDimReducObject(embeddings = latent, key = reduction.key))
  output.list <- list(latent.dr)
  names(output.list) <- new.reduction
  return(output.list)
}

attr(x = scANVIIntegration, which = 'Seurat.method') <- 'integration'
