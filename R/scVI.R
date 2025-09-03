#' @include kwargs.R
#'
NULL

#' Run scVI on Seurat's \link[SeuratObject]{Assay5} object through \code{\link[Seurat]{IntegrateLayers}}
#'
#' @description
#' A wrapper to run \code{scVI} on multi-layered Seurat V5 object.
#' Requires a conda environment with \code{scvi-tools} and necessary dependencies
#'
#' Can be called via \code{SeuratIntegrate::scVIIntegration()} or
#' \code{scVIIntegration.fix()}
#'
#' \strong{Recommendations}: use raw counts and all features
#' (\code{features = Features(object), layers = "counts"})
#'
#' @inheritParams scANVIIntegration
#' @inherit scANVIIntegration return
#' @inheritParams integration-method
#' @param groups A \bold{named} data frame with grouping information.
#' @param conda_env Path to conda environment to run scVI (should also
#' contain the scipy python module).  By default, uses the conda environment
#' registered for scVI in the conda environment manager
#' @param model.save.dir Path to a directory to save the model to. Uses
#' \code{SCVI.save()}. Does not save anndata. Note that neither the trainer
#' optimizer state nor the trainer history are saved.
#' \code{model.save.dir = NULL} (default) disables saving the model.
#' @param ndims.out Number of dimensions for \code{new.reduction} output.
#' Corresponds to \code{n_latent} argument in the original API of SCVI
#' @param latent_distribution One of the following:
#' \itemize{
#'  \item \code{normal}: Normal distribution (default)
#'  \item \code{ln}: Logistic normal distribution (Normal(0, I) transformed by softmax)
#' }
#' @param verbose.scvi Verbosity level of scVI. From quietest to talkiest:
#' CRITICAL, ERROR, WARNING, INFO (default), DEBUG, NOTSET
#' @param ... For \code{scVIIntegration()}, additional arguments to be passed to
#' \code{scvi.model.SCVI}, \code{SCVI.setup_anndata} or \code{SCVI.train} (see
#' \strong{Details} section). For \code{scVIIntegration.fix()}, all of the above
#'
#' @details
#' This wrappers calls three python functions through \pkg{reticulate}.
#' Find the \pkg{scVI}-specific arguments there:
#' \itemize{
#'   \item model initiation:
#'   \href{https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html#scvi.model.SCVI}{scvi.model.SCVI}, which relies on
#'   \href{https://docs.scvi-tools.org/en/stable/api/reference/scvi.module.VAE.html#scvi.module.VAE}{scvi.module.VAE}
#'   \item anndata setup:
#'   \href{https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html#scvi.model.SCVI.setup_anndata}{SCVI.setup_anndata}
#'   \item training:
#'   \href{https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html#scvi.model.SCVI.train}{SCVI.train}
#' }
#'
#' @importFrom reticulate use_condaenv import r_to_py py_to_r
#' @importFrom Matrix t
#' @importFrom Seurat CreateDimReducObject
#' @importFrom SeuratObject JoinLayers GetAssayData
#'
#' @export
#'@note This function requires the
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
#' obj <- IntegrateLayers(object = obj, method = scVIIntegration,
#'                        features = Features(obj), conda_env = 'scvi-tools',
#'                        layers = 'counts', groups = obj[[]], groups.name = 'Method')
#'
#' # To enable cell label-guided correction, save the model, add other
#' # 'nuisance' factors and increase number of threads used:
#' obj <- IntegrateLayers(object = obj, method = scVIIntegration,
#'                        features = Features(obj), conda_env = 'scvi-tools',
#'                        layers = 'counts', groups = obj[[]], groups.name = "Method",
#'                        labels.name = "CellType",
#'                        categorical_covariate_keys = list("Experiment"),
#'                        continuous_covariate_keys = list("percent.mito"),
#'                        ncores = 8, model.save.dir = '~/Documents/scVI.model')
#' }
#'
#' @references Lopez, R., Regier, J., Cole, M. B., Jordan, M. I. & Yosef, N.
#' Deep generative modeling for single-cell transcriptomics. Nat Methods 15,
#' 1053–1058 (2018). \href{https://doi.org/10.1038/s41592-018-0229-2}{DOI}
#'
#' @seealso \code{\link[Seurat]{IntegrateLayers}}, \code{\link[Seurat]{writing-integration}}

scVIIntegration <- function(
    object,
    groups = NULL,
    groups.name = NULL,
    labels.name = NULL,
    features = NULL,
    layers = 'counts',
    scale.layer = 'scale.data',
    conda_env = NULL,
    new.reduction = 'integrated.scVI',
    reduction.key = "scVIlatent_",
    torch.intraop.threads = 4L,
    torch.interop.threads = NULL,
    model.save.dir = NULL,
    # scvi.model.SCVI
    ndims.out = 10,
    n_hidden = 128L,
    n_layers = 1L,
    dropout_rate = 0.1,
    dispersion = c('gene', 'gene-batch', 'gene-label', 'gene-cell'),
    gene_likelihood = c("zinb", "nb", "poisson"),
    latent_distribution = c('normal', 'ln'),
    # SCVI.train
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
    if (! isValid(conda_status$current[["scvi"]], do.check = TRUE)) {
      abort(message = paste("scVI conda environment is not valid. Either",
                            "set", sQuote("conda_env"), "argument or create",
                            "the environment via the conda manager"))
    }
    message("Using conda from conda environment manager\n"[verbose], appendLF = FALSE)
    conda_env <- conda_status$current[["scvi"]][["conda.env.path"]]$value
    conda_bin <- conda_status$current[["scvi"]][["conda.bin"]]$value
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
  groups <- groups %||% CreateIntegrationGroups(object = object,
                                                layers = layers,
                                                scale.layer = scale.layer)
  if (! inherits(x = groups, what = "data.frame")) {
    # groups is supposedly a vector, a matrix or a list
    groups <- as.data.frame(groups)
  }
  groups.name <- intersect(colnames(groups), groups.name %||% colnames(groups)[1])
  if(! length(groups.name)) {
    abort(message="'groups.name' not in 'groups' data frame")
  }
  if(length(groups.name) > 1) {
    groups.name <- groups.name[1]
    warning(paste("more 'groups.name' that expected. Using the first one",
                  sQuote(x = groups.name)), call. = FALSE, immediate. = TRUE)
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

  args.call <- c(list(adata = adata, labels_key = r_to_py(labels.name),
                      layer = r_to_py(NULL),
                      batch_key = r_to_py(groups.name)),
                 varargs[intersect(names(varargs), args$scVI$setup_anndata)])
  do.call(scvi$model$SCVI$setup_anndata, args.call)

  args.call <- c(list(adata = adata, n_hidden = r_to_py(as.integer(n_hidden)),
                      n_latent = r_to_py(as.integer(ndims.out)),
                      n_layers = r_to_py(as.integer(n_layers)),
                      dropout_rate = r_to_py(dropout_rate),
                      dispersion = r_to_py(dispersion),
                      gene_likelihood = r_to_py(gene_likelihood),
                      latent_distribution = r_to_py(latent_distribution)),
                 varargs[intersect(names(varargs), args$scVI$scVI)])
  model = do.call(scvi$model$SCVI, args.call)

  max_epochs <- max_epochs %iff% as.integer(x = max_epochs)
  args.call <- c(list(max_epochs = r_to_py(max_epochs),
                      train_size = r_to_py(train_size),
                      batch_size = r_to_py(as.integer(batch_size))),
                 varargs[intersect(names(varargs), args$scVI$train)])

  do.call(model$train, args.call)

  model.save.dir %iff% model$save(dir_path = path.expand(model.save.dir),
                                  overwrite = TRUE,
                                  save_anndata = F)

  latent = model$get_latent_representation()
  latent <- as.matrix(latent)
  rownames(latent) <- py_to_r(adata$obs$index$values)
  colnames(latent) <- paste0(new.reduction, "_", 1:ncol(latent))
  suppressWarnings(latent.dr <- CreateDimReducObject(embeddings = latent, key = reduction.key))
  output.list <- list(latent.dr)
  names(output.list) <- new.reduction
  return(output.list)
}

attr(x = scVIIntegration, which = 'Seurat.method') <- 'integration'

#' @rdname scVIIntegration
#' @export
scVIIntegration.fix <- function(...) {
  scVIIntegration(...)
}
attr(x = scVIIntegration.fix, which = 'Seurat.method') <- 'integration'
