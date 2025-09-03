#' @include kwargs.R
#'
NULL

#' Run trVAE on Seurat's \link[SeuratObject]{Assay5} object through \code{\link[Seurat]{IntegrateLayers}}
#'
#' @description
#' A wrapper to run \code{trVAE} on multi-layered Seurat V5 object.
#' Requires a conda environment with \pkg{scArches} and necessary dependencies
#'
#'
#' \strong{Recommendations}: use raw counts (except for \code{recon.loss = "mse"})
#' and all features (\code{features = Features(object), layers = "counts",
#' scale.layer = NULL}).
#'
#' @inheritParams integration-method
#' @param groups A \bold{named} data frame with grouping information. Can also
#' contain surgery groups to perform surgery integration.
#' @param surgery.name Column name from \code{groups} data frame that stores
#' surgery information. If \code{surgery.name = NULL}, a one shot integration is
#' performed
#' @param surgery.sort change the order in which surgery groups are integrated.
#' By default (\code{surgery.sort = TRUE}), surgery groups are ordered by name.
#' When \code{FALSE}, each group is integrated in the order of first occurrence
#' in the column \code{surgery.name}
#' @param conda_env Path to conda environment to run trVAE (should also
#' contain the scipy python module).  By default, uses the conda environment
#' registered for trVAE in the conda environment manager
#' @param torch.intraop.threads Number of intra-op threads available to torch
#' when training on CPU instead of GPU. Set via \code{torch.set_num_threads()}.
#' @param torch.interop.threads Number of intra-op threads available to torch
#' when training on CPU instead of GPU. Set via \code{torch.set_num_interop_threads()}.
#' Can only be changed once, on first call.
#' @param model.save.dir Path to a directory to save the model(s) to. Uses
#' \code{TRVAE.save()}. Does not save anndata. \code{model.save.dir = NULL}
#' (default) disables saving the model(s).
#' @param ndims.out Number of dimensions for \code{new.reduction} output.
#' Corresponds to \code{latent_dim} argument in the original API of TRVAE from
#' \pkg{scArches}
#' @param recon.loss Definition of Reconstruction-Loss-Method. One of 'mse', 'nb'
#' or 'zinb' (hence mean squared error, negative binomial and zero-inflated
#' negative binomial respectively). Recommended to set \code{layer = "data"} for
#' 'mse' (and \code{layer = "counts"} for (zi)nb)
#' @param hidden_layer_sizes Hidden layer sizes for encoder network
#' @param dr_rate Dropout rate applied to all layers. \code{dr_rate = 0} disables
#' dropout.
#' @param use_mmd Whether an additional MMD loss is to be calculated on the
#' latent dim. (see next argument)
#' @param mmd_on Choose on which layer MMD loss will be calculated on. One of 'z'
#' for latent dim or 'y' for the first decoder layer. Only applies when
#' \code{use_mmd = TRUE}
#' @param mmd_boundary On how many groups the MMD loss should be calculated on.
#' If \code{mmd_boundary = NULL} (default), MMD is calculated on all groups.
#' Only applies when \code{use_mmd = TRUE}
#' @param beta Scaling factor for MMD loss (1 by default). Only applies when
#' \code{use_mmd = TRUE}
#' @param use_bn Whether to apply a batch normalization to layers
#' @param use_ln Whether to apply a layer normalization to layers
#' @param n_epochs Maximum number of epochs to train the model
#' @param lr Learning rate for training
#' @param eps \code{torch.optim.Adam} eps parameter to improve numerical stability
#' (see \href{https://pytorch.org/docs/stable/generated/torch.optim.Adam.html}{here})
#' @param hide.py.warn Disables some uninformative warnings from torch
#' @param ... Additional arguments to be passed to
#' \code{scarches.models.TRVAE.train}, \code{TRVAE.load_query_data}
#' or \code{TRVAE.get_latent} (see \strong{Details} section)
#'
#' @return A list containing:
#' \itemize{
#'   \item \strong{Without surgery groups:} a new DimReduc of name
#'   \code{new.reduction} (key set to \code{reduction.key}) consisting of the
#'   latent space of the model with \code{ndims.out} dimensions.
#'   \item \strong{With surgery groups:} one new DimReduc per surgery groups of
#'   name \code{new.reduction_[surgery.group]} (key set to
#'   \code{reduction.key[surgery.group]}) consisting of the latent space of the
#'   corresponding models with \code{ndims.out} dimensions,  as well as a 'full'
#'   latent representation of name \code{new.reduction_[surgery1]_[surgery2]_...}
#'   and key set to \code{reduction.keyFull-}.
#' }
#' When called via \code{\link[Seurat]{IntegrateLayers}}, a Seurat object with
#' the new reduction and/or assay is returned
#'
#' @details
#' This wrappers calls three to four python functions through \pkg{reticulate}.
#' Find the \pkg{trVAE}-specific arguments there:
#' \itemize{
#'   \item model initiation:
#'   \href{https://docs.scarches.org/en/latest/api/models.html#scarches.models.TRVAE}{scarches.models.TRVAE}
#'   \item training:
#'   \href{https://docs.scarches.org/en/latest/api/models.html#scarches.models.TRVAE.train}{TRVAE.train}, which relies on
#'   \href{https://github.com/theislab/scarches/blob/51a0294ca987dabffb6d109178e0f69a90f9c24f/scarches/trainers/trvae/trainer.py#L14}{scarches.trainers.trvae.train.Trainer}
#'   \item post-training:
#'   \href{https://github.com/theislab/scarches/blob/51a0294ca987dabffb6d109178e0f69a90f9c24f/scarches/models/base/_base.py#L285}{scarches.models.base._base.CVAELatentsMixin.get_latent}
#'   \item surgery initiation:
#'   \href{https://github.com/theislab/scarches/blob/51a0294ca987dabffb6d109178e0f69a90f9c24f/scarches/models/base/_base.py#L200}{scarches.models.base._base.SurgeryMixin.load_query_data}
#' }
#'
#' @details  Note that \code{seed.use} is passed to \code{torch.manual_seed()}.
#' If it is not sufficient to achieve full reproducibility, set
#' \code{mean = TRUE} or \code{mean_var = TRUE}
#'
#' @importFrom reticulate use_condaenv import import_builtins r_to_py py_to_r
#' @importFrom RhpcBLASctl blas_get_num_procs blas_set_num_threads omp_get_num_procs omp_set_num_threads
#' @importFrom Matrix t
#' @importFrom Seurat CreateDimReducObject
#' @importFrom SeuratObject JoinLayers GetAssayData
#'
#' @export
#' @note This function requires the
#' \href{https://docs.scarches.org/}{\pkg{scArches}} package
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
#' obj <- IntegrateLayers(object = obj, method = trVAEIntegration,
#'                        features = Features(obj), scale.layer = NULL,
#'                        layers = 'counts', groups = obj[[]],
#'                        groups.name = 'Method')
#'
#' # To enable surgery and full reproducibility and change the recon loss method:
#' obj <- IntegrateLayers(object = obj, method = trVAEIntegration,
#'                        features = Features(obj), scale.layer = NULL,
#'                        layers = 'data', groups = obj[[]],
#'                        groups.name = 'Method', surgery.name = 'Experiment',
#'                        mean_var = TRUE, recon.loss = 'mse')
#' }
#'
#' @references Lotfollahi, M., Naghipourfar, M., Theis, F. J. & Wolf, F. A.
#' Conditional out-of-distribution generation for unpaired data using transfer
#' VAE. Bioinformatics 36, i610–i617 (2020).
#' \href{https://doi.org/10.1093/bioinformatics/btaa800}{DOI}
#' @references Lotfollahi, M., Naghipourfar, M., Luecken, M. D., Khajavi, M.,
#' Büttner, M., Wagenstetter, M., Avsec, Ž., Gayoso, A., Yosef, N., Interlandi,
#' M., Rybakov, S., Misharin, A. V. & Theis, F. J. Mapping single-cell data to
#' reference atlases by transfer learning. Nat Biotechnol 40, 121–130 (2021).
#' \href{https://doi.org/10.1038/s41587-021-01001-7}{DOI}
#'
#' @seealso \code{\link[Seurat]{IntegrateLayers}}, \code{\link[Seurat]{writing-integration}}

trVAEIntegration <- function(
    object,
    orig = NULL,
    groups = NULL,
    groups.name = NULL,
    surgery.name = NULL,
    surgery.sort = TRUE,
    features = NULL,
    layers = ifelse(recon.loss == "mse", 'data', 'counts'),
    scale.layer = 'scale.data',
    conda_env = NULL,
    new.reduction = "integrated.trVAE",
    reduction.key = "trVAElatent_",
    torch.intraop.threads = 4L,
    torch.interop.threads = NULL,
    model.save.dir = NULL,
    # sca.models.TRVAE params
    ndims.out = 10L, # 'latent_dim'
    recon.loss = c('nb', 'zinb', 'mse'),
    hidden_layer_sizes = c(256L, 64L),
    dr_rate = .05,
    use_mmd = TRUE,
    mmd_on = c('z', 'y'),
    mmd_boundary = NULL,
    beta = 1,
    use_bn = FALSE,
    use_ln = TRUE,
    # sca.models.TRVAE.train params
    n_epochs=400L,
    lr=1e-3,
    eps=.01,
    hide.py.warn = T,
    seed.use =42L,
    verbose = TRUE,
    ...) {
  recon.loss <- match.arg(recon.loss)
  expected.layer <- switch (recon.loss,
                            'nb' = 'counts', 'zinb' = 'counts', 'mse' = 'data')
  ndims.out <- as.integer(ndims.out %||% 10L)
  hidden_layer_sizes <- as.integer(hidden_layer_sizes %||% c(256L, 64L))
  mmd_on <- match.arg(mmd_on)
  use_mmd <- use_mmd %||% FALSE
  use_bn <- use_bn %||% FALSE
  use_ln <- use_ln %||% FALSE
  dr_rate <- as.numeric(dr_rate %||% .05)
  beta <- as.numeric(beta %||% 1)
  n_epochs <- as.integer(n_epochs %||% 400L)
  lr <- as.numeric(lr %||% 1e-3)
  eps <- as.numeric(eps %||% .01)
  seed.use <- seed.use %iff% as.integer(x = seed.use)
  varargs <- list(...)

  conda_bin <- "auto"
  if (is.null(conda_env) || is.na(conda_env) || isFALSE(conda_env)) {
    if (! isValid(conda_status$current[["trvae"]], do.check = TRUE)) {
      abort(message = paste("trVAE conda environment is not valid. Either",
                            "set", sQuote("conda_env"), "argument or create",
                            "the environment via the conda manager"))
    }
    message("Using conda from conda environment manager\n"[verbose], appendLF = FALSE)
    conda_env <- conda_status$current[["trvae"]][["conda.env.path"]]$value
    conda_bin <- conda_status$current[["trvae"]][["conda.bin"]]$value
  }

  use_condaenv(conda_env, conda = conda_bin, required = TRUE)
  if(hide.py.warn) {
    warnings <- import("warnings", convert=FALSE)
    builtins <- import_builtins(convert=FALSE)
    warnings$simplefilter(action='ignore', category=builtins$FutureWarning)
    warnings$simplefilter(action='ignore', category=builtins$UserWarning)

    logging <- import("logging", convert = FALSE)
    logging$basicConfig(level=logging$ERROR)
  }
  sc <- import("scanpy", convert=FALSE)
  torch <- import("torch", convert=FALSE)
  sca <- import("scarches", convert=FALSE)
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
  seed.use %iff% { torch$manual_seed(seed.use) ;
    message(sprintf("Set torch's manual seed to %d\n", seed.use)[verbose],
          appendLF = FALSE) }


  message(sprintf("Using %d features\n"[verbose], length(features)), appendLF = F)

  layers <- Layers(object = object, search = layers %||% expected.layer)
  scale.layer <- scale.layer %||% "scale.data"

  groups <- groups %||% CreateIntegrationGroups(object = object,
                                                layers = layers,
                                                scale.layer = scale.layer)
  groups.name <- groups.name %||% colnames(groups)[1]
  groups.name <- intersect(colnames(groups), groups.name)
  if (! length(x = groups.name)) {
    abort(message = "'groups.name' not in 'groups' data frame")
  }
  if (length(x = groups.name) > 1) {
    groups.name <- groups.name[1]
    warning(paste("more 'groups.name' that expected. Using the first one",
                  sQuote(x = groups.name)), call. = FALSE, immediate. = TRUE)
  }

  surgery.name %||% {
    surgery.name <- "surgery.order"
    groups[,"surgery.order"] <- 1
    message("trVAE model will be trained on all samples at once\n"[verbose], appendLF = FALSE)}
  surgery.name <- intersect(colnames(groups), surgery.name)
  if (! length(x = surgery.name)) {
    abort(message = "'surgery.name' not in 'groups' data frame")
  }
  if (length(x = surgery.name) > 1) {
    surgery.name <- surgery.name[1]
    warning(paste("more 'surgery.name' that expected. Using the first one",
                  sQuote(x = surgery.name)), call. = FALSE, immediate. = TRUE)
  }
  surgery.groups <- unique(as.vector(groups[, surgery.name, drop = TRUE]))
  surgery.sort <- surgery.sort %||% FALSE
  if (surgery.sort) {
    surgery.groups <- sort(surgery.groups)
  }
  message(sprintf("%d surgery groups\n", length(surgery.groups))[verbose],
          appendLF = FALSE)

  layer <- unique(sub("\\..*", "", layers %||% expected.layer))
  if(length(layer) > 1) {
    abort(message="cannot find a consensus layer")
  }
  if (layer != expected.layer) {
    warning("With recon.loss = ", sQuote(recon.loss), ", ",
            sQuote(expected.layer), " layer is expected. layers = ",
            sQuote(layer), " was input, which might produce invalid results.",
            call. = FALSE, immediate. = TRUE)
  }

  object <- JoinLayers(object = object, layers = layer)
  adata <- sc$AnnData(
    X   = scipy$sparse$csr_matrix(
      Matrix::t( GetAssayData(object, layer = layer)[features, ] )
    ),
    obs = reticulate::r_to_py(groups[, c(groups.name, surgery.name), drop = FALSE]),
    var = r_to_py(data.frame(features = features))
  )

  adata <- sca$dataset$remove_sparsity(adata)

  # Prevent 'RuntimeError: mat1 and mat2 must have the same dtype, but got Double and Float'
  # in torch/nn/modules/linear.py
  adata$X <- adata$X$astype('float32')

  adatas <- sapply(surgery.groups, function(grp) {
    adata[adata$obs[surgery.name] == grp]
  }, simplify = FALSE, USE.NAMES = TRUE)

  iter <- 0
  trVAE.models <- list()
  latents <- list()
  latents_var <- list()
  add_var <- varargs[["mean_var"]] %||% FALSE

  for (surgery.group in surgery.groups) {
    iter <- iter + 1
    print(iter)
    if (iter == 1) {
      args.call <- c(list(adata = adatas[[surgery.group]],
                          condition_key=groups.name,
                          conditions=adatas[[surgery.group]]$obs[groups.name]$unique()$tolist(),
                          hidden_layer_sizes=r_to_py(hidden_layer_sizes),
                          latent_dim=r_to_py(ndims.out),dr_rate=r_to_py(dr_rate),
                          use_mmd=r_to_py(use_mmd), mmd_on=r_to_py(mmd_on),
                          mmd_boundary=r_to_py(mmd_boundary),
                          recon_loss=r_to_py(recon.loss), beta=r_to_py(beta),
                          use_bn=r_to_py(use_bn), use_ln=r_to_py(use_ln)),
                     varargs[intersect(names(varargs), args$trVAE$TRVAE)])
      trvae <- do.call(sca$models$TRVAE, args.call)

      args.call <- c(list(n_epochs=r_to_py(n_epochs), lr=r_to_py(lr), eps=r_to_py(eps),
                          seed=r_to_py(as.integer(seed.use %||% 2020L))),
                     varargs[intersect(names(varargs), args$trVAE$train)])
      do.call(trvae$train, args.call)
    } else {
      args.call <- c(list(adata=adatas[[surgery.group]],
                          reference_model=trVAE.models[[iter - 1]]),
                     varargs[intersect(names(varargs), args$trVAE$load_query_data)])
      trvae <- do.call(sca$models$TRVAE$load_query_data, args.call)

      args.call <- c(list(n_epochs=r_to_py(n_epochs), lr=r_to_py(lr), ep=r_to_py(eps),
                          seed=r_to_py(as.integer(seed.use %||% 2020L))),
                     varargs[intersect(names(varargs), args$trVAE$train)])
      do.call(trvae$train, args.call)
    }
    trVAE.models[[surgery.group]] <- trvae
    args.call <- c(list(x=adatas[[surgery.group]]$X,
                        c=adatas[[surgery.group]]$obs[groups.name]),
                   varargs[intersect(names(varargs), args$trVAE$get_latent)])
    latent <- do.call(trvae$get_latent, args.call)
    if (add_var && is(object = latent, class2 = "python.builtin.tuple")) {
      latents_var[[surgery.group]] <- colMeans(sqrt(py_to_r(latent[1])))
      latent <- latent[0]
    }
    latent <- py_to_r(latent)
    colnames(latent) <- paste0(sub("_$", "", reduction.key),
                               surgery.group[length(surgery.groups) > 1],
                               "_", 1:ncol(latent))
    rownames(latent) <- py_to_r(adatas[[surgery.group]]$obs$index$values)
    latents[[surgery.group]] <- latent
    if (iter == length(surgery.groups) & iter > 1) {
      args.call <- c(list(x=adata$X, c=adata$obs[groups.name]),
                     varargs[intersect(names(varargs), args$trVAE$get_latent)])
      latent <- do.call(trvae$get_latent, args.call)
      if (add_var && is(object = latent, class2 = "python.builtin.tuple")) {
        latents_var[[paste(surgery.groups, collapse = "_")]] <- colMeans(sqrt(py_to_r(latent[1])))
        latent <- latent[0]
      }
      latent <- py_to_r(latent)
      # latent <- py_to_r(trvae$get_latent(adata$X, adata$obs[groups.name]))
      colnames(latent) <- paste0(sub("_$", "Full_", reduction.key), 1:ncol(latent))
      # colnames(latent) <- paste0("trVAE.latent.full-", 1:ncol(latent))
      rownames(latent) <- py_to_r(adata$obs$index$values)
      latents[[paste(surgery.groups, collapse = "_")]] <- latent
    }
  }

  model.save.dir %iff% {
    if (!dir.exists(model.save.dir)) {
      dir.create(model.save.dir, recursive = T)
    }
    sapply(names(trVAE.models), function(dir.name) {
      trVAE.models[[dir.name]]$save(
        dir_path = path.expand(file.path(model.save.dir, dir.name)),
        overwrite = TRUE, save_anndata = F
      )
    })
  }

  blas_set_num_threads(ncores.blas.old)
  omp_set_num_threads(ncores.omp.old)

  if (length(latents) > 1) {
    names(latents) <- paste(new.reduction, names(latents), sep="_")
  } else {
    names(latents) <- new.reduction
  }
  if (add_var && length(latents_var) == length(latents)) {
    names(latents_var) <- names(latents)
  }

  output.list <- sapply(names(latents), function(reduc.name) {
    key <- sub("_1$", "_", colnames(latents[[reduc.name]])[1])
    stdev <- latents_var[[reduc.name]] %||% numeric(0)

    CreateDimReducObject(embeddings = py_to_r(latents[[reduc.name]]),
                         key = key, stdev = stdev)
  }, simplify = FALSE, USE.NAMES = TRUE)

  return(output.list)
}

attr(x = trVAEIntegration, which = 'Seurat.method') <- 'integration'

