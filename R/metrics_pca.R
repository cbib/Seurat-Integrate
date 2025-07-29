#' Score a corrected or uncorrected PCA to estimate batch mixing
#'
#' @description
#' A batch-wise kernel density estimate is computed for each dimension of the PCA.
#' Then, the global area shared by all kernels is calculated and divided by a
#' reference area (see \strong{Arguments} and \strong{Details}) as a proxy to
#' batch mixing. Finally, the aforementioned proportions are weighted by each
#' dimension's contribution to variance.
#'
#' @param object A Seurat object
#' @param batch.var The name of the batch variable (must be in the object metadata)
#' @param reduction The name of the reduction to score
#' @param dims The dimensions to consider. All dimensions are used by default
#' @param use.union Whether to use a union of densities (constructs a fake
#' density with pmax). The alternative is to use all density kernels and sum them
#' @param bw.join A function to pick a joint bandwidth. The default is to use the
#' \code{mean()}. You can also chose \code{median} for instance. If you want to
#' use a distinct bandwidth for each batch, just set \code{bw.join = c} (not recommended)
#' @param bw.kernel Bandwidth selector to use for Gaussian kernel density estimation.
#' One of 'nrd0', 'nrd', 'ucv','bcv' or 'sj' (see \link[stats:bw.nrd]{bandwidths}).
#' nrd0 is the default
#' @param weight.by one of 'var' (default) or 'stdev' (standing for variance and
#' standard deviation respectively). Use the variance or the standard deviation
#' explained by the principal components to weight the each PC's score.
#' @param assay assay to use. Passed to Seurat to automatically construct the
#' \code{batch.var} when not provided. Useless otherwise
#' @param layer layer to use. Passed to Seurat to automatically construct the
#' \code{batch.var} when not provided. Useless otherwise
#'
#' @return \code{ScoreDensityPC}: A single float corresponding to the score of
#' the given reduction
#'
#' \code{AddScoreDensityPC}: the updated Seurat \code{object} with the density
#' PCA score set for the integration.
#'
#' @importFrom stats bw.nrd0 bw.nrd bw.ucv bw.bcv bw.SJ density
#' @importFrom dplyr `%>%` group_by summarise across all_of ungroup mutate pick
#'
#' @export
#' @details
#'  The score is computed as follow :
#' \deqn{\sum_{i=1}^{p} \left ( \cfrac{ \bigcup_{b=1}^n D_{bi} }{ Ref_i } * V_i \right )}
#'
#' For a PCA with p dimensions and n batchs. \eqn{D_{bi}} is the vector of density
#' of batch b on the dimension i. \eqn{V_i} is the proportion of variance
#' explained by the \eqn{PC_i}.
#'
#' with \eqn{Ref = \bigcap_{i=1}^n D_{bi}} or \eqn{Ref = \sum_{i=1}^n D_{bi}}
#' when \code{use.union} is \code{TRUE} and \code{FALSE} respectively, for a
#' given dimension \eqn{i}.
#'
#' \code{use.union = FALSE} tends to result in lower score since the
#' reference is the sum of kernel densities instead of maximum
#'
#' @note This score is an adaptation of the principal component regression (PCR)
#' score from Luecken M.D. \emph{et al.}, 2022.
#'
#' @examples
#' \dontrun{
#' obj <- SeuratData::LoadData("pbmcsca")
#' obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj)
#'
#' score.union <- ScoreDensityPC(obj, "Method", "pca", dim = 1:30)
#' score.sum <- ScoreDensityPC(obj, "Method", "pca", dim = 1:30, use.union = FALSE)
#'
#' score.union   # ~ 0.1511
#' score.sum     # ~ 0.0319
#' }
#'
#' @references Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A.,
#' Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M.,
#' Colomé-Tatché, M. & Theis, F. J. Benchmarking atlas-level data integration in
#' single-cell genomics. Nat Methods 19, 41–50 (2021).
#' \href{https://doi.org/10.1038/s41592-021-01336-8}{DOI}
#'
#' @seealso \code{\link{ScoreRegressPC}} to compute the PCR score
#' @rdname score-densityPC

ScoreDensityPC <- function(object, batch.var=NULL, reduction = "pca", dims=NULL,
                           use.union = TRUE, bw.join = mean,
                           bw.kernel = c('nrd0', 'nrd', 'ucv','bcv', 'sj'),
                           weight.by = c("var", "stdev"),
                           assay = NULL, layer = NULL) {
  bw.kernel <- tolower(bw.kernel)
  bw.kernel <- match.arg(bw.kernel)
  bw.kernel <- switch (bw.kernel, 'nrd0' = bw.nrd0, 'nrd' = bw.nrd,
                       'ucv' = bw.ucv, 'bcv' = bw.bcv, 'SJ' = bw.SJ)
  bw.join <- bw.join %||% .Primitive("c")
  if (! is.function(bw.join)) {
    abort(message = sprintf("%s is not a function", sQuote("bw.join")))
  }

  assay <- assay %||% DefaultAssay(object)
  DefaultAssay(object) <- assay
  .prep_MetaDataBatch(object = object, batch.var = batch.var, assay = assay,
                      layer = layer)
  .prep_MetadataPC(object = object, df.mtdt = df.mtdt, idcol = idcol,
                   batch.var = batch.var, reduction = reduction, dims = dims,
                   weight.by = weight.by)

  # compute joint bandwidth
  dimred_cols <- colnames(dimred)
  joint.bws <- df.score %>% group_by(pick({{ batch.var }})) %>%
    summarise(across({{ dimred_cols }}, bw.kernel)) %>%
    ungroup() %>% mutate(across({{ dimred_cols }}, bw.join)) %>%
    as.data.frame()

  # compute densities per batch and collect `density()$y`
  dens.dfs <- sapply(colnames(dimred), function(col.nm) {
    dens.mm <- range(df.score[, col.nm, drop = TRUE], na.rm = TRUE)
    sapply(unique(df.score[,batch.var]), function(batch.nm) {
      density(df.score[df.score[,batch.var] == batch.nm, col.nm],
              bw = joint.bws[batch.nm, col.nm],
              from = dens.mm[1], to = dens.mm[2])$y
    }, simplify = "array") %>% as.data.frame()
  }, simplify = F)

  # compute proportion of common (~intersect) area
  prop.area <- sapply(dens.dfs, function(dens) {
    inter <- do.call(pmin, dens)
    union <- do.call(pmax, dens)
    if (! use.union) {
      # ue sum instead
      union <- rowSums(dens)
    }
    sum(inter) / sum(union)
  }, simplify = "array")

  return(sum(prop.area * proportions(dimvar)))
}

#' @param integration name of the integration to score
#' @export
#' @rdname score-densityPC
AddScoreDensityPC <- function(object, integration,
                              batch.var=NULL, reduction = "pca", dims=NULL,
                              use.union = TRUE, bw.join = mean,
                              bw.kernel = c('nrd0', 'nrd', 'ucv','bcv', 'sj'),
                              weight.by = c("var", "stdev"),
                              assay = NULL, layer = NULL) {
  scores <- ScoreDensityPC(object, batch.var = batch.var, reduction = reduction,
                           dims = dims, use.union = use.union, bw.join = bw.join,
                           bw.kernel = bw.kernel, weight.by = weight.by,
                           assay = assay, layer = layer)

  object <- check_misc(object)
  object <- SetMiscScore(object, integration = integration,
                         score.name = "PCA.density", score.value = scores)
  return(object)
}

#' Score a corrected or uncorrected PCA to estimate batch mixing
#'
#' @description
#' Linearly regresses principal components with batch variable as a proxy to
#' estimate batch mixing. The resulting R2 are then weighted by each dimension's
#' contribution to variance.
#'
#' @inheritParams ScoreDensityPC
#' @param adj.r2 Whether to use the adjusted R2 instead of the raw R2
#'
#' @return \code{ScoreRegressPC}: A single float corresponding to the score of
#' the given reduction
#'
#' \code{AddScoreRegressPC}: the updated Seurat \code{object} with the regression
#' PCA score set for the integration.
#'
#' @importFrom stats lm summary.lm as.formula
#' @importFrom dplyr %>% group_by summarise across all_of ungroup mutate pick
#'
#' @export
#' @details The linear regression is
#' \deqn{PC_i = Batch}
#'
#'
#' The score is computed as follow :
#' \deqn{\sum_{i=1}^{p} \left ( R^2_i * V_i \right )}
#'
#' For a PCA with p dimensions, \eqn{PC_i} is the principal component i,
#' \eqn{R^2_i} is the R squared coefficient of the linear regression for the
#' dimension i. \eqn{V_i} is the proportion of variance explained by the
#' \eqn{PC_i}.
#'
#' @examples
#' \dontrun{
#' obj <- SeuratData::LoadData("pbmcsca")
#' obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj)
#'
#' score.r2 <- ScoreRegressPC(obj, "Method", "pca", dim = 1:30)
#' score.adj.r2 <- ScoreRegressPC(obj, "Method", "pca", dim = 1:30, adj.r2 = TRUE)
#'
#' score.r2    # ~ 0.1147
#' score.adj.r2     # ~ 0.1145
#' }
#'
#' @inherit ScoreDensityPC note
#' @inherit ScoreDensityPC references
#'
#' @seealso \code{\link{ScoreDensityPC}} for an alternative and
#' \code{\link{ScoreRegressPC.CellCycle}} to regresses PCs by cell cycle scores.
#' @rdname score-regressPC

ScoreRegressPC <- function(object, batch.var=NULL, reduction = "pca", dims=NULL,
                           adj.r2 = FALSE,  weight.by = c("var", "stdev"),
                           assay = NULL, layer = NULL) {
  adj.r2 <- adj.r2 %||% FALSE
  assay <- assay %||% DefaultAssay(object)
  DefaultAssay(object) <- assay
  .prep_MetaDataBatch(object = object, batch.var = batch.var, assay = assay,
                      layer = layer)
  .prep_MetadataPC(object = object, df.mtdt = df.mtdt, idcol = idcol,
                   batch.var = batch.var, reduction = reduction, dims = dims,
                   weight.by = weight.by)

  # construct formula(s)
  formulas <- lapply(paste(colnames(dimred), batch.var, sep = " ~ "), as.formula)

  # compute linear regression per batch and collect (adj.)r.squared
  regs <- lapply(formulas, FUN = lm, data = df.score)
  r2get <- paste0("adj."[adj.r2], "r.squared")
  r2 <- sapply(lapply(regs, summary.lm), getElement, name = r2get, simplify = "numeric")

  return(sum(r2 * proportions(dimvar)))
}

#' @param integration name of the integration to score
#' @export
#' @rdname score-regressPC
AddScoreRegressPC <- function(object, integration,
                              batch.var=NULL, reduction = "pca", dims=NULL,
                              adj.r2 = FALSE,  weight.by = c("var", "stdev"),
                              assay = NULL, layer = NULL) {
  scores <- ScoreRegressPC(object, batch.var = batch.var, reduction = reduction,
                           dims = dims, adj.r2 = adj.r2, weight.by = weight.by,
                           assay = assay, layer = layer)

  object <- check_misc(object)
  object <- SetMiscScore(object, integration = integration,
                         score.name = "PCA.regression", score.value = scores)
  return(object)
}

#' Score a corrected or uncorrected PCA to estimate the contribution of S and G2M
#' scores to variance
#'
#' @description
#' Linearly regresses S and G2M scores to predict principal components. The
#' resulting R2 are then weighted by each dimension's contribution to variance.
#' Cell cycles scores and PCAs are computed for each batch independently.
#'
#' @inheritParams ScoreRegressPC
#' @inheritParams Seurat::RunPCA
#' @inheritParams Seurat::CellCycleScoring
#' @param what the slot from Seurat to score. Can be a layer or a reduction.
#' @param dims.use The dimensions from \code{what} to consider. All dimensions
#' are used by default
#' @param s.var The name of the S phase score variable (must be in the object
#' metadata)
#' @param g2m.var The name of the G2M phase score variable (must be in the
#' object metadata)
#' @param compute.cc whether to (re-)compute the cell cycle scores. Should be
#' \code{TRUE} (default), unless you have run
#' \code{\link{CellCycleScoringPerBatch}} beforehand because cell cycles scores
#' are expected to be computed per batch
#'
#' @return \code{ScoreRegressPC.CellCycle}: A 2-columns data frame with the
#' batch variable in the first one and the corresponding score in the second
#' one. It has as many rows as batches.
#'
#' \code{AddScoreRegressPC.CellCycle}: the updated Seurat \code{object} with the
#' cell cycle conservation score set for the integration.
#'
#' @importFrom stats lm summary.lm
#' @importFrom SeuratObject Layers DefaultAssay DefaultAssay<- GetAssayData Reductions Embeddings
#' @importFrom Seurat CellCycleScoring RunPCA
#' @importFrom dplyr group_by group_modify left_join summarize pick
#' @importFrom tidyr pivot_longer
#' @importFrom broom glance
#' @importFrom rlang data_sym
#'
#'
#' @export
#' @details The linear regression is
#' \deqn{PC_i = S_{score} + G2M_{score}}
#'
#' The score is computed as follow :
#' \deqn{\sum_{i=1}^{p} \left ( R^2_i * V_i \right )}
#'
#' For a PCA with p dimensions, \eqn{PC_i} is the principal component i,
#' \eqn{R^2_i} is the R squared coefficient of the linear regression for the
#' dimension i. \eqn{V_i} is the proportion of variance explained by the
#' \eqn{PC_i}.
#'
#' @examples
#' \dontrun{
#' obj <- SeuratData::LoadData("pbmcsca")
#' obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj)
#'
#' score.cc.r2 <- ScoreRegressPC.CellCycle(obj, "Method", "pca", dim.use = 1:30)
#' score.cc.adj.r2 <- ScoreRegressPC.CellCycle(obj, "Method", "pca", dim.use = 1:30, adj.r2 = TRUE)
#'
#' score.cc.r2        # ~ 0.0249
#' score.cc.adj.r2    # ~ 0.0249
#' }
#'
#' @inherit ScoreDensityPC note
#' @inherit ScoreDensityPC references
#'
#' @seealso \code{\link{CellCycleScoringPerBatch}} to compute cc scores per
#' batch. \code{\link{ScoreRegressPC}} to regresses PCs by batch.
#' @rdname score-cc

ScoreRegressPC.CellCycle <- function(object, batch.var = NULL,
                                     what = NULL,
                                     dims.use = NULL, npcs = 50L,
                                     s.var = "S.Score", g2m.var = "G2M.Score",
                                     compute.cc = TRUE,
                                     s.features = NULL, g2m.features = NULL,
                                     assay = NULL, weight.by = c("var", "stdev"),
                                     adj.r2 = FALSE,
                                     approx = FALSE) {
  assay <- assay %||% DefaultAssay(object)
  weight.by <- tolower(weight.by)
  weight.by <- match.arg(weight.by)
  what <- what %||% 'scale.data'
  DefaultAssay(object) <- assay

  if (what %in% Layers(object)) {
    ref <- GetAssayData(object, assay = assay, layer = what)
  } else if (what %in% Reductions(object)) {
    ref <- t(Embeddings(object, what))
  } else {
    abort(sprintf('%s not found in the Seurat object\'s layers and reductions',
                  sQuote(what)))
  }
  dims.use <- dims.use %||% 1:nrow(ref)
  if (max(dims.use) > nrow(ref)) {
    l <- length(dims.use)
    dims.use <- intersect(dims.use, 1:nrow(ref))
    if (length(dims.use) == 0) {
      abort("All provided dimensions are out of range")
    }
    if (length(dims.use) == 1) {
      abort("A single dimension is not enough to compute a PCA.")
    }
    warning(sprintf('dropping %d out of range dimensions (%d retained)',
                    l - length(dims.use), length(dims.use)),
            call. = FALSE, immediate. = TRUE)
  }
  if (npcs > length(dims.use)) {
    npcs <- length(dims.use) - 1
    warning(sprintf('npcs to compute reduced to %d (number of dims - 1)', npcs),
            call. = FALSE, immediate. = TRUE)
  }

  .prep_MetaDataBatch(object = object, batch.var = batch.var,
                      assay = assay)

  cols.cc <- c('S.Score', 'G2M.Score', 'Phase')
  r2get <- paste0("adj."[adj.r2], "r.squared")
  df.scores.all <- data.frame(
    setNames(list(character(0), character(0), numeric(0), numeric(0)),
             c(batch.var, 'dimred_col', r2get, 'var')))
  for (batch in unique(df.mtdt[, batch.var])) {
    cells <- df.mtdt[df.mtdt[,batch.var] == batch, idcol, drop = TRUE]
    sub.object <- subset(object, cells = cells)
    DefaultAssay(sub.object) <- assay
    if (inherits(sub.object[[assay]], "StdAssay")) {
      sub.object <- JoinLayers(sub.object)
    }
    if (compute.cc) {
      sub.object <- CellCycleScoring(sub.object, s.features = s.features,
                                     g2m.features = g2m.features, ctrl = NULL,
                                     set.ident = FALSE)
    }
    sub.object[['pca.batch']] <- RunPCA(ref[dims.use, cells, drop=F],
                                        assay = assay, npcs = npcs,
                                        verbose = FALSE,
                                        reduction.key = "boulgiboulga_",
                                        approx = approx)

    s.var <- s.var %||% "S.Score"
    g2m.var <- g2m.var %||% "G2M.Score"
    sub.df.mtdt <- sub.object[[]][cells, , drop = FALSE] %>% rownames_to_column(idcol)
    .prep_MetadataPC(object = sub.object, df.mtdt = sub.df.mtdt, idcol = idcol,
                     batch.var = batch.var, s.var = s.var, g2m.var = g2m.var,
                     reduction = 'pca.batch', dims = NULL)

    covar <- paste(c(s.var, g2m.var), collapse = " + ")
    formula_ <- as.formula(paste('dimred_val', covar, sep = " ~ "))

    df.score <- df.score %>% pivot_longer(colnames(dimred),
                                          names_to = "dimred_col",
                                          values_to = "dimred_val") %>%
      group_by(pick({{ batch.var }}, "dimred_col")) %>%
      group_modify(~ glance(lm(formula_, data = .x))) %>%
      left_join(data.frame(list('dimred_col' = colnames(sub.object[['pca.batch']]),
                                'var' = dimvar)), by = "dimred_col")

    df.scores.all <- rbind(df.scores.all,
                           df.score[,c(batch.var, 'dimred_col', r2get, 'var')])
  }
  return(df.scores.all %>% group_by(pick({{ batch.var }})) %>%
           summarize(score = sum(!!data_sym(r2get) * proportions(var))))
}

#' @param integration name of the integration to score
#' @export
#' @rdname score-cc
AddScoreRegressPC.CellCycle <- function(object, integration,
                                        batch.var = NULL,
                                        what = NULL,
                                        dims.use = NULL, npcs = 50L,
                                        s.var = "S.Score", g2m.var = "G2M.Score",
                                        compute.cc = TRUE,
                                        s.features = NULL, g2m.features = NULL,
                                        assay = NULL,
                                        weight.by = c("var", "stdev"),
                                        adj.r2 = FALSE,
                                        approx = FALSE) {
  scores <- ScoreRegressPC.CellCycle(
    object, batch.var = batch.var, what = what, dims.use = dims.use, npcs = npcs,
    s.var = s.var, g2m.var = g2m.var, compute.cc = compute.cc,
    s.features = s.features, g2m.features = g2m.features, assay = assay,
    weight.by = weight.by, adj.r2 = adj.r2, approx = approx)

  object <- check_misc(object)
  object <- SetMiscScore(object, integration = integration,
                         score.name = "cell.cycle.conservation",
                         score.value = scores,
                         class = "list")
  return(object)
}

#' @importFrom SeuratObject Reductions Embeddings DefaultAssay Layers
#' @importFrom Seurat RunPCA
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @keywords internal
#' @noRd
.prep_MetadataPC <- function(object, df.mtdt, idcol, batch.var = NULL,
                             reduction = "pca", dims = NULL,
                             weight.by = c("var", "stdev"),
                             s.var = NULL, g2m.var = NULL) {
  reduction <- reduction %||% "pca"
  if (! reduction %in% Reductions(object)) {
    abort(message = sprintf("%s reduction not in object", sQuote(reduction)))
  }
  dimred <- Embeddings(object = object, reduction = reduction)
  dims <- unique(dims %||% 1:ncol(dimred))
  if (! all(dims %in% 1:ncol(dimred))) {
    abort(message = "some dims are out of range")
  }
  dimvar <- Reductions(object = object, slot = reduction)@stdev
  if (is.null(dimvar) || length(dimvar) == 0) {
    dimred <- suppressWarnings(RunPCA(t(dimred), npcs = ncol(dimred),
                                      approx = FALSE, verbose = FALSE))
    dimvar <- dimred@stdev
    dimred <- dimred@cell.embeddings
  }
  dimred <- dimred[, dims, drop=FALSE]
  dimvar <- dimvar[dims]
  if (grepl('^var', tolower(weight.by[1]))) {
    dimvar <- dimvar^2
  }

  if ((s.var %iff% ! all(s.var %in% colnames(df.mtdt))) %||% FALSE) {
    rlang::abort(
      message = sprintf("%s with S phase score not in colnames of metadata",
                        sQuote(s.var)))
  }
  if ((g2m.var %iff% ! all(g2m.var %in% colnames(df.mtdt))) %||% FALSE) {
    abort(
      message = sprintf("%s with G2M phase score not in colnames of metadata",
                        sQuote(g2m.var)))
  }

  df.mtdt[, batch.var] <- as.factor(df.mtdt[, batch.var])
  df.mtdt <- df.mtdt %>% left_join(
    as.data.frame(dimred) %>% rownames_to_column(var = idcol),
    by = idcol
  )
  list2env(list(batch.var = batch.var, dimred = dimred, dimvar = dimvar,
              dims = dims, df.score = df.mtdt), envir = parent.frame())
}

#' Score cell cycle phases per batch
#'
#' @description
#' Assign cell cycle scores to cells. Scores are computed for each batch
#' independantly.
#'
#' @inheritParams ScoreRegressPC
#' @inheritParams Seurat::CellCycleScoring
#' @param ... Arguments to be passed to \code{\link[Seurat]{CellCycleScoring}},
#' then \code{\link[Seurat]{AddModuleScore}} (with the exception of
#' \code{set.ident} which is always \code{FALSE})
#'
#' @inherit Seurat::CellCycleScoring return
#' @importFrom stats lm summary.lm
#' @importFrom SeuratObject DefaultAssay DefaultAssay<- JoinLayers AddMetaData
#' @importFrom Seurat CellCycleScoring
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- SeuratData::LoadData("pbmcsca")
#' obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- CellCycleScoringPerBatch(obj, batch.var = 'Method',
#'                                 s.features = cc.genes.updated.2019$s.genes,
#'                                 g2m.features = cc.genes.updated.2019$g2m.genes)
#'
#' head(obj[[]])
#' }
#' @seealso \code{\link[Seurat]{CellCycleScoring}} to compute cc scores on the whole
#' dataset.

CellCycleScoringPerBatch <- function(object, batch.var = NULL,
                                     s.features,
                                     g2m.features,
                                     ctrl = NULL,
                                     assay = NULL,
                                     layer = NULL,
                                     ...) {
  assay <- assay %||% DefaultAssay(object)
  assay.old <- DefaultAssay(object)
  DefaultAssay(object) <- assay
  idcol <- "cellbarcodeid"
  #list2env(batch.var, df.mtdt)
  .prep_MetaDataBatch(object = object, batch.var = batch.var,
                      assay = assay, layer = layer, idcol = idcol)

  cols.cc <- c('S.Score', 'G2M.Score', 'Phase')
  df.cc <- data.frame(setNames(list(numeric(0), numeric(0), character(0)), cols.cc))
  for (batch in unique(df.mtdt[, batch.var])) {
    cells <- df.mtdt[df.mtdt[, batch.var] == batch, idcol]
    sub.object <- subset(object, cells = cells)
    DefaultAssay(sub.object) <- assay
    if (inherits(sub.object[[assay]], "StdAssay")) {
      sub.object <- JoinLayers(sub.object)
    }
    sub.object <- CellCycleScoring(sub.object, s.features = s.features,
                                   g2m.features = g2m.features, ctrl = ctrl,
                                   set.ident = FALSE, ...)
    df.cc <- rbind(df.cc, sub.object[[]][, cols.cc])
  }

  object <- AddMetaData(object, metadata = df.cc)
  DefaultAssay(object) <- assay.old
  return(object)
}

#' @importFrom SeuratObject DefaultAssay Layers
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @keywords internal
#' @noRd
.prep_MetaDataBatch <- function(object, batch.var = NULL, assay = NULL,
                                layer = "data", scale.layer = NULL,
                                idcol = NULL) {
  assay <- assay %||% DefaultAssay(object)
  idcol <- idcol %||% 'cellbarcodeid'
  df.mtdt <- object[[]] %>% rownames_to_column(var = idcol)

  batch.var <- batch.var %||% {
    layers <- Layers(object, search = layer, assay = assay)
    df.mtdt <- df.mtdt %>% left_join(
      CreateIntegrationGroups(object = object[[assay]],
                              layers = layers,
                              scale.layer = scale.layer) %>%
        rownames_to_column(idcol), by = idcol
    )
    "group"
  }
  if (! batch.var %in% colnames(df.mtdt)) {
    abort(message = sprintf("%s not in colnames of metadata", sQuote(batch.var)))
  }
  list2env(list(batch.var = batch.var, df.mtdt = df.mtdt, idcol = idcol),
                  envir = parent.frame())
}
