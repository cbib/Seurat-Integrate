#' @include numeric_lisi_class.R
#' @include import-standalone-types-check.R
NULL

#' Ensures a Seurat object has 'si_scores' in Misc
#' @description
#' Checks that a Seurat object has 'si_scores' in Misc, creates an empty score
#' tibble if not, and returns the Seurat object. Does not create columns for
#' LISI scores, they are created on the run.
#'
#' @importFrom SeuratObject Misc Misc<-
#' @importFrom tibble tibble
#' @keywords internal
#' @noRd
check_misc <- function(object) {
  if (Misc(object, slot = 'si_scores') %iff% F %||% T) {
    integrations <- c('Unintegrated', 'CCA', 'RPCA', 'Harmony', 'MNN', 'FastMNN',
                      'BBKNN', 'ComBat', 'Scanorama_reduction', 'Scanorama_counts',
                      'scVI', 'scANVI', 'trVAE')
    empty_dbl <- rep(NA_real_, length(integrations))
    empty_list <- sapply(integrations, function(x) NA, simplify = F, USE.NAMES = F)
    empty_scores <- tibble(
      Integration = integrations,
      PCA.regression = empty_dbl,
      PCA.density = empty_dbl,
      cell.cycle.conservation = empty_list,
      ASW = empty_dbl,
      ASW.batch = empty_dbl,
      NMI = empty_dbl,
      ARI = empty_dbl,
      Graph.connectivity = empty_dbl,
      scGraph = empty_dbl
    )
    slot(object = object, name = 'misc')[['si_scores']] <- empty_scores
  }
  return(object)
}


#' Retrieve integration scores from a Seurat object
#' @description
#' Scores are stored as tibble objects in Misc. Slot for raw scores is named
#' 'si_scores' while scaled scores are found under 'si_scaled.scores'.
#'
#' \code{GetMiscIntegrations}: Get (or search for) integration names in the
#' score tibble
#'
#' \code{GetMiscScores}:  Get (or search for) score names in the score tibble
#'
#' \code{IntegrationScores}: Get the tibble with scaled or unscaled scores
#'
#' @param object a Seurat object
#' @param search a character vector of names to search for through
#' case-insensitive exact match. \code{NULL} disables the search (default)
#' @param scaled whether to query the unscaled (default) or scaled scores (if
#' scaling has been performed)
#'
#' @return \code{GetMiscIntegrations}: a character vector of integration names,
#' or \code{NULL} when the object doesn't have 'si_scores' in Misc or when
#' search did not return any result.
#'
#' \code{GetMiscScores}: a character vector of score names, or \code{NULL} when
#' the object doesn't have 'si_scores' in Misc or when search did not return any
#' result.
#'
#' \code{IntegrationScores}: either \code{NULL} if the requested object does not
#' exist, otherwise a tibble. The first column contains the name of the
#' integrations, and each following column corresponds to a score.
#'
#' @importFrom SeuratObject Misc
#' @export
#' @rdname get-score

GetMiscIntegrations <- function(object, search = NULL) {
  integrations <- as.character(Misc(object = object, slot = "si_scores")$Integration)
  search <- as.character(search %||% integrations)
  idx <- which(tolower(integrations) %in% tolower(search))
  integrations <- integrations[idx]
  if (length(integrations) == 0) {
    integrations <- NULL
  }
  return(integrations)
}

#' @importFrom SeuratObject Misc
#' @export
#' @rdname get-score
GetMiscScores <- function(object, search = NULL) {
  scores <- colnames(Misc(object = object, slot = "si_scores"))
  search <- search %||% scores
  idx <- which(tolower(scores) %in% tolower(search))
  scores <- scores[idx]
  if (length(scores) == 0) {
    scores <- NULL
  }
  return(scores)
}

#' Add integration(s) or score(s) slot(s) to the score tibble
#' @description
#' \code{AddMiscIntegrations}: Add integration(s) slot(s) to the score tibble
#' (or make sure they are in it). This enables to save their scores later.
#'
#' @param object a Seurat object
#' @param which the name of the new integration(s) or score(s)
#'
#' @return the Seurat object with updated table of scores.
#'
#' @importFrom methods slot<- as
#' @importFrom SeuratObject Misc
#' @importFrom tibble add_row
#' @importFrom dplyr %>%
#' @importFrom rlang !!!
#' @name add-score
AddMiscIntegrations <- function(object, which) {
  integrations <- GetMiscIntegrations(object, which)
  which <- which[!tolower(which) %in% tolower(integrations)]
  if ((l <- length(which)) > 0) {
    slot(object = object, name = 'misc')[['si_scores']] <-
      Misc(object, 'si_scores') %>%
      add_row(Integration = { which }, !!! lapply(.[, -1], function(score.col) {
        rep(as(NA, class(score.col)), l)
      }))
  }
  return(object)
}

#' @description
#' \code{AddMiscScores}: Add score(s) slot(s) to the score tibble (or make sure
#' they are in it). This enables to save them later.
#'
#' @param class the class of the column, such as 'numeric' (default). For
#' instance, cell cycle conservation scores are 'list' and LISI scores are
#' 'numeric_lisi'.
#'
#' @importFrom SeuratObject Misc
#' @importFrom tibble add_column
#' @importFrom dplyr %>%
#' @importFrom rlang !!!
#' @rdname add-score
AddMiscScores <- function(object, which, class = "numeric") {
  class <- class %||% "numeric"
  scores <- GetMiscScores(object, which)
  which <- which[!tolower(which) %in% tolower(scores)]
  if (length(which) > 0) {
    l <- length(GetMiscIntegrations(object))
    empty_dbl <- as(rep(NA_real_, l), class)
    slot(object = object, name = 'misc')[['si_scores']] <-
      Misc(object, 'si_scores') %>%
      add_column(!!! sapply(which, function(new.col) empty_dbl, simplify = F))
  }
  return(object)
}


#' Set the value of a score in the score tibble
#' @description
#' Set a score value for a given integration and a given score name. If they
#' don't exist yet, they are created though \code{\link{AddMiscIntegrations}}
#' and \code{\link{AddMiscScores}} respectively.
#'
#' @param object a Seurat object
#' @param integration the name of the integration for which the score was
#' computed.
#' @param score.name the name of the computed score
#' @param score.value the value of the score
#' @param ... additional parameter to pass to \code{\link{AddMiscScores}}. Has
#' no effect if the score slot is already present in the tibble. Otherwise,
#' enable to pass the \code{class} argument to specify the class of the score
#' slot to create.
#'
#' @inherit add-score return
#'
#' @importFrom SeuratObject Misc
SetMiscScore <- function(object, integration, score.name, score.value, ...) {
  integration %||% abort('integration name cannot be null')
  if(length(integration) > 1) abort('cannot add score to multiple integrations')
  score.name %||% abort('score name cannot be null')
  if(length(score.name) > 1) abort('cannot add multiple scores')
  score.value %||% abort('score value cannot be null')

  object <- AddMiscIntegrations(object, which = integration)
  integration <- GetMiscIntegrations(object = object, search = integration)
  object <- AddMiscScores(object, which = score.name, ...)
  score.name <- GetMiscScores(object = object, search = score.name)

  scores.df <- Misc(object = object, slot = 'si_scores')

  # handle cc conservation score case
  if (inherits(scores.df[[score.name]], 'list') && ! inherits(score.value, 'list')) {
    score.value <- list(score.value)
  }
  scores.df[[score.name]][scores.df$Integration == integration] <- score.value

  slot(object = object, name = 'misc')[['si_scores']] <- scores.df
  return(object)
}


#' Compute ranks of a vector and rescale between 0 and 1
#' @description
#' Compute ranks of a vector with \code{NA}s kept in their original positions.
#' The ranks are then adjusted by adding half the number of code{NA} values.
#' Those adjusted ranks are then normalised by the length of the vector, such
#' that resulting values are bounded between 0 and 1 when there is no \code{NA}.
#' The range of the output values shrinks as the number of \code{NA}s increases.
#'
#' @param x a numerical vector
#' @keywords internal
#' @noRd
adj_rank_scaled <- function(x) {
  adj_rank <- rank(x, na.last = 'keep') + sum(is.na(x)) / 2
  return(adj_rank / length(x))
}

#' @keywords internal
#' @noRd
get.score.types <- function(col_names, batch = FALSE) {
  if(batch) {
    found.colnames <- c(
      col_names[tolower(col_names) %in% tolower(c('PCA.regression', 'PCA.density'))],
      col_names[grep("^Graph.connectivity", col_names, T)],
      col_names[grep("^ASW(?=[[:punct:][:blank:]]*batch)", col_names, T, T)],
      col_names[grep('^iLISI', col_names, T)],
      col_names[grep('^kBET', col_names, T)]
    )
  } else {
    found.colnames <- c(
      tolower(col_names)[tolower(col_names) %in% tolower('cell.cycle.conservation')],
      col_names[grep("^ARI|^NMI", col_names, T)],
      col_names[grep("^ASW(?![[:punct:][:blank:]]*batch)", col_names, T, T)],
      col_names[grep('^cLISI', col_names, T)],
      col_names[grep('^scGraph', col_names, T)]
    )
  }
  return(found.colnames)
}

#' @importFrom dplyr %>% rowwise mutate c_across ungroup case_when
#' @importFrom rlang data_syms !!!
#' @keywords internal
#' @noRd
compute.overall.scores <- function(scaled.scores, batch.scores, bio.scores,
                                   batch.coeff = .4, bio.coeff = .6) {
  scaled.scores <- scaled.scores %>% rowwise()
  l <- c(length(bio.scores), length(batch.scores))
  if (all(l == 0)) {
    scaled.scores <- scaled.scores %>%
      mutate(Overall.score = mean(c_across(!Integration), na.rm = T))
    warning("Could not discriminate between score types. Overall score is the mean of everything")
  } else {
    scaled.scores <- suppressWarnings(
      scaled.scores %>%
        mutate(Bio.conservation = mean(c(!!!data_syms(bio.scores)), na.rm = T),
               Batch.correction = mean(c(!!!data_syms(batch.scores)), na.rm = T)) %>%
        ungroup() %>% # remove rowwise
        mutate(Overall.score = case_when(
          is.na(Bio.conservation) ~ Batch.correction,
          is.na(Batch.correction) ~ Bio.conservation,
          T ~ batch.coeff * Batch.correction + bio.coeff * Bio.conservation)
        )
    )

    if (any(l == 0)) {
      i <- l == 0
      msg <- paste('Did not find any', c('bio-conservation', 'batch correction')[i],
                   'scores. Overall score is the mean of',
                   c('bio-conservation', 'batch correction')[!i], 'scores')
      warning(msg)
    } else {
      if (any((i <- is.na(scaled.scores$Bio.conservation)))) {
        msg <- paste(paste(sQuote(scaled.scores$Integration[i]), collapse = ', '),
                     'don\'t have bio-conservation scores.',
                     'Their overall score is the mean of batch correction scores')
        warning(msg)
      }
      if (any((i <- is.na(scaled.scores$Batch.correction)))) {
        msg <- paste(paste(sQuote(scaled.scores$Integration[i]), collapse = ', '),
                     'don\'t have batch correction scores.',
                     'Their overall score is the mean of bio-conservation scores')
        warning(msg)
      }
    }
  }
  return(scaled.scores)
}

#' Scale the scores in the score tibble to plot them
#' @description
#' Once the scores of interest have been computed and saved in the score tibble,
#' they can be scaled to make them comparable by bounding them between 0 and 1
#' and harmonise their direction (0 and 1 always mean bad and good performance
#' respectively). This is also a prerequisite for plotting.
#'
#' @param object a Seurat object
#' @param ref the name of the integration to use as a reference for scaling.
#' Useful for PCA regression (and density) and cell cycle conservation scores.
#' @param rescale whether to rescale each score between 0 and 1 using min-max
#' normalisation before computing overall scores. One of 'none' (default, no
#' rescaling), 'rank' or 'score', corresponding to no min-max normalisation of
#' ranks or scores respectively. Beware that the argument of
#' \code{\link{PlotScores}()} \strong{overrides the choice here} (see
#' \strong{Details} sections).
#' @param batch.coeff the weight of batch correction performance evaluation
#' scores in the overall score.
#' @param bio.coeff the weight of bio-conservation performance evaluation scores
#' in the overall score.
#'
#' @details
#' It is recommended to keep the default parameter \code{rescale = "none"} to
#' retain the possibility to plot scores without min-max rescaling later with
#' \code{\link{PlotScores}()}.
#'
#' @seealso \code{\link{PlotScores}()} for further details on rescaling
#'
#' @importFrom SeuratObject Misc
#' @importFrom dplyr %>% select arrange filter summarise across mutate bind_rows rowwise c_across ungroup case_when
#' @importFrom purrr map2 reduce
#' @importFrom rlang sym syms data_syms !! !!!
#' @importFrom scales rescale
#' @export
ScaleScores <- function(object, ref = "Unintegrated",
                        rescale = c("none", "rank", "score"),
                        batch.coeff = .4, bio.coeff = .6) {
  ref <- ref %||% "Unintegrated"
  Misc(object, slot = 'si_scores') %||%
    abort("No scores. Please compute scores before.")
  raw.scores <- Misc(object, slot = 'si_scores') %>%
    select(Integration, !Integration) # put Integration column 1st
  rescale <- rescale %||% "none"
  if (is.logical(rescale)) {
    rescale <- c("none", "rank")[rescale + 1]
  }
  rescale <- tolower(rescale)
  rescale <- match.arg(rescale)

  # ensure sum of coefficients is 1 (to keep overall scores between 0 and 1)
  sum.coeff <- batch.coeff + bio.coeff
  batch.coeff <- batch.coeff / sum.coeff
  bio.coeff <- bio.coeff / sum.coeff


  # remove cols when unintegrated is NA for scores that require reference to be scaled
  var <- c('PCA.regression', 'PCA.density', 'cell.cycle.conservation')
  if (! tolower(ref) %in% tolower(raw.scores$Integration)) {
    warning(paste(sQuote(var), collapse = ', '), ' require',
            ' a reference to scale other scores with. But provided reference ',
            dQuote(ref, '"'), ' has not been found. Cannot scale, skipping')
    raw.scores <- raw.scores[, setdiff(colnames(raw.scores), var), drop = F]
  } else {
    raw.scores <- raw.scores %>%
      arrange(tolower(Integration) != tolower(ref)) # put ref row 1st

    colsidx <- raw.scores %>%
      filter(tolower(Integration) == tolower(ref)) %>%
      summarise(across({{ var }}, is.na)) %>% unlist()
    if (any(colsidx)) {
      warning(paste(sQuote(names(colsidx)[colsidx]), collapse = ', '), ' require',
              ' a reference to scale other scores with. But ', dQuote(ref, '"'),
              ' scores have not been computed. Cannot scale, skipping')
    }
    colsidx <- names(colsidx)[colsidx]
    raw.scores <- raw.scores[, setdiff(colnames(raw.scores), colsidx), drop = F]
  }

  if (ncol(raw.scores) < 2) {
    abort("Nothing left to scale. Please compute scores before.")
  }

  # remove rows & cols with only nas (/!\ cc.conservation)
  colsidx <- colSums(is.na(select(raw.scores, !Integration))) < nrow(raw.scores)
  rowsidx <- rowSums(is.na(select(raw.scores, !Integration))) < ncol(raw.scores) - 1

  if (!any(rowsidx) | !any(colsidx)) {
    abort("Scores not computed, all missing (NAs)")
  }
  colsidx <- c("Integration", names(colsidx)[colsidx])
  raw.scores <- raw.scores[rowsidx, colsidx, drop = F]

  # remove cc.conservation if number of !NA < 2
  var <- 'cell.cycle.conservation'
  if (var %in% colnames(raw.scores)) {
    if (sum(!is.na(raw.scores[[var]])) < 2) {
      warning(sQuote(var), ' requires at least one score besides the ',
              'reference\'s score. Cannot scale, skipping')
      raw.scores <- raw.scores[, setdiff(colnames(raw.scores), var), drop = F]
    }
  }
  if (ncol(raw.scores) < 2) {
    abort("Nothing left to scale. Please compute scores before.")
  }

  identityy <- function(x, y) identity(x)
  scaling <- list(
    Integration = identityy,
    PCA.regression = function(x, y) pmax((x[1] - x) / x[1], 0),
    PCA.density = function(x, y) pmax((x - x[1]) / x, 0),
    cell.cycle.conservation = function(x, y) {
      i <- sapply(x, is.data.frame)
      res <- rep(NA, length(i))
      y <- y[i]
      res[i] <- x[i] %>%
        purrr::map2(y, ~ rename(.x, !!.y := score)) %>%
        purrr::reduce(left_join, by = colnames(.[[1]])[1]) %>%
        mutate(across(c(!!!syms(y)), ~ 1 - abs(.x - !!sym(y[1])) / !!sym(y[1]))) %>%
        summarise(across(c(!!!syms(y)), mean)) %>% t()
      res
    }
  )

  cLISI <- function(x, y) (N(x) - x) / (N(x) - 1)
  iLISI <- function(x, y) (x - 1) / (N(x) - 1)
  kBET <- function(x, y) 1 - x
  scGraph <- function(x, y) rescale(x, from = c(-1, 1), to = c(0, 1))

  scaling <- c(
    scaling,
    sapply(colnames(raw.scores)[grep("cLISI", colnames(raw.scores), T)], function(x) cLISI, simplify = F),
    sapply(colnames(raw.scores)[grep("iLISI", colnames(raw.scores), T)], function(x) iLISI, simplify = F),
    sapply(colnames(raw.scores)[grep("kBET", colnames(raw.scores), T)], function(x) kBET, simplify = F),
    sapply(colnames(raw.scores)[grep("scGraph", colnames(raw.scores), T)], function(x) scGraph, simplify = F))

  scaling <- c(
    scaling,
    sapply(colnames(raw.scores)[grep("^ARI|^NMI|^ASW", colnames(raw.scores), T)], function(x) identityy),
    sapply(colnames(raw.scores)[grep("^Graph.connectivity", colnames(raw.scores), T)], function(x) identityy)
  )

  scaling <- c(
    scaling, sapply(setdiff(colnames(raw.scores), names(scaling)), function(x) identityy)
  )

  scaled.scores <- raw.scores %>%
    purrr::map2(colnames(.), ~ scaling[[.y]](.x, raw.scores$Integration)) %>%
    bind_rows()

  if (rescale != "none") {
    warning(sprintf("only %s enables plotting non-rescaled scores with `%s`",
                    sQuote("rescale = 'none'"), "PlotScores()"),
            call. = TRUE, immediate. = TRUE)
    scaled.scores <- scaled.scores %>%
      mutate(across(!Integration, ~ switch(rescale,
                                           "rank" = adj_rank_scaled(x = .x),
                                           "score" = rescale(x = .x, to = c(0,1))
                                           )) )
  }

  bio.scores <- get.score.types(colnames(scaled.scores), batch = FALSE)
  batch.scores <- get.score.types(colnames(scaled.scores), batch = TRUE)

  scaled.scores <- compute.overall.scores(scaled.scores = scaled.scores,
                                          batch.scores = batch.scores,
                                          bio.scores = bio.scores,
                                          batch.coeff = batch.coeff,
                                          bio.coeff = bio.coeff)

  slot(object = object, name = 'misc')[['si_scaled.scores']] <- scaled.scores
  return(object)
}


#' @export
#' @importFrom SeuratObject Misc
#' @rdname get-score
IntegrationScores <- function(object, scaled = FALSE) {
  slot <- c('si_scores', 'si_scaled.scores')[scaled + 1]
  return(Misc(object, slot = slot))
}


#' Visualise and compare the performances of integration algorithms
#' @description
#' Plot the scaled integration scores to compare the obtained integrations
#'
#' @inheritParams ScaleScores
#' @param plot.type one of 'table' (default), 'radar' or 'lollipop'. Type of
#' desired plot
#' @param split.by.score.type whether to split scores by type (bio-conservation,
#' batch correction and overall scores). When set to \code{FALSE}, all scores
#' are mixed in a single figure.
#' @param order.by one of 'score' (default), 'name' or 'asis'. Determines the
#' order of integrations in the legend (and on the y-axis for lolliplop and
#' table plots). Scores are ordered by decreasing overall score by default,
#' by name or by row-order when setting 'name' and 'asis' respectively.
#' @param hide.zeros whether zero(ed) scores should be visible on the plot.
#' @param include.integration name of the integration(s) to include. The
#' default value (\code{NULL}) enable to include them all.
#' @param exclude.integration name of the integration(s) to exclude. The default
#' value (\code{NULL}) enable to include them all.
#' @param include.score name of the score(s) to include. The default value
#' (\code{NULL}) enable to include them all.
#' @param exclude.score name of the score(s) to exclude. The default value
#' (\code{NULL}) enable to include them all.
#' @param recompute.overall.scores whether to recompute overall scores. Useful
#' when some scores or integrations are excluded. When \code{FALSE},
#' coefficient parameters have no impact.
#' @param rescale whether to rescale each score between 0 and 1 using min-max
#' normalisation before computing overall scores. One of 'rank' (default),
#' 'score' or 'none'. The first two enable rescaling on ranks and scores
#' respectively, while 'none' disables rescaling. Rescaling ensures that each
#' metric equally contributes to the overall scores. 'rank' is more stable to
#' changes in integrations included. Has no effect when
#' \code{recompute.overall.scores = FALSE}.
#' @param point.max.size inoperative unless \code{plot.type = 'table'} and
#' \code{use.ggforce = FALSE}. Determine the maximum size of the points
#' (only achieved for a score of 1) to fit the plotting area (handled
#' automaticaly when ggforce is used).
#' @param use.ggforce for \code{plot.type = 'table'}, enable or disable the use
#' of \pkg{ggforce} to draw the circles. Used by default when the package is
#' installed
#'
#' @return a ggplot object
#'
#' @details
#' The \code{rescale} parameter controls how scores are rescaled prior to
#' computing overall scores. Thus, when \code{rescale != 'none'},
#' \code{recompute.overall.scores} should be \code{TRUE} (default). Otherwise,
#' the overall scores will no longer match the scores displayed on the plot.
#'
#' Albeit the min-max rescaling of scores in Luecken M.D. \emph{et al.}, 2022 is
#' achievable with \code{rescale = 'score'}, min-max rescaling on ranks
#' (\code{rescale = 'rank'}) is the default. Indeed, it results in more stable
#' rankings of overall scores when the list of integrations to plot is
#' modified.
#'
#' To plot non-rescaled scores, make sure to use
#' \code{\link{ScaleScores}(..., rescale = "none")} beforehand. To plot rescaled
#' scores however, the \code{rescale} parameter of \code{\link{ScaleScores}} has
#' no impact on the scores obtained for plotting.
#'
#' @seealso \code{\link{ScaleScores}()}
#'
#' @references Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A.,
#' Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M.,
#' Colomé-Tatché, M. & Theis, F. J. Benchmarking atlas-level data integration in
#' single-cell genomics. Nat Methods 19, 41–50 (2021).
#' \href{https://doi.org/10.1038/s41592-021-01336-8}{DOI}
#'
#' @importFrom rlang is_installed !!
#' @importFrom SeuratObject Misc
#' @importFrom dplyr %>% mutate across where case_when desc select filter
#' @importFrom scales rescale
#' @importFrom forcats fct_reorder fct_relabel
#' @importFrom tidyr pivot_longer
#' @include numeric_lisi_class.R
#' @export
PlotScores <- function(object, plot.type = c('dot', 'radar', 'lollipop'),
                       split.by.score.type = TRUE,
                       order.by = c('score', 'name', 'asis'),
                       hide.zeros = FALSE,
                       include.integration = NULL,
                       exclude.integration = NULL,
                       include.score = NULL,
                       exclude.score = NULL,
                       recompute.overall.scores = TRUE,
                       rescale = c("rank", "score", "none"),
                       batch.coeff = .4, bio.coeff = .6,
                       point.max.size = 20L,
                       use.ggforce = is_installed('ggforce')) {
  Misc(object, slot = 'si_scaled.scores') %||% abort('Scale scores first')
  plot.type <- tolower(plot.type)
  plot.type <- match.arg(plot.type)
  order.by <- tolower(order.by %||% 'asis')
  order.by <- match.arg(order.by)
  rescale <- rescale %||% "none"
  if (is.logical(rescale)) {
    rescale <- c("none", "rank")[rescale + 1]
  }
  rescale <- tolower(rescale)
  rescale <- match.arg(rescale)

  scaled.scores <- Misc(object, slot = 'si_scaled.scores')

  exclude.integration <- if(isFALSE(exclude.integration)) NULL else tolower(exclude.integration)
  include.integration <- if(isTRUE(include.integration)) NULL else include.integration
  include.integration <- setdiff(
    tolower(include.integration %||% scaled.scores$Integration),
    exclude.integration)

  exclude.score <- if(isFALSE(exclude.score)) NULL else tolower(exclude.score)
  include.score <- if(isTRUE(include.score)) NULL else include.score
  include.score <- setdiff(
    tolower(include.score %||% colnames(scaled.scores)),
    exclude.score)
  include.score <- unique(c("Integration", include.score, "Overall.score"))

  scaled.scores <- scaled.scores %>%
    select(which(tolower(colnames(.)) %in% include.score)) %>%
    filter(tolower(Integration) %in% !!include.integration)
  if(nrow(scaled.scores) < 1) {
    abort(paste('All integrations filtered out. Nothing left to plot.',
                'Consider less harsh exclusion or broader inclusion criteria'))
  }
  if(ncol(scaled.scores) < 2) {
    abort(paste('All scores filtered out. Nothing left to plot.',
                'Consider less harsh exclusion or broader inclusion criteria'))
  }

  if (rescale != 'none') {
    scaled.scores <- scaled.scores %>%
      mutate(across(!Integration, ~ switch(rescale,
                                           "rank" = adj_rank_scaled(x = .x),
                                           "score" = rescale(x = .x, to = c(0,1))
                                           )) )
  }

  bio.scores <- get.score.types(colnames(scaled.scores), batch = FALSE)
  if (recompute.overall.scores) {
    # ensure sum of coefficients is 1 (to keep overall scores between 0 and 1)
    sum.coeff <- batch.coeff + bio.coeff
    batch.coeff <- batch.coeff / sum.coeff
    bio.coeff <- bio.coeff / sum.coeff

    batch.scores <- get.score.types(colnames(scaled.scores), batch = TRUE)
    scaled.scores <- compute.overall.scores(scaled.scores = scaled.scores,
                                            batch.scores = batch.scores,
                                            bio.scores = bio.scores,
                                            batch.coeff = batch.coeff,
                                            bio.coeff = bio.coeff)
  }

  pretty_string <- function(x) {
    gsub(pattern = '[\\._]+', replacement = ' ', x = x)
  }

  lo.cap <- c(0, NA)[hide.zeros + 1]
  scaled.scores <- scaled.scores %>%
    mutate(Integration = if(order.by == "asis") {
      factor(Integration, levels = unique(Integration))
    } else {
      factor(Integration)
    }) %>%
    mutate(across(where(is.numeric_lisi), as.numeric)) %>%
    # cap between [0, 1], or ]0,1] when plot is a table (NA hides circles)
    mutate(across(!Integration, ~ case_when(.x <= 0 ~ !!lo.cap, .x >= 1 ~ 1, T ~ .x)))
  if (order.by == 'score') {
    scaled.scores <- scaled.scores %>%
      mutate(Integration = fct_reorder(Integration, desc(Overall.score)))
    # .desc = plot.type != 'lollipop'))
  }

  scaled.scores <- scaled.scores %>%
    pivot_longer(!Integration, names_to = 'Score', values_to = 'y', ) %>%
    mutate(Integration = fct_relabel(Integration, pretty_string)) %>%
    mutate(Score.type = case_when(
      Score %in% bio.scores ~ 'Bio-conservation',
      Score %in% c('Overall.score', 'Bio.conservation', 'Batch.correction') ~ 'Overall',
      T ~ 'Batch correction'
    )) %>% mutate(Score = factor(Score)) %>%
    mutate(Score = forcats::fct_relabel(Score, pretty_string))

  ggp <- switch(plot.type,
                'radar' = PlotScoresRadar(scaled.scores = scaled.scores,
                                          split.by.score.type = split.by.score.type),
                'lollipop' = PlotScoresLillipop(scaled.scores = scaled.scores,
                                                split.by.score.type = split.by.score.type),
                'dot' = PlotScoresTable(scaled.scores = scaled.scores,
                                          split.by.score.type = split.by.score.type,
                                          vanilla = !use.ggforce,
                                          point.max.size = point.max.size))
  ggp
}


#' @importFrom rlang arg_match
#' @importFrom ggplot2 ggproto ggplot aes geom_point geom_polygon theme element_line element_blank scale_y_continuous facet_grid CoordRadial
#' @importFrom dplyr %>% arrange
#' @importFrom cowplot theme_cowplot
#' @importFrom grid unit
#' @keywords internal
#' @noRd
PlotScoresRadar <- function(scaled.scores, split.by.score.type) {
  coord_radar <- function (theta = c("x", "y"), start = 0, end = NULL,
                           expand = TRUE, direction = 1, clip = "off",
                           r.axis.inside = NULL, rotate.angle = FALSE,
                           inner.radius = 0) {
    theta <- tolower(theta)
    theta <- arg_match(theta)
    r <- c('x', 'y')[(theta == 'x') + 1]
    check_bool(r.axis.inside, allow_null = TRUE)
    check_bool(expand)
    check_bool(rotate.angle)
    check_number_decimal(start, allow_infinite = FALSE)
    check_number_decimal(end, allow_infinite = FALSE, allow_null = TRUE)
    check_number_decimal(inner.radius, min = 0, max = 1, allow_infinite = FALSE)
    end <- end %||% (start + 2 * pi)
    if (start > end) {
      n_rotate <- ((start - end)%/%(2 * pi)) + 1
      start <- start - n_rotate * 2 * pi
    }
    r.axis.inside <- r.axis.inside %||% !(abs(end - start) >=
                                            1.999 * pi)
    ggproto(NULL, CoordRadial, theta = theta, r = r, arc = c(start, end),
            expand = expand, direction = sign(direction),
            r_axis_inside = r.axis.inside, rotate_angle = rotate.angle,
            inner_radius = c(inner.radius, 1) * 0.4, clip = clip,
            is_linear = function(coord) TRUE)
  }

  RADAR <- scaled.scores %>% arrange(Score) %>%
    ggplot(aes(x = Score, y = y, colour = Integration, group = Integration)) +
    geom_point(size = 2.5) +
    geom_polygon(fill = NA, show.legend = F, linewidth = 1.2) +
    theme_cowplot() +
    coord_radar(start = - .25 * pi, end = 1.75 * pi, clip = 'off', expand = T,
                r.axis.inside = T, rotate.angle = T, inner.radius = .05) +
    theme(panel.grid.major.y = element_line(colour = 'grey70', linetype = 3, linewidth = 1),
          panel.grid.major.x = element_line(colour = 'grey60', linewidth = 1),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.r = element_blank(),
          axis.ticks.theta = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          plot.margin = unit(rep(2, 4), "lines")) +
    scale_y_continuous(breaks = c(0,.25, .5, .75, 1),
                       limits = c(0,1),
                       expand = rep(0,4))
  if (split.by.score.type) {
    RADAR <- RADAR + facet_grid(~ Score.type, scales = 'free_x')
  }
  RADAR
}


#' @importFrom ggplot2 ggplot aes geom_point geom_errorbarh position_dodge theme element_line element_rect element_blank scale_x_continuous xlab facet_grid
#' @importFrom dplyr %>%
#' @importFrom cowplot theme_cowplot
#' @keywords internal
#' @noRd
PlotScoresLillipop <- function(scaled.scores, split.by.score.type) {
  LOLLIPOP <- scaled.scores %>%
    ggplot(aes(x = y, y = Score, colour = Integration, group = Integration)) +
    geom_point(size = 2.5, position = position_dodge(width = .7)) +
    geom_errorbarh(aes(xmin = 0, xmax = y),
                   height = 0,
                   position = position_dodge(width = .7)) +
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(colour = 'grey70', linetype = 2, linewidth = 1),
          panel.grid.minor.x = element_line(colour = 'grey80', linetype = 3, linewidth = .8),
          strip.background = element_rect(colour="black", fill="white",
                                          linewidth=1.5, linetype="solid"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) +
    # guides(colour = guide_legend(reverse = TRUE)) +
    scale_x_continuous(breaks = c(0, .5, 1), limits = c(0,1), expand = c(0, 0, .1, 0),
                       minor_breaks = c(.25, .75)) +
    xlab("Score")
  if (split.by.score.type) {
    LOLLIPOP <- LOLLIPOP + facet_grid(Score.type ~ ., scales = 'free_y', space = 'free')
  }
  LOLLIPOP
}


#' @importFrom ggplot2 ggplot aes geom_point scale_size_area geom_rect guides theme element_blank element_rect element_text element_line facet_grid
#' @importFrom dplyr %>% mutate group_by
#' @importFrom forcats fct_rev
#' @importFrom cowplot theme_cowplot
#' @keywords internal
#' @noRd
PlotScoresTable <- function(scaled.scores, split.by.score.type, vanilla, point.max.size = 20) {
  TABLE <- scaled.scores %>%
    mutate(Integration = fct_rev(Integration)) %>%
    mutate(RECT_y = as.integer(droplevels(Integration))) %>%
    group_by(if (split.by.score.type) {Score.type} else {NULL}) %>%
    mutate(RECT_x = as.integer(droplevels(Score))) %>%
    ggplot(aes(x = Score, y = Integration, group = Integration))
  TABLE <- if (vanilla) {
    TABLE + geom_point(aes(colour = y, size = y), na.rm = T) +
      scale_size_area(max_size = point.max.size)
  } else {
    TABLE + ggforce::geom_circle(aes(x0 = RECT_x, y0 = RECT_y, fill = y, r = y / 2.3),
                                 color = "#808080", linewidth = 2, na.rm = T)
  }
  TABLE <- TABLE +
    geom_rect(aes(xmin = RECT_x - .5, xmax = RECT_x + .5,
                  ymin = RECT_y - .5, ymax = RECT_y + .5), colour = "black",
              fill = NA, inherit.aes = T) +
    guides(size = 'none') +
    theme_cowplot() +
    theme(axis.line = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 1),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          panel.grid.minor = element_line(colour = "black", linewidth = 1),
          strip.background = element_rect(colour="black", fill="white",
                                          linewidth=1.5, linetype="solid"))
  if (split.by.score.type) {
    TABLE <- TABLE + facet_grid(~ Score.type, scales = 'free_x', space = 'free_x')
  }
  TABLE + fullsize_colorbar(vanilla = vanilla)
}

fullsize_colorbar <- function(vanilla = T) structure(list(vanilla.circles = vanilla), class = "fullsizebar")

#' @method ggplot_add fullsizebar
#' @exportS3Method ggplot2::ggplot_add
#' @importFrom ggplot2 ggplot_add ggplotGrob scale_fill_distiller scale_colour_distiller guide_colorbar theme element_text
#' @importFrom grid unitType unit
#' @keywords internal
#' @noRd
ggplot_add.fullsizebar <- function(object, plot, object_name = "fullsizebar", ...) {
  h <- ggplotGrob(plot)$heights
  panel <- which(unitType(h) == "null")
  panel_height <- unit(1, "npc") - sum(h[-panel])

  scale_distiller <- scale_fill_distiller
  if (object$vanilla.circles) {
    scale_distiller <- scale_colour_distiller
  }

  plot +
    scale_distiller(palette = 'RdBu', limits = c(0, 1),
                    guide = guide_colorbar(barheight = panel_height,
                                           title.position = "right",
                                           title = "Score", nbin = 1000,
                                           ticks.colour = "black",
                                           frame.colour = "black",
                                           barwidth = 1.5)) +
    theme(legend.title = element_text(angle = -90, hjust = 0.5))
}
