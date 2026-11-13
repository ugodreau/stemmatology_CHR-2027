## evaluate_trees.R
## Accuracy / resolution / recall evaluation of inferred trees against a ground truth.
## Requires: ape, TreeTools, TreeDist, Quartet (ggplot2 for the plots only).
##
## Definitions (clustering information, TreeDist):
##   M  = mutual clustering information of pred and gt
##   Hp = clustering entropy of pred,  Hg = clustering entropy of gt
##   accuracy   = M / Hp   (share of the information claimed by pred that is correct)
##   resolution = Hp / Hg  (information claimed by pred relative to the truth; > 1 = over-resolved)
##   recall     = M / Hg   (= accuracy * resolution; share of the truth recovered)
## Quartet outputs are kept as a lightweight secondary view:
##   quartet_div      = d / (s + d)   (error among quartets resolved in both trees)
##   quartet_resolved = share of all quartets resolved in pred

## ---- Which QuartetStatus column means "resolved in gt only"? ----------------
## Detected at load time, so the script does not depend on the r1/r2 convention.
.quartet_cols <- local({
  resolved <- ape::read.tree(text = "((a,b),(c,d));")
  star     <- ape::read.tree(text = "(a,b,c,d);")
  qs <- Quartet::QuartetStatus(star, cf = resolved)  # the only quartet is resolved in cf only
  gt_only <- if (qs[1, "r1"] > 0) "r1" else "r2"
  list(gt_only = gt_only, pred_only = setdiff(c("r1", "r2"), gt_only))
})

.metric_names <- c("quartet_div", "quartet_resolved", "cid", "accuracy",
                   "resolution", "recall", "M", "Hp", "Hg")

## ---- Score one prediction against the ground truth --------------------------
## pred: phylo, multiPhylo (scores are averaged over trees) or NULL (failed run -> NA row)
## collapse_tol: NULL (default) = keep the tree as is; a number = collapse internal edges
##   of length <= tol into polytomies (note: this also collapses negative NJ edges).
evaluate_tree <- function(pred, gt, collapse_tol = NULL) {
  na_row <- setNames(rep(NA_real_, length(.metric_names)), .metric_names)
  if (is.null(pred) || is.null(gt)) return(na_row)
  
  if (inherits(pred, "multiPhylo")) {
    rows <- vapply(pred, evaluate_tree, numeric(length(.metric_names)),
                   gt = gt, collapse_tol = collapse_tol)
    out <- rowMeans(rows, na.rm = TRUE)   # NA trees (e.g. stars) are ignored in the mean
    out[is.nan(out)] <- NA_real_
    return(setNames(out, .metric_names))
  }
  
  if (!setequal(pred$tip.label, gt$tip.label)) {
    warning("Tip labels differ between prediction and ground truth: returning NA")
    return(na_row)
  }
  if (ape::is.rooted(pred)) pred <- ape::unroot(pred)
  if (ape::is.rooted(gt))   gt   <- ape::unroot(gt)
  pred <- TreeTools::RenumberTips(pred, gt$tip.label)
  
  if (!is.null(collapse_tol) && !is.null(pred$edge.length)) {
    pred <- ape::di2multi(pred, tol = collapse_tol)
  }
  
  ## Clustering-information framework
  M   <- as.numeric(TreeDist::MutualClusteringInfo(pred, gt))
  Hp  <- as.numeric(TreeDist::ClusteringEntropy(pred))
  Hg  <- as.numeric(TreeDist::ClusteringEntropy(gt))
  cid <- as.numeric(TreeDist::ClusteringInfoDistance(pred, gt, normalize = TRUE))
  if (Hp + Hg > 0 &&
      !isTRUE(all.equal(cid, (Hp + Hg - 2 * M) / (Hp + Hg), tolerance = 1e-6))) {
    warning("Normalized CID differs from (Hp + Hg - 2M) / (Hp + Hg): check TreeDist's normalization",
            call. = FALSE)
  }
  
  ## Quartets (secondary view)
  qs <- Quartet::QuartetStatus(pred, cf = gt)
  s <- as.numeric(qs[1, "s"]); d <- as.numeric(qs[1, "d"])
  pred_only <- as.numeric(qs[1, .quartet_cols$pred_only])
  n_quartets <- choose(length(gt$tip.label), 4)
  
  c(quartet_div      = as.numeric(
    Quartet::QuartetDivergence(qs, similarity = FALSE)
    ), #TODO: fix in paper
    quartet_resolved = (s + d + pred_only) / n_quartets,
    cid              = cid,
    accuracy         = if (Hp > 0) min(M / Hp, 1) else NA_real_,   # clamp rounding error
    resolution       = if (Hg > 0) Hp / Hg else NA_real_,
    recall           = if (Hg > 0) min(M / Hg, 1) else NA_real_,
    M = M, Hp = Hp, Hg = Hg)
}

## ---- Score all methods for one tradition -> long data frame -----------------
## preds: named list, method name -> phylo / multiPhylo / NULL
score_methods <- function(preds, gt, tradition, collapse_tol = NULL) {
  rows <- lapply(names(preds), function(m) {
    r <- evaluate_tree(preds[[m]], gt, collapse_tol)
    data.frame(tradition = tradition, method = m, t(r), row.names = NULL)
  })
  do.call(rbind, rows)
}

## ---- Summary per method: medians + number of undefined accuracies -----------
summarise_methods <- function(df,
                              metrics = c("quartet_div", "quartet_resolved", "cid",
                                          "accuracy", "resolution", "recall")) {
  do.call(rbind, lapply(split(df, df$method), function(d) {
    meds <- vapply(metrics, function(m) median(d[[m]], na.rm = TRUE), numeric(1))
    data.frame(method = d$method[1], n = nrow(d),
               n_NA_accuracy = sum(is.na(d$accuracy)),   # e.g. star-shaped predictions
               as.list(setNames(meds, paste0("median_", metrics))),
               row.names = NULL)
  }))
}

## ---- Paired comparison of two methods on one metric -------------------------
paired_diff <- function(df, metric, a, b) {
  wide <- merge(df[df$method == a, c("tradition", metric)],
                df[df$method == b, c("tradition", metric)],
                by = "tradition", suffixes = c("_a", "_b"))
  d <- wide[[paste0(metric, "_a")]] - wide[[paste0(metric, "_b")]]
  d <- d[!is.na(d)]
  list(n = length(d), median_diff = median(d),
       wins_a = sum(d > 0), ties = sum(d == 0), wins_b = sum(d < 0),
       p_wilcoxon = suppressWarnings(wilcox.test(d)$p.value))
}



## ---- Accuracy (y) vs resolution (x), with iso-recall curves -----------------
## df: output of score_methods(); if it has a `threshold` column, medians are
## joined into one path per method (e.g. bootstrap consensus thresholds).
plot_accuracy_resolution <- function(df, recall_levels = c(.2, .4, .6, .8)) {
  xmax <- max(1.5, max(df$resolution, na.rm = TRUE))
  iso <- do.call(rbind, lapply(recall_levels, function(r) {
    res <- seq(r, xmax, length.out = 200)          # accuracy <= 1 requires resolution >= recall
    data.frame(recall = r, resolution = res, accuracy = r / res)
  }))
  med <- aggregate(cbind(resolution, accuracy) ~ method, data = df, FUN = median)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = resolution, y = accuracy, colour = method)) +
    ggplot2::geom_line(data = iso, ggplot2::aes(x = resolution, y = accuracy, group = recall),
                       colour = "grey80", inherit.aes = FALSE) +
    ggplot2::geom_point(alpha = .15, na.rm = TRUE) +
    ggplot2::stat_ellipse(level = .5, na.rm = TRUE) +
    ggplot2::geom_point(data = med, size = 3) +     # 2D medians per method
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Resolution (Hp / Hg)", y = "Accuracy (M / Hp)")
  
  if ("threshold" %in% names(df)) {
    path <- aggregate(cbind(resolution, accuracy) ~ method + threshold, data = df, FUN = median)
    path <- path[order(path$method, path$threshold), ]
    p <- p + ggplot2::geom_path(data = path,
                                ggplot2::aes(x = resolution, y = accuracy, colour = method),
                                inherit.aes = FALSE)
  }
  p
}

## ---- Example of use inside the simulation loop ------------------------------
## results <- NULL
## for (i in seq_along(traditions)) {
##   ...
##   preds <- list(NJ       = treeNJ,
##                 UNJ      = treeUNJ,
##                 MP_SC    = ape::consensus(mm_tree, p = 1),
##                 MP_mean  = mm_tree,               # multiPhylo: scores are averaged
##                 ML_NNI   = fit_opt_NNI$tree,
##                 ML_Ratch = fit_opt_ratch$tree,    # NULL if the ratchet failed -> NA row
##                 ML_NNI_BS = fit_opt_NNIBS)
##   results <- rbind(results, score_methods(preds, gt, tradition = i))
## }
## summarise_methods(results)
## paired_diff(results, "cid", "ML_NNI", "MP_SC")
## plot_accuracy_resolution(results)

## Column name in the results table -> label shown on the x axis, in display order.
## Columns absent from the data are skipped, so commented-out methods cause no error.
default_methods <- c(
  NJ        = "NJ",
  UNJ       = "UNJ",
  MP        = "MP",
  NJ_BS     = "NJ",
  UNJ_BS    = "UNJ",
  MP_BS     = "MP",
  ML_BS     = "ML",
  MP_mean   = "MP (mean)",
  MP_SC     = "MP (SC)",
  MP_IW_SC  = "MP (IW, SC)",
  ML        = "ML",
  FBD_MCC   = "FBD (MCC)",
  FBD_MRC   = "FBD (MRC)"
  # , FBD   = "BDSky"
)


## Shared drawing code. `long` has a factor `method` and a numeric `value`.
## The number of NA values per method is added to its label, since boxplots
## silently drop them (which can flatter methods that sometimes fail).
.boxplot_core <- function(long, ylab, title, subtitle) {
  n_na <- tapply(is.na(long$value), long$method, sum)
  n_na[is.na(n_na)] <- 0
  lv <- levels(long$method)
  levels(long$method) <- paste0(lv, ifelse(n_na[lv] > 0, paste0("\n(", n_na[lv], " NA)"), ""))
  
  ggplot2::ggplot(long, ggplot2::aes(x = method, y = value)) +
    ggplot2::geom_boxplot(na.rm = TRUE) +
    ggplot2::labs(x = "Method", y = ylab, title = title, subtitle = subtitle)
}

## ---- Wide input: one column per method (your results_grf / results_quartet) ----
plot_method_boxplots <- function(res, ylab, title = NULL, subtitle = NULL,
                                 methods = default_methods) {
  res <- as.data.frame(res)
  methods <- methods[names(methods) %in% names(res)]
  long <- data.frame(
    method = factor(rep(unname(methods), each = nrow(res)), levels = unname(methods)),
    value  = unlist(res[names(methods)], use.names = FALSE)
  )
  .boxplot_core(long, ylab, title, subtitle)
}

## ---- Long input: output of score_methods() (columns `method` and a metric) ------
## metric: name of the column to plot, e.g. "cid", "quartet_div", "accuracy".
plot_metric_boxplots <- function(df, metric, ylab = metric, title = NULL,
                                 subtitle = NULL, methods = NULL) {
  df <- as.data.frame(df)
  if (is.null(methods)) methods <- setNames(unique(df$method), unique(df$method))
  methods <- methods[names(methods) %in% df$method]
  df <- df[df$method %in% names(methods), ]
  long <- data.frame(
    method = factor(unname(methods[df$method]), levels = unname(methods)),
    value  = df[[metric]]
  )
  .boxplot_core(long, ylab, title, subtitle)
}



