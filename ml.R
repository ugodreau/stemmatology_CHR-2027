## exhaustive_ml.R
## Exhaustive maximum-likelihood tree search for small taxon sets (8 taxa: 10,395 unrooted topologies).
## Requires: phangorn, ape, parallel. mclapply forks, so ncores > 1 works on Linux/macOS only.
##
## Two stages:
##   1. Screen every topology with a quick likelihood at fixed (equal) branch lengths.
##   2. Optimise branch lengths (optim.pml, topology fixed) for the best `screen_k` topologies,
##      plus any extra trees you supply (e.g. your NNI result), and keep the best.
## With screen_k = Inf stage 1 is skipped and ALL topologies are optimised (slower, fully exhaustive).
## Validate the screening on a few traditions: exhaustive_ml(d, screen_k = Inf)$best_ll should equal
## exhaustive_ml(d, screen_k = 300)$best_ll.

library(phangorn)
library(ape)

exhaustive_ml <- function(d, ncores = 4, screen_k = 300, tie_tol = 1e-4,
                          extra_trees = list(), start_length = 0.1) {
  n <- length(d)
  if (n > 9) stop("Exhaustive search is only practical up to 9 taxa (got ", n, ")")
  
  ## All unrooted binary topologies, with equal starting branch lengths
  trees <- allTrees(n, rooted = FALSE, tip.label = names(d))
  trees <- lapply(trees, function(tr) ape::compute.brlen(tr, start_length))
  n_topologies <- length(trees)
  
  ## Stage 1: quick screen at fixed branch lengths
  if (is.finite(screen_k) && screen_k < n_topologies) {
    screen_ll <- unlist(parallel::mclapply(trees, function(tr) pml(tr, d)$logLik,
                                           mc.cores = ncores))
    keep  <- order(screen_ll, decreasing = TRUE)[seq_len(screen_k)]
    trees <- trees[keep]
  }
  
  ## Always include the trees supplied by the caller (made binary and unrooted)
  extra <- lapply(extra_trees, function(tr) {
    ape::compute.brlen(ape::multi2di(ape::unroot(tr)), start_length)
  })
  cand <- c(trees, extra)
  
  ## Stage 2: optimise branch lengths on each candidate topology (topology stays fixed)
  fits <- parallel::mclapply(cand, function(tr) {
    fit <- try(optim.pml(pml(tr, d), optEdge = TRUE,
                         control = pml.control(trace = 0)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    list(ll = fit$logLik, tree = fit$tree)
  }, mc.cores = ncores)
  fits <- Filter(Negate(is.null), fits)
  
  ll   <- vapply(fits, function(f) f$ll, numeric(1))
  best <- which.max(ll)
  tied <- which(ll >= ll[best] - tie_tol)   # topologies tied with the best (zero-length edges)
  
  list(best_ll      = ll[best],
       best_tree    = fits[[best]]$tree,
       tied_trees   = structure(lapply(fits[tied], `[[`, "tree"), class = "multiPhylo"),
       n_tied       = length(tied),
       n_evaluated  = length(cand),
       n_topologies = n_topologies)
}

## ---- Usage inside the point-estimate loop -------------------------------------------------
## d_k: alignment recoded with the same number of states as BEAST (levels = 0:(k - 1)).
## ex <- exhaustive_ml(d_k, ncores = 8, screen_k = 300, extra_trees = list(fit_opt$tree))
## ml_exh_tree <- ex$best_tree          # one best tree
## # if several topologies tie, a conservative alternative is their strict consensus:
## # ml_exh_tree <- if (ex$n_tied > 1) ape::consensus(ex$tied_trees, p = 1) else ex$best_tree
## nni_missed <- ex$best_ll - fit_opt$logLik   # >= 0: how much the NNI search left on the table