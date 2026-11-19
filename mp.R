## mp_exact.R
## Exact parsimony searches for small taxon sets (<= 9 taxa), with equal or implied weights.
## Requires: ape, phangorn, TreeSearch (implied weights only).
##
## Design: enumerate all unrooted topologies once (10,395 for 8 taxa), then score them with the SAME
## engine under every weighting scheme. This removes the confound of comparing different search
## algorithms, and the optimal sets are exact (no missed tied trees).

library(ape)
library(phangorn)

as_multiPhylo <- function(x) {
  if (inherits(x, "multiPhylo")) x else structure(list(x), class = "multiPhylo")
}

## ---- 1. All unrooted binary topologies (practical up to 9 taxa) ---------------------------------
all_topologies <- function(data) {
  n <- length(data)
  if (n > 9) stop("Enumerating all topologies is only practical up to 9 taxa (got ", n, ")")
  phangorn::allTrees(n, rooted = FALSE, tip.label = names(data))
}

## ---- 2. Equal weights: all most parsimonious trees by enumeration -------------------------------
mpt_enumerate <- function(data, trees = all_topologies(data)) {
  sc  <- phangorn::parsimony(trees, data)        # Fitch score of every topology
  out <- trees[sc == min(sc)]
  attr(out, "score") <- min(sc)
  out
}


## Sanity check: two sets of trees contain the same topologies (compares Robinson-Foulds = 0)
same_tree_sets <- function(a, b) {
  a <- as_multiPhylo(a); b <- as_multiPhylo(b)
  covers <- function(x, y) all(vapply(x, function(t)
    any(vapply(y, function(u) phangorn::RF.dist(t, u) == 0, logical(1))), logical(1)))
  length(a) == length(b) && covers(a, b) && covers(b, a)
}

## ---- 4. Collapse branches that carry no change (zero length under ACCTRAN) ----------------------
## Removes spurious resolution inside each tree. TNT's collapsing rule is based on the minimum
## possible length, so this is an approximation; check on a few trees before relying on it.
collapse_zero <- function(trees, data, tol = 1e-8) {
  trees <- as_multiPhylo(trees)
  out <- lapply(trees, function(tr) {
    rooted <- ape::root(tr, outgroup = 1, resolve.root = TRUE)   # acctran needs a rooted binary tree
    with_len <- try(phangorn::acctran(rooted, data), silent = TRUE)
    if (inherits(with_len, "try-error")) return(tr)
    ape::unroot(ape::di2multi(with_len, tol = tol))
  })
  structure(out, class = "multiPhylo")
}

## ---- 5. Consensus summaries of a set of trees ----------------------------------------------------
consensus_set <- function(trees) {
  trees <- as_multiPhylo(trees)
  if (length(trees) == 1) return(list(strict = trees[[1]], majority = trees[[1]], n_trees = 1L))
  list(strict   = ape::consensus(trees, p = 1),
       majority = ape::consensus(trees, p = 0.5),
       n_trees  = length(trees))
}

## ---- 6. Implied weights, exact (enumeration) ----------------------------------------------------
## data must be prepared the way TreeSearch expects (e.g. with your clean_phydat_attributes()).
## k = Inf gives equal weights with the same engine, which is the right control.
## delta > 0 keeps slightly suboptimal trees (score within delta of the best).
## Check ?TreeSearch::TreeLength: here a LOWER score is assumed to be better.
iw_enumerate <- function(data, k, 
                         trees = TreeTools::RootTree(all_topologies(data), 1), 
                         delta = 0) {
  sc <- TreeSearch::TreeLength(trees, dataset = data, concavity = k)
  out <- trees[sc <= min(sc) + delta + 1e-9]
  attr(out, "score") <- min(sc)
  out
}

## Sweep over concavity constants, all on the same enumerated topologies
mpt_sweep <- function(data, 
                      trees=TreeTools::RootTree(all_topologies(data), 1),
                      k_grid = c(2, 3, 5, 10, 20, Inf), 
                      delta = 0) {
  res <- lapply(k_grid, function(k) iw_enumerate(data, k, trees, delta))
  setNames(res, paste0("k", k_grid))
}

## ---- 7. Implied weights, heuristic (for more than 9 taxa) ---------------------------------------
iw_search <- function(data, k, seed = 1, ratchIter = 30L, tbrIter = 6L, maxHits = 30L) {
  set.seed(seed)
  TreeSearch::MaximizeParsimony(data, concavity = k, ratchIter = ratchIter,
                                tbrIter = tbrIter, maxHits = maxHits, verbosity = 0L)
}

## ---- Example ----------------------------------------------------------------------------------
## trees <- all_topologies(data)
## ew_enum <- mpt_enumerate(data, trees)               # equal weights, exact
## ew_bab  <- mpt_bab(data)                            # same set by branch and bound
## same_tree_sets(ew_enum, ew_bab)                     # should be TRUE
## ew_collapsed <- collapse_zero(ew_enum, data)        # remove unsupported resolution
## cs <- consensus_set(ew_collapsed)                   # cs$strict, cs$majority
## sweep <- mpt_sweep(data_iw)                         # data_iw: data prepared for TreeSearch
## sapply(sweep, length)                               # number of optimal trees per concavity
## lapply(sweep, consensus_set)                        # compare with accuracy / resolution plots