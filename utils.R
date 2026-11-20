library(ape)
library(stringr)
library(phytools)

### Cleaning or conversion functions
clean_phydat_attributes <- function(dataset) {
  # Strip NA entries from the hidden level attributes
  if (!is.null(attr(dataset, "levels"))) {
    attr(dataset, "levels") <- attr(dataset, "levels")[!is.na(attr(dataset, "levels"))]
  }
  if (!is.null(attr(dataset, "allLevels"))) {
    attr(dataset, "allLevels") <- attr(dataset, "allLevels")[!is.na(attr(dataset, "allLevels"))]
  }
  return(dataset)
}
## Read a ground-truth stemma in Newick with node metadata such as
##   ms_38[state=False:birth_time=133.07:length=2]
## and return a tip-labelled tree on the surviving witnesses (state=True).
##
## Rules
##  * Metadata brackets are read for the state flag, then removed.
##  * The label prefix (default "ms_") is stripped, so tip labels match the alignment ("38", not "ms_38").
##  * A surviving witness that is an ancestor of other witnesses (internal node with state=True)
##    is represented by an extra leaf carrying its label:
##      attach = "star":  the leaf is one more child of its own node  ->  (C1, C2, X)
##                        (no split is invented between the copies of X; recommended)
##      attach = "above": the leaf is the sister of the whole descendant bundle -> ((C1, C2), X)
##                        (what ape::bind.tip(where = node) does: it resolves the polytomy)
##  * Lost nodes (state=False) that keep a single surviving child are collapsed; lost leaves are dropped.
## The labels of converted ancestors are stored in attr(tree, "ancestors").
## Requires: ape.
convert_gt_internal_nodes <- function(file = "", 
                                      strip_prefix = "ms_", 
                                      attach = c("star", "above")) {
  attach <- match.arg(attach)
  txt <- paste(readLines(file, warn = FALSE), collapse = "")
  
  ## 1. Node states: "label[state=True...]" or "label[state=False...]"
  m <- regmatches(txt, gregexpr("[0-9A-Za-z_]+\\[state=(?:True|False)[^\\]]*\\]", txt, perl = TRUE))[[1]]
  if (length(m) == 0) stop("No [state=...] metadata found in ", file)
  lab   <- gsub(strip_prefix, "", sub("\\[.*$", "", m), fixed = TRUE)
  state <- setNames(grepl("state=True", m, fixed = TRUE), lab)
  
  ## 2. Plain Newick: drop every bracketed metadata block and the label prefix
  plain <- gsub("\\[[^\\]]*\\]", "", txt, perl = TRUE)
  plain <- gsub(strip_prefix, "", plain, fixed = TRUE)
  phy   <- ape::read.tree(text = plain)
  
  ## 3. Rebuild the topology from the edge table
  nt     <- ape::Ntip(phy)
  kids   <- split(phy$edge[, 2], phy$edge[, 1])
  labels <- c(phy$tip.label, phy$node.label)
  anc    <- character(0)
  
  rec <- function(n) {
    me <- labels[n]
    ch <- kids[[as.character(n)]]
    if (is.null(ch)) return(if (isFALSE(unname(state[me]))) NULL else me)   # leaf
    parts <- Filter(Negate(is.null), lapply(ch, rec))
    if (isTRUE(unname(state[me]))) {                                         # surviving ancestor
      anc <<- c(anc, me)
      if (length(parts) == 0) return(me)
      if (attach == "star") parts <- c(parts, me)
      else parts <- list(paste0("(", paste(unlist(parts), collapse = ","), ")"), me)
    }
    if (length(parts) == 0) return(NULL)
    if (length(parts) == 1) return(parts[[1]])                               # collapse single-child nodes
    paste0("(", paste(unlist(parts), collapse = ","), ")")
  }
  
  nwk <- rec(nt + 1)
  if (is.null(nwk)) stop("No surviving witness in ", file)
  out <- ape::read.tree(text = paste0(nwk, ";"))
  
  ## 4. Sanity check: tips must be exactly the surviving witnesses
  expected <- names(state)[state]
  if (!setequal(out$tip.label, expected)) {
    warning("Tip labels do not match the state=True nodes in ", file)
  }
  attr(out, "ancestors") <- unique(anc)
  out
}

## ---- Example -----------------------------------------------------------------------------------
## gt <- convert_gt_internal_nodes("simulations_data/.../tree_0.newick")          # star attachment
## gt_above <- convert_gt_internal_nodes("simulations_data/.../tree_0.newick", attach = "above")
## phangorn::RF.dist(gt, gt_above)       # 0 if the two conventions give the same topology
## length(attr(gt, "ancestors"))         # number of surviving witnesses that are ancestors
## setdiff(gt$tip.label, names(dates_vector))   # should be character(0)


convert_gt_internal_nodes_old <- function(file="") {
  # 0. Read file
  newick_str = readLines(con = file, warn = FALSE)
  
  # 1. Strip branch lengths inside [] metadata tags so ape can parse it temporarily
  # Keeps [state=True] or [state=False] clean
  clean_str <- str_replace_all(newick_str, ":length=[0-9]+", "")
  
  # 2. Extract node states before stripping labels
  # Pattern to find labels and their state flags
  nodes_to_convert <- str_match_all(clean_str, "([0-9a-zA-Z_]+)\\[state=True\\]")[[1]][,2]
  
  # 3. Strip all NHX metadata brackets to leave a standard Newick string ape can read
  basic_newick <- str_replace_all(clean_str, "\\[state=(True|False)\\]", "")
  
  # Parse tree into ape
  phy <- read.tree(text = basic_newick)
  
  # 4. For internal nodes marked state=True, attach a 0-length tip in ape safely
  if (!is.null(phy$node.label)) {
    for (node_name in nodes_to_convert) {
      # Check if this state=True label sits on an internal node
      node_idx <- match(node_name, phy$node.label)
      
      if (!is.na(node_idx)) {
        # Convert internal node index to ape's internal indexing scheme
        internal_node_id <- Ntip(phy) + node_idx
        
        # Bind a new tip with 0 length safely using ape
        phy <- bind.tip(phy, tip.label = node_name, where = internal_node_id, edge.length = 0)
        
        # Clear label from the internal parent node to avoid duplication
        phy$node.label[node_idx] <- ""
      }
    }
  }
  
  return(phy)
}

### Loading functions

load_simu_data <-  function(file="", return_dates=FALSE){
  out = list()
  # Load data and gt tree
  ## load, clean, write, read each nexus file
  ### Data first
  raw <- readLines(file, warn = FALSE)
  if(return_dates){
    # create the dates vector for bayesian methods
    matches <- regmatches(raw, regexec("ms_([0-9]+)\\|([0-9]+)", raw))
    matched_matrix <- do.call(rbind, lapply(matches, function(x) if(length(x) == 3) x[2:3] else NULL))
    # Create a lookup vector: c("305" = 246, "108" = 120, ...)
    out$dates_vector <- setNames(as.numeric(matched_matrix[,2]), matched_matrix[,1])
  }
  # and now clean the names to avoid issues with nexus reading
  clean_raw <- gsub("ms\\_|[|]\\S*", "", raw)
  tmp <- tempfile(fileext = ".nex")
  writeLines(clean_raw, tmp)
  data <- ape::read.nexus.data(tmp)
  unlink(tmp)
  m <- do.call(rbind, data)
  max_state <- max(as.integer(m), na.rm = TRUE)
  out$data = phyDat(data, type = "USER", levels = 0:max_state, ambiguity=NA)
  return(out)
}


import.TEIApparatus <- function (file = "", appTypes = c("substantive")) 
{
  doc = xml2::read_xml(file)
  xml2::xml_ns_strip(doc)
  myWits = xml2::xml_text(xml2::xml_find_all(doc, "/descendant::listWit/descendant::witness/@xml:id"))
  myWits = sort(myWits)
  tableVariantes = matrix(nrow = 0, ncol = length(myWits), 
                          dimnames = list(NULL, myWits))
  if (is.null(appTypes)) {
    myApps = xml2::xml_find_all(doc, "/descendant::app")
  }
  else {
    myTypes = paste("@type = '", appTypes, "'", sep = "", 
                    collapse = " or ")
    myXpath = paste("/descendant::app[", myTypes, "]", sep = "")
    myApps = xml2::xml_find_all(doc, xpath = myXpath)
  }
  for (i in seq_len(length(myApps))) {
    attempt <- tryCatch({
      myRdgs = xml2::xml_find_all(myApps[i], "./rdg")
      myRdgWits = xml2::xml_text(xml2::xml_find_all(myApps[i], 
                                                    "./rdg/@wit"))
      myRdgWits = gsub("#", " ", myRdgWits)
      myRdgWits = sub("^\\s+|\\s+$", "", myRdgWits)
      myRdgWits = strsplit(myRdgWits, "\\s+")
      myRdgVars = NULL
      for (j in seq_len(length(myRdgs))) {
        if (xml2::xml_text(myRdgs[j], "normalize-space(.)") == 
            "") {
          myRdgVars = c(myRdgVars, 0)
        }
        else {
          myRdgVars = c(myRdgVars, j)
        }
      }
      myRow = NULL
      for (j in seq_len(length(myWits))) {
        myWitReading = myRdgVars[which(sapply(myRdgWits, 
                                              FUN = function(X) myWits[j] %in% X))]
        if (length(myWitReading) == 0) {
          myWitReading = NA
        }
        if (length(myWitReading) > 1) {
          myWitReading = paste(as.character(myWitReading), 
                               collapse = ",")
        }
        myRow = c(myRow, myWitReading)
      }
      myLabel = xml2::xml_text(xml2::xml_find_all(myApps[i], 
                                                  "@xml:id"))
      if (length(myLabel) == 0) {
        myLabel = as.character(i)
      }
      tableVariantes = rbind(tableVariantes, myRow)
      rownames(tableVariantes)[i] = myLabel
    }, error = function(e) {
      # Print a warning showing which file failed
      warning(paste("Importing failed for app", 
                    xml2::xml_text(xml2::xml_find_all(myApps[i], "./@xml:id"))))
      return(NULL)
    })
    # --- 2. Skip to next iteration if loading failed ---
    if (is.null(attempt)) {
      next
    }
  }
  return(tableVariantes)
}
