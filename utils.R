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

convert_gt_internal_nodes <- function(file="") {
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
