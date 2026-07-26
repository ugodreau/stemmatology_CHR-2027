library(ape)
library(stringr)
library(phytools)

convert_gt_internal_nodes <- function(newick_str) {
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