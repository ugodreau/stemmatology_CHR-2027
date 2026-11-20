## Notes
## - Uses (exact) tip dates by default
## - NJ initialisation (TODO)
source("generate_beast.R")

fbd <- function(data, 
                dates_vector,
                template_file = "beast_templates/FBD_NJinit_serial.xml", 
                nj_init = TRUE,
                out = "test", 
                chain_length=1000000, 
                bin_folder = "~/.local/share/beast/bin/"){
  
  ## Clean potential files from previous run
  output_files <- c(
    paste0(out, ".log"),
    paste0(out, "-dummyname.trees"),
    paste0(out, ".xml.state")
  )
  file.remove(output_files[file.exists(output_files)])
  
  # and now start
  xml_out = paste(out, ".xml", sep = "")
  # Generate config file
  generate_xml_from_template(
    data,
    dates_vector,
    nj_init = nj_init,
    template_file,
    xml_out,
    chain_length = chain_length
  )
  # Run Beast2
  beast_bin <- paste(bin_folder, "beast", sep="")
  system2(
    command = beast_bin, 
    args = c("-overwrite -seed 42 -prefix", sub("[^/]+$", "", out), 
             xml_out), 
    stdout = paste(out, ".txt", sep=""), stderr = paste(out, ".txt", sep="")
    )
  
  if (!file.exists(output_files[2])) {
    warning("BEAST2 run failed for: ", out)
    return(NULL) # Return NULL so scoring functions know it failed
  }
  
  # Parse & Process Posterior Trees
  tree_file <- paste(out, "-dummyname.trees", sep="")
  
  BI = return_bi_trees(tree_file, bin_folder)
  
  return(BI)
}

return_bi_trees <- function(tree_file, bin_folder, out="temp.txt", strip_prefix = "ms_"){
  trees <- ape::read.nexus(tree_file)
  
  # Remove 10% burn-in manually
  n_trees <- length(trees)
  burnin_cutoff <- floor(0.1 * n_trees)
  trees_after_burnin <- trees[(burnin_cutoff + 1):n_trees]
  
  mrc_tree = consensus(trees_after_burnin, p = 0.5)
  
  # Handle Sampled Ancestors in Consensus Trees
  treeannotator_bin <- paste(bin_folder, "treeannotator", sep="")
  system2(command = treeannotator_bin, 
          args = c("-burnin 10 -height median", tree_file, "mcc_tree.tree"),
          stdout = out, stderr = out)
  # Read the MCC tree with annotations preserved
  # and clean labels…
  raw <- readLines("mcc_tree.tree", warn = FALSE)
  clean_raw <- gsub("ms\\_|[|]\\S*", "", raw)
  tmp <- tempfile(fileext = ".nex")
  writeLines(clean_raw, tmp)
  mcc_tree <- read.nexus(tmp)
  unlink(tmp)
  
  mrc_tree$tip.label = gsub("ms\\_|[|]\\S*", "", mrc_tree$tip.label)
  
  BI = list()
  BI$mcc_tree = mcc_tree
  BI$mrc_tree = mrc_tree
  return(BI)
}
