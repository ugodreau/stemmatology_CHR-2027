generate_beast2_mk_xml <- function(phy_data, xml_outfile, chain_length = 1000000, sample_every = 1000) {
  
  # 1. Extract matrix and taxa names
  mat <- as.character(phy_data)
  taxa <- gsub("[;| ]", "_", rownames(mat))
  n_taxa <- length(taxa)
  n_sites <- ncol(mat)
  
  # Find unique numeric state levels (e.g., "0", "1", "2", ..., "6")
  unique_states <- sort(unique(as.vector(mat)))
  unique_states <- unique_states[!unique_states %in% c("-", "?")] # Exclude gaps
  k_states <- length(unique_states)
  state_code_str <- paste(unique_states, collapse = " ")
  
  # 2. Build Base XML Structure
  doc <- xml2::xml_new_root("beast", 
                      version = "2.7", 
                      namespace = "beast.pkgannotator:beast.base.inference:beast.base.evolution.alignment:beast.base.evolution.tree:beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel")
  
  # Add UserDataType definition (for discrete numeric states)
  user_data_type <- xml2::xml_add_child(doc, "userDataType", 
                                  id = "customDataType", 
                                  spec = "beast.base.evolution.datatype.UserDataType", 
                                  states = as.character(k_states), 
                                  codeMap = state_code_str)
  
  # Add Sequence Alignment
  alignment <- xml2::xml_add_child(doc, "alignment", id = "alignment", spec = "Alignment", userDataType = "@customDataType")
  for (i in seq_len(n_taxa)) {
    seq_str <- paste(mat[i, ], collapse = "")
    xml2::xml_add_child(alignment, "sequence", id = paste0("seq_", taxa[i]), taxon = taxa[i], totalcount = as.character(k_states), value = seq_str)
  }
  
  # 3. Add MCMC Block
  mcmc <- xml2::xml_add_child(doc, "run", id = "mcmc", spec = "MCMC", chainLength = as.character(chain_length))
  
  # State Parameters (Tree + Clock)
  state <- xml2::xml_add_child(mcmc, "state", id = "state", storeEvery = "5000")
  xml2::xml_add_child(state, "tree", id = "Tree.t:alignment", spec = "Tree", name = "stateNode")
  
  # Init Random Tree
  init <- xml2::xml_add_child(mcmc, "init", id = "RandomTree.t:alignment", spec = "RandomTree", estimate = "false", initial = "@Tree.t:alignment", taxa = "@alignment")
  xml2::xml_add_child(init, "populationModel", id = "ConstantPopulation0", spec = "ConstantPopulation")
  
  # Posterior Distribution
  posterior <- xml2::xml_add_child(mcmc, "distribution", id = "posterior", spec = "CompoundDistribution")
  prior <- xml2::xml_add_child(posterior, "distribution", id = "prior", spec = "CompoundDistribution")
  
  # Yule Tree Prior
  yule <- xml2::xml_add_child(prior, "distribution", id = "YuleModel.t:alignment", spec = "beast.base.evolution.speciation.YuleModel", tree = "@Tree.t:alignment")
  xml2::xml_add_child(yule, "parameter", id = "birthRate.t:alignment", spec = "parameter.RealParameter", name = "birthRate", value = "1.0")
  
  # Likelihood Calculation (Lewis Mk + Gamma)
  likelihood <- xml2::xml_add_child(posterior, "distribution", id = "likelihood", spec = "CompoundDistribution")
  tree_lik <- xml2::xml_add_child(likelihood, "distribution", id = "treeLikelihood", spec = "TreeLikelihood", data = "@alignment", tree = "@Tree.t:alignment")
  
  # Site Model (JC69 as Mk model + Gamma Rate Categories)
  site_model <- xml2::xml_add_child(tree_lik, "siteModel", id = "SiteModel", spec = "SiteModel", gammaCategoryCount = "4")
  xml2::xml_add_child(site_model, "substModel", id = "LewisMk", spec = "JukesCantor")
  
  # Operators (Tree Topology Proposals)
  xml2::xml_add_child(mcmc, "operator", id = "YuleBirthRateScaler", spec = "kernel.base.ScalerOperator", parameter = "@birthRate.t:alignment", weight = "3.0")
  xml2::xml_add_child(mcmc, "operator", id = "SubtreeSlide", spec = "operator.kernel.SubtreeSlide", tree = "@Tree.t:alignment", weight = "15.0")
  xml2::xml_add_child(mcmc, "operator", id = "NarrowExchange", spec = "operator.kernel.Exchange", tree = "@Tree.t:alignment", weight = "15.0")
  xml2::xml_add_child(mcmc, "operator", id = "WideExchange", spec = "operator.kernel.Exchange", isNarrow = "false", tree = "@Tree.t:alignment", weight = "3.0")
  xml2::xml_add_child(mcmc, "operator", id = "WilsonBalding", spec = "operator.kernel.WilsonBalding", tree = "@Tree.t:alignment", weight = "3.0")
  
  # Loggers
  logger_file <- xml2::xml_add_child(mcmc, "logger", id = "tracelog", fileName = paste0(xml_outfile, ".log"), logEvery = as.character(sample_every), spec = "Logger")
  xml2::xml_add_child(logger_file, "log", idref = "posterior")
  xml2::xml_add_child(logger_file, "log", idref = "likelihood")
  
  logger_tree <- xml2::xml_add_child(mcmc, "logger", id = "treelog", fileName = paste0(xml_outfile, ".trees"), logEvery = as.character(sample_every), spec = "Logger")
  xml2::xml_add_child(logger_tree, "log", idref = "Tree.t:alignment")
  
  # Write XML to file
  xml2::write_xml(doc, xml_outfile)
  message("BEAST2 XML successfully written to: ", xml_outfile)
}