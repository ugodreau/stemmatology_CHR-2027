library(phangorn)
library(xml2)
library(tracerer)
library(ape)

get_node <- function(doc, xpath) {
  node <- xml_find_first(doc, xpath)
  if (inherits(node, "xml_missing")) stop("Template is missing node: ", xpath, call. = FALSE)
  node
}

set_param <- function(doc, id, value) {
  node <- get_node(doc, sprintf("//parameter[@id='%s']", id))
  xml_text(node) <- as.character(value)
}

generate_xml_from_template <- function(phy_data, 
                                       dates_vector,
                                       nj_init = TRUE,
                                       template_xml_path, 
                                       output_xml_path, 
                                       chain_length = 10000000) {
  
  # 1. Read Template XML First
  doc <- read_xml(template_xml_path)
  
  # 2. Parse Input Data
  mat <- as.character(phy_data)
  mat[is.na(mat)] <- "?" # BEAST's missing-data symbol
  
  taxa <- gsub("[;| ]", "_", rownames(mat)) # Clean taxon names (e.g. "305")
  n_taxa <- length(taxa)
  # Determine max state number dynamically
  #unique_states <- sort(unique(as.vector(mat)))
  #unique_states <- unique_states[!unique_states %in% c("-", "?")]
  unique_states <- setdiff(unique(as.vector(mat)), c("-", "?"))
  max_states <- max(as.integer(unique_states)) + 1 # some states can
  # be missing, due to lost manuscripts, so we avoid length(unique_states)
  sep <- if (max(nchar(unique_states)) > 1) "," else ""   # commas only when needed
  
  
  # 3. Match Dates from dates_vector
  # Convert calendar years to "Years Before Present" (where youngest = 0)
  # and check that none is missing
  missing <- setdiff(taxa, names(dates_vector))
  if (length(missing) > 0) stop("No date for taxa: ", paste(missing, collapse = ", "))
  calendar_dates <- as.numeric(dates_vector[taxa])
  #calendar_dates <- sapply(taxa, function(tx) {
  #  if (tx %in% names(dates_vector)) as.numeric(dates_vector[[tx]]) else 0
  #})
  
  max_year <- max(calendar_dates, na.rm = TRUE)
  tip_ages <- max_year - calendar_dates  # Youngest becomes 0.0!
  
  trait_pairs <- paste0(taxa, "=", tip_ages)
  trait_string <- paste(trait_pairs, collapse = ",")
  
  # Max age is now the oldest manuscript's age:
  max_age <- max(tip_ages, na.rm = TRUE)
  # Set the safe origin at max age + 125, 
  # and bound it between max_age and max_age + 250
  origin_slack <- 250   #TODO: add in paper
  origin_node  <- get_node(doc, "//parameter[@id='originFBD.t:dummyname']")
  xml_text(origin_node) <- as.character(max_age + origin_slack / 2)
  
  origin_prior <- get_node(doc, "//prior[@id='originPriorFBD.t:dummyname']/Uniform")
  xml_attr(origin_prior, "lower") <- as.character(max_age)
  xml_attr(origin_prior, "upper") <- as.character(max_age + origin_slack)

  # 5. Update Chain Length
  mcmc_node <- get_node(doc, "//run[@id='mcmc']")
  xml_attr(mcmc_node, "chainLength") <- format(chain_length, scientific = FALSE)
  # increase number of initialisation attempts
  xml_attr(mcmc_node, "numInitializationAttempts") <- "1000"
  
  # 6. Populate <data id="dummyname">
  data_node <- get_node(doc, "//data[@id='dummyname']")
  xml_replace(data_node, xml_new_root("data", id = "dummyname", spec = "Alignment", dataType = "standard", name = "alignment"))
  data_node <- get_node(doc, "//data[@id='dummyname']")
  
  # Add userDataType mapping
  xml_add_child(
    data_node, 
    "userDataType", 
    id = "StandardData.0", 
    spec = "beast.base.evolution.datatype.StandardData", 
    nrOfStates = as.character(max_states)
  )
  
  # Append Sequences
  for (i in seq_len(n_taxa)) {
    seq_str <- paste(mat[i, ], collapse = sep)
    xml_add_child(data_node, "sequence", 
                  id = paste0("seq_", taxa[i]), 
                  taxon = taxa[i], 
                  totalcount = as.character(max_states), 
                  value = seq_str)
  }
  
  # 7. Update Date Trait Set
  trait_node <- get_node(doc, "//trait[@id='dateTrait.t:dummyname']")
  xml_attr(trait_node, "value") <- trait_string
  xml_attr(trait_node, "traitname") <- "date-backward"
  
  # 8. Update LewisMK State Number Dynamically
  mk_node <- get_node(doc, "//substModel[@id='LewisMK.s:dummyname']")
  xml_attr(mk_node, "stateNumber") <- as.character(max_states)
  
  # 9. Keep or remove NJ initialisation
  #if (!nj_init) xml_remove(get_node(doc, "//init[@id='NJTree.t:dummyname']"))
  
  # 10. Set a reasonable initial value to some parameters
  #set_param(doc, "clockRate.c:dummyname",                 0.002)  # ~0.5 changes/site over ~250 years
  set_param(doc, "diversificationRateFBD.t:dummyname",    0.02)   # lambda - mu from Camps et al.
  set_param(doc, "turnoverFBD.t:dummyname",               0.45)   # mu / lambda
  set_param(doc, "samplingProportionFBD.t:dummyname",     0.05)
  set_param(doc, "originFBD.t:dummyname",                 max_age + 50)
  
  # 11. Write Executable XML
  write_xml(doc, output_xml_path)
  message("BEAST2 XML successfully adapted and saved to: ", output_xml_path)
}

## Diagnostics
#div <- rowSums(as.character(data) != "0")      # changes from the root state, per taxon
#cor(div, calendar_dates)           # tip dates in the original calendar

ll_on_topology <- function(tree, d_k) {
  tr <- ape::compute.brlen(ape::multi2di(ape::unroot(tree)), 0.05)
  optim.pml(pml(tr, d_k), optEdge = TRUE, control = pml.control(trace = 0))$logLik
}

bi_diagnostics <- function(i, burnin = 0.1, d_k) {
  # i: the index of the bayesian run
  # d_k the data loaded in the same fashion as used by BI
  # Compute ML to compare with BI
  ml_ll <- optim.pml(
    pml(NJ(dist.hamming(d_k)), d_k), 
    optNni = TRUE, optEdge = TRUE, 
    control = pml.control(trace = 0))$logLik
  
  lg <- read.table(sprintf("rbeast_out/sequences_%d.log", i), header = TRUE, comment.char = "#")
  lg <- lg[-seq_len(floor(burnin * nrow(lg))), ]
  
  L_topo         <- ll_on_topology(
    ape::read.tree(sprintf("rbeast_out/mcc_%d.tree", i)), d_k)
  mean_lik <- mean(lg$likelihood)
  
  c(
    tradition  = i, 
    mean_lik   = mean(lg$likelihood),
    ess_post   = coda::effectiveSize(lg$posterior),
    ess_lik    = coda::effectiveSize(lg$likelihood),
    ess_height = coda::effectiveSize(lg$Tree.height),
    sa_count   = mean(lg$SACountFBD),
    ml_diff    = ml_ll - mean_lik,
    topology_cost  = ml_ll - L_topo,
    branchlen_cost = L_topo - mean_lik
    )
}
