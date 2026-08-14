library(phangorn)
library(xml2)
library(tracerer)
library(ape)

generate_xml_from_template <- function(phy_data, 
                                       dates_vector,
                                       template_xml_path, 
                                       output_xml_path, 
                                       chain_length = 10000000) {
  
  # 1. Read Template XML First
  doc <- read_xml(template_xml_path)
  
  # 2. Parse Input Data
  mat <- as.character(phy_data)
  taxa <- gsub("[;| ]", "_", rownames(mat)) # Clean taxon names (e.g. "305")
  n_taxa <- length(taxa)
  
  # Determine max state number dynamically
  unique_states <- sort(unique(as.vector(mat)))
  unique_states <- unique_states[!unique_states %in% c("-", "?")]
  max_states <- length(unique_states)
  
  # 3. Match Dates from dates_vector
  # Convert calendar years to "Years Before Present" (where youngest = 0)
  calendar_dates <- sapply(taxa, function(tx) {
    if (tx %in% names(dates_vector)) as.numeric(dates_vector[[tx]]) else 0
  })
  
  max_year <- max(calendar_dates, na.rm = TRUE)
  tip_ages <- max_year - calendar_dates  # Youngest becomes 0.0!
  
  trait_pairs <- paste0(taxa, "=", tip_ages)
  trait_string <- paste(trait_pairs, collapse = ",")
  
  # Max age is now the oldest manuscript's age:
  max_age <- max(tip_ages, na.rm = TRUE)
  safe_origin <- max_age + 100.0

  # Update origin parameter value in XML
  origin_node <- xml_find_first(doc, "//parameter[@id='originFBD.t:dummyname']")
  if (!is.na(origin_node)) {
    xml_text(origin_node) <- as.character(safe_origin)
  }
  
  # 5. Update Chain Length
  mcmc_node <- xml_find_first(doc, "//run[@id='mcmc']")
  if (!is.na(mcmc_node)) {
    xml_attr(mcmc_node, "chainLength") <- as.character(chain_length)
  }
  # increase number of initialisation attempts
  xml_attr(mcmc_node, "numInitializationAttempts") <- "1000"
  
  # 6. Populate <data id="dummyname">
  data_node <- xml_find_first(doc, "//data[@id='dummyname']")
  xml_replace(data_node, xml_new_root("data", id = "dummyname", spec = "Alignment", dataType = "standard", name = "alignment"))
  data_node <- xml_find_first(doc, "//data[@id='dummyname']")
  
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
    seq_str <- paste(mat[i, ], collapse = "")
    xml_add_child(data_node, "sequence", 
                  id = paste0("seq_", taxa[i]), 
                  taxon = taxa[i], 
                  totalcount = as.character(max_states), 
                  value = seq_str)
  }
  
  # 7. Update Date Trait Set
  trait_node <- xml_find_first(doc, "//trait[@id='dateTrait.t:dummyname']")
  if (!is.na(trait_node)) {
    xml_attr(trait_node, "value") <- trait_string
  }
  
  # 8. Update LewisMK State Number Dynamically
  mk_node <- xml_find_first(doc, "//substModel[@id='LewisMK.s:dummyname']")
  if (!is.na(mk_node)) {
    xml_attr(mk_node, "stateNumber") <- as.character(max_states)
  }
  
  # 9. Write Executable XML
  write_xml(doc, output_xml_path)
  message("BEAST2 XML successfully adapted and saved to: ", output_xml_path)
}