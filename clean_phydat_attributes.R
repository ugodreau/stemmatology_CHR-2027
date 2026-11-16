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