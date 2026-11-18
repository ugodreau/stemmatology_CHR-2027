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
