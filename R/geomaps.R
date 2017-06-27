## Contains functions useful for geographical maps.

##' Convert GPS coordinates
##'
##' Convert GPS coordinates into the \href{https://en.wikipedia.org/wiki/Keyhole_Markup_Language}{KML} format.
##' @param x data frame with point locations in rows and at least columns named "name", "longitude" and "latitude"; if columns "description" and "altitude" are present, they will also be used
##' @return vector of characters (one item per line of the KML file)
##' @author Timothee Flutre
##' @examples
##' coords <- data.frame(name=c("Paris", "Montpellier", "Toulouse"),
##'                      longitude=c(2.352222, 3.876716, 1.444209),
##'                      latitude=c(48.856614, 43.610769, 43.604652))
##' kml <- gps2kml(x=coords)
##'
##' \dontrun{## save to a file
##' out.file <- "coords.kml"
##' cat(kml, file=out.file, sep="\n")
##' }
##' @export
gps2kml <- function(x){
  stopifnot(is.data.frame(x),
            ncol(x) >= 3,
            all(c("name", "longitude", "latitude") %in% colnames(x)))

  if(! "altitude" %in% colnames(x))
    x$altitude <- NA

  out <- c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
  out <- c(out, "<kml xmlns=\"http://www.opengis.net/kml/2.2\">")
  out <- c(out, "<Document>")

  for(i in 1:nrow(x)){
    out <- c(out, "<Placemark>")
    txt <- paste0("<name>", x$name[i], "</name>")
    out <- c(out, txt)
    if("description" %in% colnames(x)){
      txt <- paste0("<description>", x$description[i], "</description>")
      out <- c(out, txt)
    }
    out <- c(out, "<Point>")
    txt <- paste0("<coordinates>", x$longitude[i], ",", x$latitude[i], ",",
                  x$altitude[i], "</coordinates>")
    out <- c(out, txt)
    out <- c(out, "</Point>")
    out <- c(out, "</Placemark>")
  }

  out <- c(out, "</Document>")
  out <- c(out, "</kml>")

  return(out)
}
