% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geomaps.R
\name{gps2kml}
\alias{gps2kml}
\title{Convert GPS coordinates}
\usage{
gps2kml(x)
}
\arguments{
\item{x}{data frame with point locations in rows and at least columns named "name", "longitude" and "latitude"; if columns "description" and "altitude" are present, they will also be used}
}
\value{
vector of characters (one item per line of the KML file)
}
\description{
Convert GPS coordinates into the \href{https://en.wikipedia.org/wiki/Keyhole_Markup_Language}{KML} format.
}
\examples{
coords <- data.frame(name=c("Paris", "Montpellier", "Toulouse"),
                     longitude=c(2.352222, 3.876716, 1.444209),
                     latitude=c(48.856614, 43.610769, 43.604652))
kml <- gps2kml(x=coords)

\dontrun{## save to a file
out.file <- "coords.kml"
cat(kml, file=out.file, sep="\n")
}
}
\author{
Timothee Flutre
}
