\name{graphPack}
\alias{graphPack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for producing stepwise graphs of stomatal conductance (gsw) and net photosynthesis (A) from cleaned Licor data
}
\description{
This function graphs the values of gsw and A at varrying CO2 concentrations (400ppm, 800 ppm, 100 ppm).
}
\usage{
graphPack(cleandata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cleandata}{
A data frame produced by \code{cleandata} containing gsw and A data to graph. 
}
}

\value{
\code{graphPack} returns two graphs.
}
\author{
Robert Twohey III
}

\examples{
#import example data
data1 <- system.file("extdata", "func23_example.csv", package = "licordataanalysis")

#run function
graphs <- graphPack(data1)
}
