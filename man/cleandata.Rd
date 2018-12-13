\name{cleandata}
\alias{cleandata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for filtering and trimming raw Licor data. 
}
\description{
This function will take a raw Licor data csv file and remove extra information in the file above the logged data, remove unnecessary columns for other functions in the package, and remove any logged data that is an error due to machine matching.  
}
\usage{
cleandata(data, csvname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A numeric data frame with logs in rows and variables in columns
}
  \item{csvname}{
A file name indicating what the csv file produced by the function will be called
}
}
\value{
\code{cleandata} returns a data frame of numeric values.
}
\author{
Robert Twohey III
}
\examples{
#import example data
data2 <- system.file("extdata", "func1_example.csv", package = "licordataanalysis")

#run function
cleaned <- cleandata(data2)
}

