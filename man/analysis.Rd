\name{analysis}
\alias{analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to look for significant difference between steady states between two genotypes
}
\description{
This function takes the two given data frames of licor data from two different genotypes and runs ttests at each CO2 level (400ppm, 800ppm, 100ppm) to look for significant differences between gsw and A values. This function also produces graphs showing averages and standard error bars of gsw and A for both genotypes. 
}
\usage{
analysis(genotype1, genotype2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotype1}{
A data frame of licor data produced using \code{cleandata}
}
  \item{genotype2}{
A data frame of licor data produced using \code{cleandata}
}
}

\value{
\code{analysis} returns six sets of \code{t.test} results and two bar graphs showing averages and standard error.
}

\author{
Robert Twohey III
}

\examples{
#import example data
data <- system.file("extdata", "func23_example.csv", package = "licordataanalysis")

#run function
analize <- analysis(data)
}

