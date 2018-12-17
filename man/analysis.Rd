\name{analysis}
\alias{analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to look for significant difference between steady states between two genotypes
}
\description{
This function takes the two given data frames of licor data from two different genotypes and runs ttests at each CO2 level 400ppm, 800ppm, 100ppm to look for significant differences between gsw and A values. This function also produces graphs showing averages and standard error bars of gsw and A for both genotypes. 
}
\usage{
analysis(genotype1, genotype2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotype1}{
A data frame output by \code{\link{cleandata}} containing measurements relevant for graphing CO2 steady states of one genotype replication. 
}
  \item{genotype2}{
A data frame output by \code{\link{cleandata}} containing measurements relevant for graphing CO2 steady states of one genotype replication. 
}
}

\value{
\item{ttest1}{A t-test comparing the last 10 genotype 1 gsw data points of the 400ppm CO2 stead-state and the last 10 gsw steady-state points at 400ppm from genotype 2.}  
\item{ttest2}{A t-test comparing the last 10 genotype 1 gsw data points of the 800ppm CO2 stead-state and the last 10 gsw steady-state points at 800ppm from genotype 2.}  
\item{ttest3}{A t-test comparing the last 10 genotype 1 gsw data points of the 100ppm CO2 stead-state and the last 10 gsw steady-state points at 100ppm from genotype 2.}  
\item{ttest4}{A t-test comparing the last 10 genotype 1 A data points of the 400ppm CO2 stead-state and the last 10 A steady-state points at 400ppm from genotype 2.}  
\item{ttest5}{A t-test comparing the last 10 genotype 1 A data points of the 800ppm CO2 stead-state and the last 10 A steady-state points at 800ppm from genotype 2.}  
\item{ttest6}{A t-test comparing the last 10 genotype 1 A data points of the 100ppm CO2 stead-state and the last 10 A steady-state points at 100ppm from genotype 2.}  
\item{a_graph}{ A bar graph with Net Photosynthesis on the y-axis and the two genotypes being shown at three different CO2 levels (400ppm, 800ppm, 100ppm) on the x-axis. Error bars indicate the standard error for each average.}
\item{gsw_graph}{ A bar graph with stomatal conductance on the y-axis and the two genotypes being shown at three different CO2 levels (400ppm, 800ppm, 100ppm) on the x-axis. Error bars indicate the standard error for each average.}

}

\author{
Robert Twohey III
}

\examples{
#import example data
data1 <- system.file("extdata", "func3_1_example.csv", package = "licordataanalysis")
data2 <- system.file("extdata", "func3_2_example.csv", package = "licordataanalysis")

#run function
data1 <- read.csv(data1)
data2 <- read.csv(data2)
analize <- analysis(data1, data2)
analize[[1]]
analize[[2]]
}

