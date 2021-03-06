% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/bestSubset.r
\name{squared}
\alias{squared}
\title{Squared likelihood function}
\usage{
squared(fitness_metrics)
}
\arguments{
\item{fitness_metrics}{These metrics will be provided automatically by the fitnessMetric parameter in bestSubset()}
}
\value{
The function returns a vector with the squared values of the given vector.
}
\description{
The squared function is for use as the "likelihoodFn" in bestSubset().  It will be greatest with larger input values and is useful for maximizing model descriptive statistics such as AUC, area under the ROC curve, or R^2. NOTE: minimize_fitness_metric in bestSubset() should be set to FALSE when using this likelihood function.
}
\examples{
See help for bestSubset()
}

