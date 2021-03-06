% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assign_arguments.R
\name{assign_arguments}
\alias{assign_arguments}
\title{Assigns input arguments}
\usage{
assign_arguments(arguments, options)
}
\arguments{
\item{arguments}{character vector, with inputted arguments names.}

\item{options}{character vector, with all possible options.}
}
\value{
invisible NULL
}
\description{
From a list of available options, assigns all given arguments as TRUE
and the rest as FALSE, and gives error if any of the inputted arguments
doesn't match the possible ones.
}
