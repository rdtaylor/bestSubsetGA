#!/usr/bin/Rscript

library(testthat, quietly=T)
library(bestSubsetGA, quietly=T)
library(doMC, quietly=T)
registerDoMC(cores=2)


modelSigma <- function(fitted_model) {
  summary(fitted_model)$sigma
}

generateIVs <- function(n) {
  data.frame( A=seq(1, n)
            , B=runif(n, 0, 10)
            , C=runif(n, 0, 10)
            , D=runif(n, -100, 100)
            , E=runif(n, 0, 10)
            , F=rnorm(n, 0, 1000)
            , H=runif(n, 0, 10)
            , I=runif(n, -10, 10)
            , J=rnorm(n, 0, 10)
            , K=runif(n, 0, 10)
            )
}

generateDV <- function(ivs) {
  29 + (3 * ivs$A) - (2 * ivs$C) + (3.5 * ivs$D) + rnorm(nrow(ivs), 0, 1)
}


set.seed(1)

ivs <- generateIVs(1000)
dv <- generateDV(ivs)

test_that("Minimizes linear model's random error", {
  min_subset <- bestSubset(dv, ivs, buildModel=lm, fitnessMetric=modelSigma,
                           likelihoodFn=invLog, n_generations=250, n_genes=50,
                           mutation_rate=0.005)
  expect_equal(min_subset, c(1, 2, 3, 4, 5, 8, 9))
})

test_that("Maximizes linear model's random error", {
  max_subset <- bestSubset(dv, ivs, buildModel=lm, fitnessMetric=modelSigma,
                           likelihoodFn=squared, minimize_fitness_metric=FALSE,
                           n_generations=250, n_genes=50, mutation_rate=0.005)
  expect_equal(max_subset, c(2, 3, 5, 6, 7, 8, 9, 10))
})

print("All tests passed.")
