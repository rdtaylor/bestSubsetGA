#' Best modeling subset
#'
#' This function attempts to optimize a subset of variables for modeling based on given criteria with a genetic algorithm.  The use case is intended for situations where stepwise variable selection is either not likely to result in the best model or is impractical.
#'
#' @param dv A vector or matrix of dependent variable(s) for the model.
#' @param ivs A matrix containing the indepdendent variables on which the model is to be based.
#' @param buildModel A function which actually creates the model such as lm().
#' @param fitnessMetric A function taking the model built by buildModel and returns a numeric value specifying a model's fitness -- e.g. the model's variance, R^2, or AUC.  In Darwinian terms, metric will be used to determine the "fitness" of the model to survive and procreate.
#' @param likelihoodFn A function which takes the fitness metrics created by fitnessMetric and provides the probability of that fitness metric being selected for the next generation.  The two likelihood functions invLog() and squared() are examples provided to choose smaller and larger values, respectively.
#' @param minimize_fitness_metric This parameter specifies whether the fitness metric is to be minimized (TRUE) or maximized (FALSE).
#' @param n_generations The maximum number of generations to run the simulation.
#' @param n_genes The number of independent "organisms" to create in every generation for crossover or breeding.
#' @param mutation_rate The percentage chance for every variable to be flipped after crossover.  Higher values give higher volatility between generations and vice versa.
#' @param elitism_rate The rate at which best performers are to be kept, in order to prevent regression between generations.  The top elitism_rate percent will be preserved unchanged.
#' @param VERBOSE This parameter will provide output as to the progress.  The default, FALSE, will output none.
#' @param ... Any parameters needed to be passed to the given buildModel function can be added here.
#'
#' @return The function returns a vector containing the columns of the independent variables to use in the best subset found.
#'
#' @examples 
#' modelSigma <- function(fitted_model) {
#'     summary(fitted_model)$sigma
#' }
#' 
#' # Minimize the error variance of a linear model for given data
#' min_subset <- bestSubset(my_response, my_data, buildModel=lm,
#'                          fitnessMetric=modelSigma, likelihoodFn=invLog,
#'                          n_generations=250, n_genes=50, mutation_rate=0.005)
#' 
#' # Take only the columns specified by bestSubset()
#' selected_data <- my_data[, min_subset]
#' best_model <- lm(dv ~ ., data=selected_data)
#'
#' @export
bestSubset <- function(dv, ivs, buildModel, fitnessMetric, likelihoodFn,
                       minimize_fitness_metric=T, n_generations=100,
                       n_genes=50, mutation_rate=0.01, elitism_rate=0.01,
                       VERBOSE=FALSE, ...) {

  init_generation <- generateInitGeneration(buildModel, fitnessMetric, n_genes,
                                            dv, ivs, ...)
  genes <- init_generation[['genes']]
  fitness_metrics <- init_generation[['fitness_metrics']]

  # TODO -- check gradients to stop automatically
  for (i in seq(1, n_generations)) {
    genes <- selectGenes(genes, likelihoodFn, fitness_metrics,
                         minimize_fitness_metric, mutation_rate,
                         elitism_rate)
    fitness_metrics <- geneMetrics(buildModel, fitnessMetric, genes, dv, ivs, ...)

    if (VERBOSE) print(paste0("Generation ", i, " done."))
  }

  best_subset_index <- if (minimize_fitness_metric) which.min(fitness_metrics)
                       else which.max(fitness_metrics)
  best_subset_cols <- which(genes[best_subset_index, ] == 1)

  best_subset_cols
}

generateInitGeneration <- function(buildModel, fitnessMetric, n_genes, dv, ivs, ...) {

  generateGene <- function() {
    rbinom(p, 1, 0.5)
  }

  p <- ncol(ivs)
  genes <- generateGene()

  for (i in seq(2, n_genes)) {
    genes <- rbind(genes, generateGene())
  }

  fitness_metrics <- geneMetrics(buildModel, fitnessMetric, genes, dv, ivs, ...)

  list(genes=genes, fitness_metrics=fitness_metrics)
}

geneMetrics <- function(buildModel, fitnessMetric, genes, dv, ivs, ...) {
  gene_row_iter <- seq(1, nrow(genes))

  foreach(gene_row_index=gene_row_iter,
          .combine=c, .inorder=TRUE) %dopar% {
    used_variables <- genes[gene_row_index, ]
    fitness(dv, ivs, buildModel, fitnessMetric, used_variables, ...)
  }
}

selectGenes <- function(genes, likelihoodFn, fitness_metrics,
                        minimize_fitness_metric, mutation_rate,
                        elitism_rate) {

  n <- nrow(genes)
  p <- ncol(genes)

  n_elites <- ceiling(n * elitism_rate)

  new_genes <- if (n_elites > 0) {
                   order_decreasing <- if (minimize_fitness_metric) FALSE
                                       else TRUE
                   elite_rows <- head(order(fitness_metrics,
                                            decreasing=order_decreasing),
                                      n=n_elites)
                   genes[elite_rows, ]
                 }
                 else NULL

  selection_probabilities <- selectionProbabilities(likelihoodFn,
                                                    fitness_metrics)

  for (i in seq(n_elites + 1, n)) {
    parent_indexes <- sample(seq(1, n), 2, replace=T,
                             prob=selection_probabilities)
    offspring <- crossover(genes[parent_indexes[1], ],
                           genes[parent_indexes[2], ],
                           mutation_rate)
    new_genes <- if (is.null(new_genes)) offspring
                 else rbind(new_genes, offspring)
  }

  new_genes
}

selectionProbabilities <- function(likelihoodFn, fitness_metrics) {
  likelihoods <- likelihoodFn(fitness_metrics)
  total_likelihoods <- sum(likelihoods)
  selection_probabilities <- likelihoods / total_likelihoods

  selection_probabilities
}

crossover <- function(father, mother, mutation_rate) {
  p <- length(father)
  parents <- list(father, mother)

  offspring <- sapply(
    seq(1, p),
    function(i) {
      donor_selection <- rbinom(1, 1, 0.5) + 1
      parents[[donor_selection]][i]
    }
  )

  mutation_vector <- rbinom(p, 1, mutation_rate)
  offspring <- (offspring + mutation_vector) %% 2

  offspring
}

fitness <- function(dv, ivs, buildModel, fitnessMetric, used_variables, ...) {

  used_ivs <- subset(ivs, select=which(used_variables == 1))

  no_chosen_variables <- sum(used_variables) == 0
  fitted_model <-
    if (no_chosen_variables) {
      buildModel(dv ~ 1, ...)
    }
    else {
      buildModel(dv ~ ., data=used_ivs, ...)
    }

  fitnessMetric(fitted_model)
}

#' Inverse log likelihood function
#'
#' The inverse log function is for use as the "likelihoodFn" in bestSubset().  It will be greatest with smaller input values and is useful for minimizing model descriptive statistics such as error variance or prediction error. NOTE: minimize_fitness_metric in bestSubset() should be set to TRUE (the default) when using this likelihood function.
#'
#' @param fitness_metrics These metrics will be provided automatically by the fitnessMetric parameter in bestSubset()
#'
#' @return The function returns a vector with the inverse log of the given vector.
#'
#' @examples 
#' See help for bestSubset()
#'
#' @export
invLog <- function(fitness_metrics) {
  1 / log(fitness_metrics)
}

#' Squared likelihood function
#'
#' The squared function is for use as the "likelihoodFn" in bestSubset().  It will be greatest with larger input values and is useful for maximizing model descriptive statistics such as AUC, area under the ROC curve, or R^2. NOTE: minimize_fitness_metric in bestSubset() should be set to FALSE when using this likelihood function.
#'
#' @param fitness_metrics These metrics will be provided automatically by the fitnessMetric parameter in bestSubset()
#'
#' @return The function returns a vector with the squared values of the given vector.
#'
#' @examples 
#' See help for bestSubset()
#'
#' @export
squared <- function(fitness_metrics) {
  fitness_metrics ^ 2
}
