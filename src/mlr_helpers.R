# extracts coefficients from model
extr_glm_coef <- function(l, ...) {
  l$model %>%
    coef(...) %>%
    map(~ .x[which(.x != 0), 1]) %>%
    map(~ if (is.null(names(.x))) set_names(.x, "") else .x) %>%
    imap_dfr(~ data.frame(class = .y, feature = names(.x), coef = .x, row.names = NULL))
}

coef_summary <- function(df) {
  n_iter <- max(df$iteration)
  df %>%
    filter(feature != "(Intercept)") %>%
    group_by(class, feature) %>%
    summarize(n = n(), median_coef = median(coef), .groups = "drop") %>%
    mutate(prop = n / n_iter)
}

multiclass_auc <- function(cl) {
  R6::R6Class(paste0("MeasureAUC.", cl),
    inherit = mlr3::MeasureClassif,
    public = list(
      initialize = function() {
        super$initialize(
          id = paste0("auc.", cl),
          packages = "mlr3measures",
          properties = character(),
          predict_type = "prob",
          range = c(0, 1),
          minimize = FALSE
        )
      }
    ),
    private = list(
      # custom scoring function operating on the prediction object
      .score = function(prediction, ...) {
        mlr3measures::auc(factor(prediction$truth == cl), prediction$prob[, cl], "TRUE")
      }
    )
  )
}

multiclass_mcr <- function(cl) {
  R6::R6Class(paste0("MeasureMCR", cl),
    inherit = mlr3::MeasureClassif,
    public = list(
      initialize = function() {
        super$initialize(
          id = paste0("mcr.", cl),
          packages = "mlr3measures",
          properties = character(),
          predict_type = "response",
          range = c(0, 1),
          minimize = FALSE
        )
      }
    ),
    private = list(
      # custom scoring function operating on the prediction object
      .score = function(prediction, ...) {
        1 - mlr3measures::acc(
          factor(prediction$truth == cl, levels = c("FALSE", "TRUE")),
          factor(prediction$response == cl, levels = c("FALSE", "TRUE"))
        )
      }
    )
  )
}


set_task_name <-  . %>%
  mutate(task_name = case_when(
    task_id == "all_rnaseq" ~ "All genes (RNA-Seq)",
    task_id == "immune_rnaseq" ~ "Only immune-specific genes (RNA-Seq)",
    task_id == "non_immune_rnaseq" ~ "Non immune-specific genes (RNA-Seq)",
    task_id == "sc_scores" ~ "Cell type signature scores",
    task_id == "cancer_scores" ~ "Cancer-specific signature scores",
    task_id == "hallmark_scores" ~ "Hallmark signature scores"
  )) %>%
  mutate(task_name = factor(task_name,
                            levels = c("All genes (RNA-Seq)",
                                       "Only immune-specific genes (RNA-Seq)",
                                       "Non immune-specific genes (RNA-Seq)",
                                       "Cell type signature scores",
                                       "Cancer-specific signature scores",
                                       "Hallmark signature scores")))
