if (!requireNamespace("fastglm", quietly = TRUE)) {
  stop("The 'fastglm' package is required for this benchmark.", call. = FALSE)
}

library(plsRglm)

benchmark_reps <- as.integer(Sys.getenv("PLSRGLM_BENCH_REPS", "5"))

bench_case <- function(expr, reps = benchmark_reps) {
  elapsed <- replicate(reps, system.time(force(expr))[["elapsed"]])
  c(
    min = min(elapsed),
    median = stats::median(elapsed),
    mean = mean(elapsed),
    max = max(elapsed)
  )
}

data(Cornell, package = "plsRglm")
data(aze_compl, package = "plsRglm")

results <- rbind(
  cbind(
    case = "Cornell gaussian",
    backend = "stats",
    t(bench_case(plsRglm(
      Y ~ ., data = Cornell, nt = 2,
      modele = "pls-glm-gaussian",
      fit_backend = "stats",
      verbose = FALSE
    )))
  ),
  cbind(
    case = "Cornell gaussian",
    backend = "fastglm",
    t(bench_case(plsRglm(
      Y ~ ., data = Cornell, nt = 2,
      modele = "pls-glm-gaussian",
      fit_backend = "fastglm",
      verbose = FALSE
    )))
  ),
  cbind(
    case = "aze_compl logistic",
    backend = "stats",
    t(bench_case(plsRglm(
      y ~ ., data = aze_compl, nt = 2,
      modele = "pls-glm-logistic",
      fit_backend = "stats",
      verbose = FALSE
    )))
  ),
  cbind(
    case = "aze_compl logistic",
    backend = "fastglm",
    t(bench_case(plsRglm(
      y ~ ., data = aze_compl, nt = 2,
      modele = "pls-glm-logistic",
      fit_backend = "fastglm",
      verbose = FALSE
    )))
  )
)

print(as.data.frame(results), row.names = FALSE)
