pls_multi_validate_positive_integer <- function(value, name) {
if (length(value) != 1L || is.na(value)) {
pls_multi_stop(paste(name, "must be a single positive integer"))
}

value_int <- as.integer(value)
if (!identical(as.numeric(value_int), as.numeric(value)) || value_int < 1L) {
pls_multi_stop(paste(name, "must be a single positive integer"))
}

value_int
}

pls_multi_validate_limited_numeric <- function(value, name) {
if (length(value) != 1L || is.na(value) || !is.numeric(value) || !is.finite(value)) {
pls_multi_stop(paste(name, "must be a single finite numeric value"))
}

as.numeric(value)
}

pls_multi_reject_cv_unsupported <- function(modele, family, weights_supplied,
                                            keepMclassed, EstimXNA,
                                            pvals.expli, alpha.pvals.expli,
                                            MClassed, sparse,
                                            sparseStop, naive) {
if (!identical(modele, "pls")) {
pls_multi_stop("cv.plsRmulti currently supports only modele = \"pls\"")
}
if (!is.null(family)) {
pls_multi_stop("cv.plsRmulti does not support GLM families in this experimental release")
}
if (weights_supplied) {
pls_multi_stop("cv.plsRmulti does not support weights in this experimental release")
}
if (!identical(keepMclassed, FALSE) || !identical(MClassed, FALSE)) {
pls_multi_stop("cv.plsRmulti does not support classification diagnostics in this experimental release")
}
if (!identical(EstimXNA, FALSE)) {
pls_multi_stop("cv.plsRmulti does not support missing-data estimation in this experimental release")
}
if (!identical(pvals.expli, FALSE) || !identical(alpha.pvals.expli, .05)) {
pls_multi_stop("cv.plsRmulti does not support predictor p-values in this experimental release")
}
if (!identical(sparse, FALSE) || !identical(sparseStop, FALSE)) {
pls_multi_stop("cv.plsRmulti does not support sparse extraction in this experimental release")
}
if (!identical(naive, FALSE)) {
pls_multi_stop("cv.plsRmulti does not support naive degrees of freedom options in this experimental release")
}
}

pls_multi_wrap_partition <- function(grouplist, NK, K) {
if (is.null(grouplist)) {
return(NULL)
}

if (!is.list(grouplist)) {
pls_multi_stop("grouplist must be a list when supplied")
}

if (length(grouplist) == K && all(vapply(grouplist, function(x) !is.list(x), logical(1)))) {
return(list(grouplist))
}

if (length(grouplist) != NK) {
pls_multi_stop("grouplist must contain one partition per NK repetition")
}

grouplist
}

pls_multi_partition_from_grouplist <- function(partition, n, K) {
if (!is.list(partition) || length(partition) != K) {
pls_multi_stop("Each grouplist partition must be a list of K groups")
}

partition <- lapply(partition, function(group) {
group <- sort(unique(as.integer(group)))
group[!is.na(group)]
})

is_holdout_partition <- function(groups) {
covered <- unlist(groups, use.names = FALSE)
identical(sort(covered), seq_len(n))
}

if (is_holdout_partition(partition)) {
return(partition)
}

lapply(partition, function(group) {
setdiff(seq_len(n), group)
})
}

pls_multi_make_cv_groups <- function(n, K, NK, random, grouplist) {
wrapped_grouplist <- pls_multi_wrap_partition(grouplist, NK = NK, K = K)
groups_all <- vector("list", NK)

for (nnkk in seq_len(NK)) {
if (!is.null(wrapped_grouplist)) {
groups_all[[nnkk]] <- pls_multi_partition_from_grouplist(
partition = wrapped_grouplist[[nnkk]],
n = n,
K = K
)
next
}

indices <- if (random) sample(seq_len(n), replace = FALSE) else seq_len(n)
group_ids <- rep(seq_len(K), length.out = n)
groups_all[[nnkk]] <- unname(split(indices, group_ids))
}

groups_all
}

pls_multi_cv_fit_fold <- function(dataY, dataX, holdout_idx, nt,
                                  scaleX, scaleY, tol_Xi,
                                  keepcoeffs, verbose) {
n <- nrow(dataX)
all_idx <- seq_len(n)
if (length(holdout_idx) == n) {
train_idx <- all_idx
test_idx <- all_idx
} else {
train_idx <- setdiff(all_idx, holdout_idx)
test_idx <- holdout_idx
}

fold_fit <- plsRmulti(
object = dataY[train_idx, , drop = FALSE],
dataX = dataX[train_idx, , drop = FALSE],
nt = nt,
scaleX = scaleX,
scaleY = scaleY,
tol_Xi = tol_Xi,
verbose = verbose
)

scaled_test <- pls_multi_scale_newdata(dataX[test_idx, , drop = FALSE], fold_fit)
pred_scaled_path <- pls_multi_predict_response_path_scaled(
x_scaled = scaled_test,
object = fold_fit,
max_comps = fold_fit$computed_nt
)
pred_path <- pls_multi_backtransform_response_path(
y_scaled_path = pred_scaled_path,
object = fold_fit
)

list(
predictions = pred_path,
observed = dataY[test_idx, , drop = FALSE],
coeffs = if (keepcoeffs) pls_multi_standardized_boot_vector(fold_fit) else NULL,
computed_nt = fold_fit$computed_nt,
train_idx = train_idx
)
}

pls_multi_cv_min_components <- function(cv_object) {
min(vapply(cv_object$results_kfolds, function(one_nk) {
min(vapply(one_nk, function(one_fold) {
dim(one_fold)[3]
}, integer(1L)))
}, integer(1L)))
}

pls_multi_cv_press <- function(cv_object, max_comps) {
press_total <- vector("list", length(cv_object$results_kfolds))
press_by_response <- vector("list", length(cv_object$results_kfolds))

for (nnkk in seq_along(cv_object$results_kfolds)) {
press_total[[nnkk]] <- rep(0, max_comps)
press_by_response[[nnkk]] <- matrix(
0,
nrow = max_comps,
ncol = cv_object$ny,
dimnames = list(pls_multi_component_labels(max_comps), cv_object$response_names)
)

for (ii in seq_along(cv_object$results_kfolds[[nnkk]])) {
pred_path <- cv_object$results_kfolds[[nnkk]][[ii]][, , seq_len(max_comps), drop = FALSE]
observed <- cv_object$dataY_kfolds[[nnkk]][[ii]]
for (comp in seq_len(max_comps)) {
err_sq <- (pred_path[, , comp] - observed)^2
press_by_response[[nnkk]][comp, ] <- press_by_response[[nnkk]][comp, ] + colSums(err_sq)
press_total[[nnkk]][comp] <- press_total[[nnkk]][comp] + sum(err_sq)
}
}
}

list(total = press_total, by_response = press_by_response)
}

pls_multi_cv_summary_matrix <- function(press_total, press_by_response,
                                        reference_fit, limQ2set) {
computed_nt <- length(press_total)
limQ2 <- rep(limQ2set, computed_nt)
q2 <- 1 - press_total / reference_fit$RSS_total[-1]
q2cum <- rep(NA_real_, computed_nt)
for (comp in seq_len(computed_nt)) {
q2cum[comp] <- 1 - prod(press_total[seq_len(comp)]) / prod(reference_fit$RSS_total[2:(comp + 1L)])
}

summary_mat <- cbind(
AIC = rep(NA_real_, computed_nt + 1L),
Q2cum_Y = c(NA_real_, q2cum),
LimQ2_Y = c(NA_real_, limQ2),
Q2_Y = c(NA_real_, q2),
PRESS_Y = c(NA_real_, press_total),
RSS_Y = reference_fit$RSS_total[seq_len(computed_nt + 1L)],
R2_Y = c(NA_real_, reference_fit$R2_total[seq_len(computed_nt)]),
AIC.std = rep(NA_real_, computed_nt + 1L)
)

for (resp_idx in seq_along(reference_fit$response_names)) {
resp_name <- reference_fit$response_names[resp_idx]
    summary_mat <- cbind(
      summary_mat,
      stats::setNames(
        data.frame(
          PRESS = c(NA_real_, press_by_response[, resp_idx]),
          RSS = reference_fit$RSS_by_response[seq_len(computed_nt + 1L), resp_idx],
          Q2 = c(
            NA_real_,
            1 - press_by_response[, resp_idx] / reference_fit$RSS_by_response[-1, resp_idx]
          ),
          R2 = c(NA_real_, reference_fit$R2_by_response[seq_len(computed_nt), resp_idx]),
          check.names = FALSE
        ),
        c(
          paste("PRESS_", resp_name, sep = ""),
          paste("RSS_", resp_name, sep = ""),
          paste("Q2_", resp_name, sep = ""),
          paste("R2_", resp_name, sep = "")
        )
      )
    )
  }

summary_mat <- as.matrix(summary_mat)
rownames(summary_mat) <- paste("Nb_Comp_", 0:computed_nt, sep = "")
attr(summary_mat, "computed_nt") <- computed_nt
summary_mat
}
