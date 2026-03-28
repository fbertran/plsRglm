pls_multi_stop <- function(message) {
stop(message, call. = FALSE)
}

pls_multi_check_dots <- function(dots) {
if (!length(dots)) {
return(invisible(NULL))
}

dot_names <- names(dots)
if (is.null(dot_names) || any(dot_names == "")) {
pls_multi_stop("plsRmulti does not support unnamed extra arguments in this experimental release")
}

pls_multi_stop(
paste(
"plsRmulti does not support these extra arguments in this experimental release:",
paste(dot_names, collapse = ", ")
)
)
}

pls_multi_reject_unsupported <- function(modele, family, typeVC, weights_supplied,
                                         limQ2set, dataPredictY_supplied,
                                         EstimXNA, pvals.expli,
                                         alpha.pvals.expli, MClassed,
                                         sparse, sparseStop, naive) {
if (!identical(modele, "pls")) {
pls_multi_stop("plsRmulti currently supports only modele = \"pls\"")
}
if (!is.null(family)) {
pls_multi_stop("plsRmulti does not support GLM families in this experimental release")
}
if (!identical(typeVC, "none")) {
pls_multi_stop("plsRmulti does not support cross-validation arguments in this experimental release")
}
if (weights_supplied) {
pls_multi_stop("plsRmulti does not support weights in this experimental release")
}
if (!identical(limQ2set, .0975)) {
pls_multi_stop("plsRmulti does not support limQ2set in this experimental release")
}
if (dataPredictY_supplied) {
pls_multi_stop("plsRmulti does not support dataPredictY; fit the model first and use predict()")
}
if (!identical(EstimXNA, FALSE)) {
pls_multi_stop("plsRmulti does not support missing-data estimation in this experimental release")
}
if (!identical(pvals.expli, FALSE)) {
pls_multi_stop("plsRmulti does not support p-values for predictors in this experimental release")
}
if (!identical(alpha.pvals.expli, .05)) {
pls_multi_stop("plsRmulti does not support alpha.pvals.expli in this experimental release")
}
if (!identical(MClassed, FALSE)) {
pls_multi_stop("plsRmulti does not support classification diagnostics in this experimental release")
}
if (!identical(sparse, FALSE) || !identical(sparseStop, FALSE)) {
pls_multi_stop("plsRmulti does not support sparse extraction in this experimental release")
}
if (!identical(naive, FALSE)) {
pls_multi_stop("plsRmulti does not support naive degrees of freedom options in this experimental release")
}
}

pls_multi_validate_nt <- function(nt) {
if (length(nt) != 1L || is.na(nt)) {
pls_multi_stop("nt must be a single positive integer")
}
nt_int <- as.integer(nt)
if (!identical(as.numeric(nt_int), as.numeric(nt)) || nt_int < 1L) {
pls_multi_stop("nt must be a single positive integer")
}
nt_int
}

pls_multi_prepare_response <- function(dataY, tol_Xi) {
if (is.data.frame(dataY)) {
non_numeric <- !vapply(dataY, is.numeric, logical(1))
if (any(non_numeric)) {
pls_multi_stop("plsRmulti requires a numeric multivariate response")
}
dataY <- as.matrix(dataY)
} else {
dataY <- as.matrix(dataY)
}

if (!is.numeric(dataY)) {
pls_multi_stop("plsRmulti requires a numeric multivariate response")
}
if (is.null(dim(dataY)) || ncol(dataY) < 2L) {
pls_multi_stop("plsRmulti requires at least two response columns; use plsR() for univariate responses")
}
if (nrow(dataY) < 2L) {
pls_multi_stop("plsRmulti requires at least two observations")
}
if (anyNA(dataY)) {
pls_multi_stop("plsRmulti does not support missing values in Y in this experimental release")
}

storage.mode(dataY) <- "double"
if (is.null(colnames(dataY))) {
colnames(dataY) <- paste("Y", seq_len(ncol(dataY)), sep = ".")
}
if (is.null(rownames(dataY))) {
rownames(dataY) <- seq_len(nrow(dataY))
}

if (any(apply(dataY, 2, function(col) max(col) - min(col) <= tol_Xi))) {
pls_multi_stop("plsRmulti does not support constant or near-constant response columns in this experimental release")
}

dataY
}

pls_multi_prepare_predictors <- function(dataX, n_expected, tol_Xi) {
if (is.data.frame(dataX)) {
non_numeric <- !vapply(dataX, is.numeric, logical(1))
if (any(non_numeric)) {
pls_multi_stop("plsRmulti default interface requires numeric predictors")
}
dataX <- as.matrix(dataX)
} else {
dataX <- as.matrix(dataX)
}

if (!is.numeric(dataX)) {
pls_multi_stop("plsRmulti default interface requires numeric predictors")
}
if (nrow(dataX) != n_expected) {
pls_multi_stop("dataY and dataX must have the same number of rows")
}
if (ncol(dataX) < 1L) {
pls_multi_stop("plsRmulti requires at least one predictor column")
}
if (anyNA(dataX)) {
pls_multi_stop("plsRmulti does not support missing values in X in this experimental release")
}

storage.mode(dataX) <- "double"
if (is.null(colnames(dataX))) {
colnames(dataX) <- paste("X", seq_len(ncol(dataX)), sep = ".")
}
if (is.null(rownames(dataX))) {
rownames(dataX) <- seq_len(nrow(dataX))
}

if (any(apply(dataX, 2, function(col) max(col) - min(col) <= tol_Xi))) {
pls_multi_stop("plsRmulti does not support constant or near-constant predictor columns in this experimental release")
}

dataX
}

pls_multi_scale_matrix <- function(x, scale_matrix) {
x <- as.matrix(x)
if (isTRUE(scale_matrix)) {
center <- colMeans(x)
centered <- sweep(x, 2, center, "-")
scale_vec <- sqrt(colSums(centered^2) / (nrow(x) - 1))
scaled <- sweep(centered, 2, scale_vec, "/")
} else {
center <- rep(0, ncol(x))
scale_vec <- rep(1, ncol(x))
scaled <- x
}

attr(scaled, "scaled:center") <- center
attr(scaled, "scaled:scale") <- scale_vec
scaled
}

pls_multi_component_step <- function(x_resid, y_resid, tol_Xi) {
cross_cov <- crossprod(x_resid, y_resid)
if (sum(cross_cov^2) <= tol_Xi) {
return(NULL)
}

sv <- svd(cross_cov, nu = 0L, nv = 1L)
ww_raw <- drop(cross_cov %*% sv$v[, 1, drop = FALSE])
norm_ww <- sqrt(sum(ww_raw^2))
if (!is.finite(norm_ww) || norm_ww <= tol_Xi) {
return(NULL)
}

ww_norm <- ww_raw / norm_ww
sign_idx <- which.max(abs(ww_norm))
sign_value <- if (ww_norm[sign_idx] < 0) -1 else 1
ww_raw <- ww_raw * sign_value
ww_norm <- ww_norm * sign_value

t_score <- drop(x_resid %*% ww_norm)
t_norm <- drop(crossprod(t_score))
if (!is.finite(t_norm) || t_norm <= tol_Xi) {
return(NULL)
}

p_loading <- drop(crossprod(x_resid, t_score) / t_norm)
c_loading <- drop(crossprod(y_resid, t_score) / t_norm)

list(
ww = ww_raw,
wwnorm = ww_norm,
t = t_score,
p = p_loading,
c = c_loading
)
}

pls_multi_compute_wwetoile <- function(wwnorm, pp) {
if (!ncol(wwnorm)) {
return(matrix(numeric(0), nrow = nrow(wwnorm), ncol = 0L))
}

wwnorm %*% solve(crossprod(pp, wwnorm))
}

pls_multi_scaled_coefficients <- function(object, comps = object$computed_nt) {
if (comps < 1L) {
return(matrix(0, nrow = object$nc, ncol = object$ny))
}

wwnorm <- object$wwnorm[, seq_len(comps), drop = FALSE]
pp <- object$pp[, seq_len(comps), drop = FALSE]
coeff_c <- object$CoeffC[, seq_len(comps), drop = FALSE]

wwnorm %*% solve(crossprod(pp, wwnorm)) %*% t(coeff_c)
}

pls_multi_original_coefficients <- function(object, comps = object$computed_nt) {
b_scaled <- pls_multi_scaled_coefficients(object, comps = comps)
b_original <- sweep(b_scaled, 1, attr(object$ExpliX, "scaled:scale"), "/")
b_original <- sweep(b_original, 2, attr(object$RepY, "scaled:scale"), "*")
colnames(b_original) <- object$response_names
rownames(b_original) <- object$predictor_names
b_original
}

pls_multi_intercept <- function(object, comps = object$computed_nt) {
b_original <- pls_multi_original_coefficients(object, comps = comps)
drop(attr(object$RepY, "scaled:center") - crossprod(attr(object$ExpliX, "scaled:center"), b_original))
}

pls_multi_predict_scores_matrix <- function(x_scaled, object, comps) {
x_scaled %*% object$wwetoile[, seq_len(comps), drop = FALSE]
}

pls_multi_predict_response_scaled <- function(x_scaled, object, comps) {
x_scaled %*% pls_multi_scaled_coefficients(object, comps = comps)
}

pls_multi_backtransform_response <- function(y_scaled, object) {
y_original <- sweep(y_scaled, 2, attr(object$RepY, "scaled:scale"), "*")
y_original <- sweep(y_original, 2, attr(object$RepY, "scaled:center"), "+")
colnames(y_original) <- object$response_names
y_original
}

pls_multi_component_labels <- function(ncomp) {
paste("Comp_", seq_len(ncomp), sep = "")
}

pls_multi_predict_response_path_scaled <- function(x_scaled, object,
                                                   max_comps = object$computed_nt) {
if (length(max_comps) != 1L || is.na(max_comps)) {
pls_multi_stop("max_comps must be a single positive integer")
}
max_comps <- as.integer(max_comps)
if (max_comps < 1L || max_comps > object$computed_nt) {
pls_multi_stop("Cannot predict using more components than extracted")
}

pred_path <- array(
0,
dim = c(nrow(x_scaled), object$ny, max_comps),
dimnames = list(
rownames(x_scaled),
object$response_names,
pls_multi_component_labels(max_comps)
)
)

for (comp in seq_len(max_comps)) {
pred_path[, , comp] <- pls_multi_predict_response_scaled(
x_scaled = x_scaled,
object = object,
comps = comp
)
}

pred_path
}

pls_multi_backtransform_response_path <- function(y_scaled_path, object) {
path_dims <- dim(y_scaled_path)
y_original_path <- array(
0,
dim = path_dims,
dimnames = dimnames(y_scaled_path)
)

for (comp in seq_len(path_dims[3])) {
y_original_path[, , comp] <- pls_multi_backtransform_response(
y_scaled_path[, , comp, drop = FALSE][, , 1],
object = object
)
}

y_original_path
}

pls_multi_compute_fit_metrics <- function(repY, fitted_path_scaled) {
computed_nt <- dim(fitted_path_scaled)[3]
rss_by_response <- matrix(
0,
nrow = computed_nt + 1L,
ncol = ncol(repY),
dimnames = list(
paste("Nb_Comp_", 0:computed_nt, sep = ""),
colnames(repY)
)
)
rss_by_response[1, ] <- colSums(repY^2)

for (comp in seq_len(computed_nt)) {
rss_by_response[comp + 1L, ] <- colSums((repY - fitted_path_scaled[, , comp])^2)
}

rss_total <- rowSums(rss_by_response)
r2_by_response <- 1 - sweep(
rss_by_response[-1, , drop = FALSE],
2,
rss_by_response[1, ],
"/"
)
r2_total <- 1 - rss_total[-1] / rss_total[1]

list(
rss_by_response = rss_by_response,
rss_total = rss_total,
r2_by_response = r2_by_response,
r2_total = r2_total
)
}

pls_multi_vectorize_coefficients <- function(coef_matrix, intercept,
                                             response_names, predictor_names) {
coef_matrix <- as.matrix(coef_matrix)
intercept <- as.numeric(intercept)
if (length(intercept) != length(response_names)) {
pls_multi_stop("intercept must match the number of responses")
}

coef_with_intercept <- rbind(Intercept = intercept, coef_matrix)
coef_names <- rownames(coef_with_intercept)
if (is.null(coef_names)) {
coef_names <- c("Intercept", predictor_names)
rownames(coef_with_intercept) <- coef_names
}

vec <- matrix(
as.vector(coef_with_intercept),
ncol = 1L,
dimnames = list(
unlist(lapply(response_names, function(resp) {
paste(resp, coef_names, sep = ":")
}), use.names = FALSE),
NULL
)
)
vec
}

pls_multi_standardized_boot_vector <- function(object) {
pls_multi_vectorize_coefficients(
coef_matrix = object$Std.Coeffs,
intercept = rep(0, object$ny),
response_names = object$response_names,
predictor_names = object$predictor_names
)
}

pls_multi_boot_thresholds <- function(coef_vector, stabvalue) {
matrix(
stabvalue * pmax(abs(as.numeric(coef_vector)), 1),
ncol = 1L,
dimnames = dimnames(coef_vector)
)
}

pls_multi_fit_core <- function(dataY, dataX, nt, scaleX, scaleY, tol_Xi,
                               verbose = TRUE) {
dataY <- pls_multi_prepare_response(dataY, tol_Xi = tol_Xi)
dataX <- pls_multi_prepare_predictors(dataX, n_expected = nrow(dataY), tol_Xi = tol_Xi)
nt <- pls_multi_validate_nt(nt)

RepY <- pls_multi_scale_matrix(dataY, scale_matrix = scaleY)
ExpliX <- pls_multi_scale_matrix(dataX, scale_matrix = scaleX)

response_names <- colnames(RepY)
predictor_names <- colnames(ExpliX)

x_resid <- ExpliX
y_resid <- RepY

ww <- matrix(numeric(0), nrow = ncol(ExpliX), ncol = 0L)
wwnorm <- matrix(numeric(0), nrow = ncol(ExpliX), ncol = 0L)
tt <- matrix(numeric(0), nrow = nrow(ExpliX), ncol = 0L)
pp <- matrix(numeric(0), nrow = ncol(ExpliX), ncol = 0L)
CoeffC <- matrix(numeric(0), nrow = ncol(RepY), ncol = 0L)

computed_nt <- 0L
for (kk in seq_len(nt)) {
step <- pls_multi_component_step(x_resid, y_resid, tol_Xi = tol_Xi)
if (is.null(step)) {
if (verbose) {
cat("Warning only ", computed_nt, " components could thus be extracted\n", sep = "")
}
break
}

computed_nt <- kk
ww <- cbind(ww, step$ww)
wwnorm <- cbind(wwnorm, step$wwnorm)
tt <- cbind(tt, step$t)
pp <- cbind(pp, step$p)
CoeffC <- cbind(CoeffC, step$c)

x_resid <- x_resid - tcrossprod(step$t, step$p)
y_resid <- y_resid - tcrossprod(step$t, step$c)
}

if (computed_nt < 1L) {
pls_multi_stop("plsRmulti could not extract any component from the supplied data")
}

colnames(ww) <- colnames(wwnorm) <- colnames(tt) <- colnames(pp) <- colnames(CoeffC) <-
paste("Comp_", seq_len(computed_nt), sep = "")
rownames(ww) <- rownames(wwnorm) <- rownames(pp) <- predictor_names
rownames(CoeffC) <- response_names

wwetoile <- pls_multi_compute_wwetoile(wwnorm, pp)
colnames(wwetoile) <- paste("Comp_", seq_len(computed_nt), sep = "")
rownames(wwetoile) <- predictor_names

fitted_scaled <- pls_multi_predict_response_scaled(ExpliX, list(
wwnorm = wwnorm,
pp = pp,
CoeffC = CoeffC,
computed_nt = computed_nt,
nc = ncol(ExpliX),
ny = ncol(RepY)
), comps = computed_nt)

fit_object <- list(
nr = nrow(ExpliX),
nc = ncol(ExpliX),
ny = ncol(RepY),
nt = nt,
computed_nt = computed_nt,
ww = ww,
wwnorm = wwnorm,
wwetoile = wwetoile,
tt = tt,
pp = pp,
CoeffC = CoeffC,
RepY = RepY,
YChapeau = NULL,
residY = y_resid,
ExpliX = ExpliX,
Coeffs = NULL,
Std.Coeffs = NULL,
CoeffConstante = NULL,
ValsPredictY = NULL,
Std.ValsPredictY = NULL,
typeVC = "none",
dataX = dataX,
dataY = dataY,
na.miss.X = FALSE,
na.miss.Y = FALSE,
modele = "pls",
response_names = response_names,
predictor_names = predictor_names,
verbose = verbose
)

fit_object$fitted_path_scaled <- pls_multi_predict_response_path_scaled(
x_scaled = ExpliX,
object = fit_object,
max_comps = computed_nt
)
fit_object$fitted_path <- pls_multi_backtransform_response_path(
y_scaled_path = fit_object$fitted_path_scaled,
object = fit_object
)
fit_metrics <- pls_multi_compute_fit_metrics(
repY = RepY,
fitted_path_scaled = fit_object$fitted_path_scaled
)

fit_object$Std.Coeffs <- pls_multi_scaled_coefficients(fit_object, comps = computed_nt)
fit_object$Coeffs <- pls_multi_original_coefficients(fit_object, comps = computed_nt)
fit_object$CoeffConstante <- pls_multi_intercept(fit_object, comps = computed_nt)
fit_object$Std.ValsPredictY <- fit_object$fitted_path_scaled[, , computed_nt, drop = FALSE][, , 1]
fit_object$ValsPredictY <- fit_object$fitted_path[, , computed_nt, drop = FALSE][, , 1]
fit_object$YChapeau <- fit_object$ValsPredictY
fit_object$RSS_by_response <- fit_metrics$rss_by_response
fit_object$RSS_total <- fit_metrics$rss_total
fit_object$R2_by_response <- fit_metrics$r2_by_response
fit_object$R2_total <- fit_metrics$r2_total

colnames(fit_object$Std.Coeffs) <- response_names
rownames(fit_object$Std.Coeffs) <- predictor_names
names(fit_object$CoeffConstante) <- response_names
colnames(fit_object$ValsPredictY) <- response_names
colnames(fit_object$Std.ValsPredictY) <- response_names
rownames(fit_object$ValsPredictY) <- rownames(dataX)
rownames(fit_object$Std.ValsPredictY) <- rownames(dataX)

fit_object
}

pls_multi_build_matrix_from_formula <- function(formula, data, subset,
                                                weights, contrasts) {
if (missing(data)) {
data <- environment(formula)
}

mf_args <- list(
formula = formula,
data = data,
drop.unused.levels = TRUE,
na.action = na.pass
)
if (!is.null(subset)) {
mf_args$subset <- subset
}
if (!is.null(weights)) {
mf_args$weights <- weights
}
mf <- do.call(model.frame, mf_args)

weights_value <- as.vector(model.weights(mf))
if (!is.null(weights_value)) {
pls_multi_stop("plsRmulti does not support weights in this experimental release")
}

terms_full <- attr(mf, "terms")
terms_x <- stats::delete.response(terms_full)
dataY <- model.response(mf, "numeric")
dataX <- model.matrix(terms_x, mf, contrasts.arg = contrasts)
if (attr(terms_x, "intercept") > 0L) {
dataX <- dataX[, -1, drop = FALSE]
}

list(
dataY = dataY,
dataX = dataX,
terms_x = terms_x,
contrasts = attr(dataX, "contrasts"),
xlevels = lapply(mf[vapply(mf, is.factor, logical(1))], levels)
)
}

pls_multi_prepare_newdata <- function(object, newdata) {
if (missing(newdata) || is.null(newdata)) {
return(object$ExpliX)
}

if (!is.null(object$terms_x)) {
mf <- model.frame(
object$terms_x,
data = newdata,
na.action = na.pass,
xlev = object$xlevels,
drop.unused.levels = TRUE
)
newdata_matrix <- model.matrix(object$terms_x, mf, contrasts.arg = object$contrasts)
if (attr(object$terms_x, "intercept") > 0L) {
newdata_matrix <- newdata_matrix[, -1, drop = FALSE]
}
} else {
if (is.data.frame(newdata)) {
non_numeric <- !vapply(newdata, is.numeric, logical(1))
if (any(non_numeric)) {
pls_multi_stop("plsRmulti default predict() requires numeric predictors")
}
newdata_matrix <- as.matrix(newdata)
} else {
newdata_matrix <- as.matrix(newdata)
}
}

if (!is.numeric(newdata_matrix)) {
pls_multi_stop("plsRmulti predict() requires numeric predictors")
}
if (ncol(newdata_matrix) != object$nc) {
pls_multi_stop("newdata must contain the same number of predictor columns used at fit time")
}
if (anyNA(newdata_matrix)) {
pls_multi_stop("plsRmulti predict() does not support missing values in this experimental release")
}

if (!is.null(colnames(newdata_matrix)) && !is.null(object$predictor_names)) {
if (!all(object$predictor_names %in% colnames(newdata_matrix))) {
pls_multi_stop("newdata is missing predictor columns required by the fitted plsRmulti model")
}
newdata_matrix <- newdata_matrix[, object$predictor_names, drop = FALSE]
}

storage.mode(newdata_matrix) <- "double"
newdata_matrix
}

pls_multi_scale_newdata <- function(newdata_matrix, object) {
if (isTRUE(object$scaleX)) {
sweep(
sweep(newdata_matrix, 2, attr(object$ExpliX, "scaled:center")),
2,
attr(object$ExpliX, "scaled:scale"),
"/"
)
} else {
newdata_matrix
}
}
