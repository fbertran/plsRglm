.plsRglm_backend_warnings <- new.env(parent = emptyenv())

pls_glm_normalize_fit_backend <- function(fit_backend = "stats") {
match.arg(fit_backend, c("stats", "fastglm"))
}

pls_glm_warn_once <- function(key, message) {
if (!exists(key, envir = .plsRglm_backend_warnings, inherits = FALSE)) {
assign(key, TRUE, envir = .plsRglm_backend_warnings)
warning(message, call. = FALSE)
}
}

pls_glm_resolve_fit_backend <- function(fit_backend = "stats",
                                        modele,
                                        pvals.expli,
                                        has_missing,
                                        has_weights) {
fit_backend <- pls_glm_normalize_fit_backend(fit_backend)
if (fit_backend != "fastglm") {
return(fit_backend)
}

reasons <- character()
if (!requireNamespace("fastglm", quietly = TRUE)) {
reasons <- c(reasons, "`fastglm` is not installed")
}
if (identical(modele, "pls-glm-polr")) {
reasons <- c(reasons, "`pls-glm-polr` is not supported")
}
if (identical(modele, "pls")) {
reasons <- c(reasons, "`pls` does not use a GLM backend")
}
if (isTRUE(pvals.expli)) {
reasons <- c(reasons, "`pvals.expli = TRUE` is not validated for `fastglm`")
}
if (isTRUE(has_missing)) {
reasons <- c(reasons, "missing-data fits use the compatibility path")
}
if (isTRUE(has_weights)) {
reasons <- c(reasons, "weighted fits use the compatibility path")
}

if (length(reasons) == 0L) {
return("fastglm")
}

key <- paste("fit_backend", modele, paste(reasons, collapse = "|"), sep = "::")
pls_glm_warn_once(
key,
paste0(
"Falling back to `fit_backend = \"stats\"` because ",
paste(reasons, collapse = "; "),
"."
)
)
"stats"
}

pls_glm_build_design_matrix <- function(tt = NULL, xcol = NULL) {
if (is.null(tt)) {
base <- matrix(1, nrow = length(xcol), ncol = 1)
} else {
tt <- as.matrix(tt)
base <- cbind(`(Intercept)` = 1, tt)
}

if (is.null(xcol)) {
base
} else {
cbind(base, xcol)
}
}

pls_glm_fit_matrix <- function(x, y, family, fit_backend = "stats") {
if (identical(fit_backend, "fastglm")) {
return(fastglm::fastglm(x = x, y = y, family = family))
}

fit <- stats::glm.fit(x = x, y = y, family = family)
class(fit) <- c("glm", "lm")
fit
}

pls_glm_format_score_matrix <- function(tt) {
tt <- as.matrix(tt)
if (is.null(colnames(tt))) {
colnames(tt) <- paste("Comp_", seq_len(ncol(tt)), sep = "")
}
tt
}

pls_glm_score_design_matrix <- function(tt) {
tt <- pls_glm_format_score_matrix(tt)
x <- cbind(`(Intercept)` = 1, tt)
colnames(x)[-1] <- paste0("tt.", colnames(tt))
x
}

pls_glm_fit_intercept_model <- function(y, family, fit_backend = "stats") {
x <- matrix(1, nrow = length(y), ncol = 1)
colnames(x) <- "(Intercept)"
fit <- pls_glm_fit_matrix(x = x, y = as.numeric(y), family = family, fit_backend = fit_backend)
fit$pls_score_names <- character()
fit
}

pls_glm_fit_score_model <- function(y, tt, family, fit_backend = "stats") {
tt <- pls_glm_format_score_matrix(tt)
fit <- pls_glm_fit_matrix(
x = pls_glm_score_design_matrix(tt),
y = as.numeric(y),
family = family,
fit_backend = fit_backend
)
fit$pls_score_names <- colnames(tt)
fit
}

pls_glm_predict_score_model <- function(fit, tt, type = c("link", "response", "terms")) {
type <- match.arg(type)
tt <- pls_glm_format_score_matrix(tt)
design <- pls_glm_score_design_matrix(tt)
beta <- stats::coef(fit)
beta <- beta[colnames(design)]
eta <- drop(design %*% beta)

if (type == "link") {
return(eta)
}
if (type == "response") {
return(fit$family$linkinv(eta))
}

terms_matrix <- sweep(tt, 2, beta[-1], `*`)
colnames(terms_matrix) <- colnames(design)[-1]
attr(terms_matrix, "constant") <- unname(beta[1])
terms_matrix
}

pls_glm_refit_compatible_model <- function(y, tt, family) {
tt <- pls_glm_format_score_matrix(tt)
stats::glm(y ~ ., data = data.frame(y = as.numeric(y), tt = tt), family = family)
}

pls_glm_fit_weight_column <- function(tt = NULL,
                                      xcol,
                                      y,
                                      family,
                                      fit_backend = "stats",
                                      compute_pvalue = FALSE) {
y <- as.numeric(y)
ok <- !is.na(xcol) & !is.na(y)
design <- pls_glm_build_design_matrix(
tt = if (is.null(tt)) NULL else as.matrix(tt)[ok, , drop = FALSE],
xcol = xcol[ok]
)
fit <- pls_glm_fit_matrix(design, y[ok], family = family, fit_backend = fit_backend)
coef_index <- ncol(design)
coef_value <- unname(stats::coef(fit)[coef_index])

if (!compute_pvalue) {
return(list(coef = coef_value, fit = fit))
}

summary_coef <- summary(fit)$coefficients
list(
coef = coef_value,
p_value = summary_coef[coef_index, ncol(summary_coef), drop = TRUE],
fit = fit
)
}
