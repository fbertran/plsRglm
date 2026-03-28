pls_glm_stop_on_all_missing <- function(x, margin, message) {
if (any(apply(is.na(x), MARGIN = margin, all))) {
stop(message, call. = FALSE)
}
}

pls_glm_apply_sparse_filter <- function(tempww, tempvalpvalstep,
                                        alpha.pvals.expli, sparse,
                                        sparseStop) {
temppvalstep <- (tempvalpvalstep < alpha.pvals.expli)
break_nt_sparse <- FALSE

if (sparse) {
if (sum(temppvalstep) > 1L) {
tempww[!temppvalstep] <- 0
}
}

if (sparseStop) {
if (sum(temppvalstep) == 0L) {
break_nt_sparse <- TRUE
}
}

list(
tempww = tempww,
tempvalpvalstep = tempvalpvalstep,
temppvalstep = temppvalstep,
break_nt_sparse = break_nt_sparse
)
}

pls_glm_as_polr_response <- function(y) {
if (is.matrix(y) || is.data.frame(y)) {
y <- y[, 1]
}
if (is.ordered(y)) {
return(as.ordered(y))
}
as.factor(y)
}

pls_glm_validate_inputs <- function(dataY, dataX, dataPredictY,
                                    weights, weights_missing,
                                    method, method_missing,
                                    modele, family,
                                    sparse, sparseStop,
                                    pvals.expli, naive, naive_missing,
                                    fit_backend = "stats",
                                    verbose,
                                    weighted_naive_message = "Only naive DoF can be used with weighted PLS\n",
                                    eval_env = parent.frame()) {
PredYisdataX <- identical(dataPredictY, dataX)

pls_glm_stop_on_all_missing(
dataX,
margin = 2,
message = "One of the columns of dataX is completely filled with missing data"
)
pls_glm_stop_on_all_missing(
dataX,
margin = 1,
message = "One of the rows of dataX is completely filled with missing data"
)
if (!PredYisdataX) {
pls_glm_stop_on_all_missing(
dataPredictY,
margin = 1,
message = "One of the rows of dataPredictY is completely filled with missing data"
)
}

if (weights_missing) {
NoWeights <- TRUE
} else {
NoWeights <- all(weights == rep(1, length(dataY)))
}
if (method_missing) {
method <- "logistic"
}

na.miss.X <- any(is.na(dataX))
na.miss.Y <- any(is.na(dataY))
na.miss.PredictY <- any(is.na(dataPredictY))

if (is.null(modele)) {
naive <- FALSE
} else {
if (modele == "pls") {
naive <- FALSE
} else {
if (!naive_missing) {
if (verbose) {
cat("Only naive DoF can be used with PLS GLM\n")
}
}
naive <- TRUE
}
}

if (na.miss.X | na.miss.Y) {
naive <- TRUE
if (verbose) {
cat("Only naive DoF can be used with missing data\n")
}
if (!NoWeights) {
if (verbose) {
cat("Weights cannot be used with missing data\n")
}
}
if (sparse) {
if (verbose) {
cat("sparse option cannot be used with missing data\n")
}
sparse <- FALSE
}
}

if (!NoWeights) {
naive <- TRUE
if (verbose) {
cat(weighted_naive_message)
}
} else {
NoWeights <- TRUE
}

if (sparse) {
sparseStop <- TRUE
}
if (sparseStop) {
pvals.expli <- TRUE
}

if (!is.data.frame(dataX)) {
dataX <- data.frame(dataX)
}

if (is.null(modele) & !is.null(family)) {
modele <- "pls-glm-family"
}
allowed_modeles <- c(
"pls", "pls-glm-logistic", "pls-glm-family", "pls-glm-Gamma",
"pls-glm-gaussian", "pls-glm-inverse.gaussian", "pls-glm-poisson",
"pls-glm-polr"
)
if (!(modele %in% allowed_modeles)) {
print(modele)
stop("'modele' not recognized")
}
if (!(modele %in% "pls-glm-family") & !is.null(family)) {
stop("Set 'modele=pls-glm-family' to use the family option")
}
if (modele == "pls") {
family <- NULL
}
if (modele == "pls-glm-Gamma") {
family <- Gamma(link = "inverse")
}
if (modele == "pls-glm-gaussian") {
family <- gaussian(link = "identity")
}
if (modele == "pls-glm-inverse.gaussian") {
family <- inverse.gaussian(link = "1/mu^2")
}
if (modele == "pls-glm-logistic") {
family <- binomial(link = "logit")
}
if (modele == "pls-glm-poisson") {
family <- poisson(link = "log")
}
if (modele == "pls-glm-polr") {
family <- NULL
}
fit_backend <- pls_glm_resolve_fit_backend(
fit_backend = fit_backend,
modele = modele,
pvals.expli = pvals.expli,
has_missing = na.miss.X | na.miss.Y,
has_weights = !NoWeights
)
if (!is.null(family)) {
if (is.character(family)) {
family <- get(family, mode = "function", envir = eval_env)
}
if (is.function(family)) {
family <- family()
}
if (is.language(family)) {
family <- eval(family, envir = eval_env)
}
}

if (modele %in% c(
"pls-glm-family", "pls-glm-Gamma", "pls-glm-gaussian",
"pls-glm-inverse.gaussian", "pls-glm-logistic", "pls-glm-poisson"
)) {
if (verbose) {
print(family)
}
}
if (modele %in% c("pls-glm-polr")) {
if (verbose) {
cat("\nModel:", modele, "\n")
cat("Method:", method, "\n\n")
}
}
if (modele == "pls") {
if (verbose) {
cat("\nModel:", modele, "\n\n")
}
}

list(
PredYisdataX = PredYisdataX,
NoWeights = NoWeights,
method = method,
na.miss.X = na.miss.X,
na.miss.Y = na.miss.Y,
na.miss.PredictY = na.miss.PredictY,
naive = naive,
sparse = sparse,
sparseStop = sparseStop,
pvals.expli = pvals.expli,
dataX = dataX,
modele = modele,
family = family,
fit_backend = fit_backend
)
}

pls_glm_preprocess_inputs <- function(dataY, dataX, dataPredictY,
                                      scaleX, scaleY,
                                      weights, NoWeights,
                                      PredYisdataX, modele) {
scaleY <- NULL
if (is.null(scaleY)) {
if (!(modele %in% c("pls"))) {
scaleY <- FALSE
} else {
scaleY <- TRUE
}
}

if (scaleY) {
if (NoWeights) {
RepY <- scale(dataY)
} else {
meanY <- weighted.mean(dataY, weights)
stdevY <- sqrt(
((length(dataY) - 1) / length(dataY)) *
weighted.mean((dataY - meanY)^2, weights)
)
RepY <- (dataY - meanY) / stdevY
attr(RepY, "scaled:center") <- meanY
attr(RepY, "scaled:scale") <- stdevY
}
} else {
RepY <- dataY
attr(RepY, "scaled:center") <- 0
attr(RepY, "scaled:scale") <- 1
}

if (scaleX) {
if (NoWeights) {
ExpliX <- scale(dataX)
} else {
meanX <- apply(dataX, 2, weighted.mean, weights)
stdevX <- sqrt(
((length(dataY) - 1) / length(dataY)) *
apply((sweep(dataX, 2, meanX))^2, 2, weighted.mean, weights)
)
ExpliX <- sweep(sweep(dataX, 2, meanX), 2, stdevX, "/")
attr(ExpliX, "scaled:center") <- meanX
attr(ExpliX, "scaled:scale") <- stdevX
}
if (PredYisdataX) {
PredictY <- ExpliX
} else {
PredictY <- sweep(
sweep(dataPredictY, 2, attr(ExpliX, "scaled:center")),
2,
attr(ExpliX, "scaled:scale"),
"/"
)
}
} else {
ExpliX <- dataX
attr(ExpliX, "scaled:center") <- rep(0, ncol(dataX))
attr(ExpliX, "scaled:scale") <- rep(1, ncol(dataX))
PredictY <- dataPredictY
}

if (is.null(colnames(ExpliX))) {
colnames(ExpliX) <- paste("X", 1:ncol(ExpliX), sep = ".")
}
if (is.null(rownames(ExpliX))) {
rownames(ExpliX) <- 1:nrow(ExpliX)
}

XXNA <- !is.na(ExpliX)
YNA <- !is.na(RepY)
if (PredYisdataX) {
PredictYNA <- XXNA
} else {
PredictYNA <- !is.na(PredictY)
}

XXwotNA <- as.matrix(ExpliX)
XXwotNA[!XXNA] <- 0

YwotNA <- as.matrix(RepY)
YwotNA[!YNA] <- 0

if (PredYisdataX) {
PredictYwotNA <- XXwotNA
} else {
PredictYwotNA <- as.matrix(PredictY)
PredictYwotNA[is.na(PredictY)] <- 0
}

list(
scaleY = scaleY,
RepY = RepY,
ExpliX = ExpliX,
PredictY = PredictY,
XXNA = XXNA,
YNA = YNA,
PredictYNA = PredictYNA,
XXwotNA = XXwotNA,
YwotNA = YwotNA,
PredictYwotNA = PredictYwotNA
)
}

pls_glm_compute_glm_weights <- function(modele, XXwotNA, XXNA, YwotNA,
                                        tt, family, method, kk,
                                        pvals.expli, alpha.pvals.expli,
                                        sparse, sparseStop,
                                        fit_backend = "stats") {
supported_modeles <- c(
"pls-glm-family", "pls-glm-Gamma", "pls-glm-gaussian",
"pls-glm-inverse.gaussian", "pls-glm-logistic", "pls-glm-poisson",
"pls-glm-polr"
)
if (!(modele %in% supported_modeles)) {
stop("Unsupported 'modele' for GLM-specific weighting", call. = FALSE)
}

tempww <- rep(0, ncol(XXwotNA))
XXglm <- XXwotNA
XXglm[!XXNA] <- NA
yglm <- drop(YwotNA)

if (modele %in% c(
"pls-glm-family", "pls-glm-Gamma", "pls-glm-gaussian",
"pls-glm-inverse.gaussian", "pls-glm-logistic", "pls-glm-poisson"
)) {
if (!pvals.expli) {
for (jj in 1:ncol(XXglm)) {
tempww[jj] <- pls_glm_fit_weight_column(
tt = tt,
xcol = XXglm[, jj],
y = yglm,
family = family,
fit_backend = fit_backend,
compute_pvalue = FALSE
)$coef
}
return(list(tempww = tempww, break_nt_sparse = FALSE))
}

tempvalpvalstep <- rep(0, ncol(XXglm))
for (jj in 1:ncol(XXglm)) {
tmww <- pls_glm_fit_weight_column(
tt = tt,
xcol = XXglm[, jj],
y = yglm,
family = family,
fit_backend = "stats",
compute_pvalue = TRUE
)
tempww[jj] <- tmww$coef
tempvalpvalstep[jj] <- tmww$p_value
}
filtered <- pls_glm_apply_sparse_filter(
tempww = tempww,
tempvalpvalstep = tempvalpvalstep,
alpha.pvals.expli = alpha.pvals.expli,
sparse = sparse,
sparseStop = sparseStop
)
return(filtered)
}

YwotNA <- pls_glm_as_polr_response(YwotNA)
requireNamespace("MASS")
tts <- tt

if (!pvals.expli) {
for (jj in 1:ncol(XXglm)) {
tempww[jj] <- -1 * MASS::polr(
YwotNA ~ cbind(tts, XXglm[, jj]),
na.action = na.exclude,
method = method
)$coef[kk]
}
return(list(tempww = tempww, break_nt_sparse = FALSE))
}

tempvalpvalstep <- rep(0, ncol(XXglm))
for (jj in 1:ncol(XXglm)) {
tmww <- -1 * summary(MASS::polr(
YwotNA ~ cbind(tts, XXglm[, jj]),
na.action = na.exclude,
Hess = TRUE,
method = method
))$coefficients[kk, ]
tempww[jj] <- tmww[1]
tempvalpvalstep[jj] <- 2 * pnorm(-abs(tmww[3]))
}

pls_glm_apply_sparse_filter(
tempww = tempww,
tempvalpvalstep = tempvalpvalstep,
alpha.pvals.expli = alpha.pvals.expli,
sparse = sparse,
sparseStop = sparseStop
)
}

pls_glm_extract_component <- function(res, tempww, XXwotNA, XXNA,
                                      PredictYNA, PredictYwotNA,
                                      PredYisdataX,
                                      na.miss.X, na.miss.Y,
                                      na.miss.PredictY,
                                      tol_Xi, kk, sparse,
                                      break_nt_sparse,
                                      break_nt_sparse1,
                                      alpha.pvals.expli,
                                      verbose) {
if (break_nt_sparse & (kk == 1L)) {
if (verbose) {
cat(paste("No significant predictors (<", alpha.pvals.expli, ") found\n", sep = ""))
cat("Warning only one standard component (without sparse option) was thus extracted\n")
}
break_nt_sparse1 <- TRUE
}
if (break_nt_sparse & !(kk == 1L)) {
res$computed_nt <- kk - 1
if (!break_nt_sparse1) {
if (verbose) {
cat(paste("No more significant predictors (<", alpha.pvals.expli, ") found\n", sep = ""))
cat(paste("Warning only ", res$computed_nt, " components were thus extracted\n", sep = ""))
}
}
return(list(
res = res,
break_nt = FALSE,
break_nt_sparse1 = break_nt_sparse1,
should_break = TRUE
))
}

component_cpp <- pls_component_step_cpp(
xxwotna_r = as.matrix(XXwotNA),
xxna_r = 1 * as.matrix(XXNA),
tempww_r = as.numeric(tempww),
prev_pp_r = if (is.null(res$pp)) NULL else as.matrix(res$pp),
predict_na_r = if ((!PredYisdataX) && (na.miss.PredictY & !na.miss.Y) && (sparse == FALSE)) {
1 * as.matrix(PredictYNA)
} else {
NULL
},
tol_xi = tol_Xi,
check_xx = isTRUE(na.miss.X & !na.miss.Y & (sparse == FALSE)),
check_predict = isTRUE((!PredYisdataX) && (na.miss.PredictY & !na.miss.Y) && (sparse == FALSE))
)
tempwwnorm <- as.numeric(component_cpp$tempwwnorm)
temptt <- as.matrix(component_cpp$temptt)
temppp <- as.numeric(component_cpp$temppp)
res$residXX <- as.matrix(component_cpp$residXX)

break_nt <- FALSE
if (na.miss.X & !na.miss.Y) {
if (sparse == FALSE && component_cpp$bad_xx_row > 0L) {
ii <- component_cpp$bad_xx_row
break_nt <- TRUE
res$computed_nt <- kk - 1
if (verbose) {
cat(paste(
"Warning : reciprocal condition number of t(cbind(res$pp,temppp)[XXNA[",
ii,
",],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[",
ii,
",],,drop=FALSE] < 10^{-12}\n",
sep = ""
))
cat(paste("Warning only ", res$computed_nt, " components could thus be extracted\n", sep = ""))
}
}
}

if ((!PredYisdataX) && (!break_nt)) {
if (na.miss.PredictY & !na.miss.Y) {
if (sparse == FALSE && component_cpp$bad_predict_row > 0L) {
ii <- component_cpp$bad_predict_row
break_nt <- TRUE
res$computed_nt <- kk - 1
if (verbose) {
cat(paste(
"Warning : reciprocal condition number of t(cbind(res$pp,temppp)[PredictYNA[",
ii,
",,drop=FALSE],])%*%cbind(res$pp,temppp)[PredictYNA[",
ii,
",,drop=FALSE],] < 10^{-12}\n",
sep = ""
))
cat(paste("Warning only ", res$computed_nt, " components could thus be extracted\n", sep = ""))
}
}
}
}

if (break_nt) {
return(list(
res = res,
break_nt = break_nt,
break_nt_sparse1 = break_nt_sparse1,
should_break = TRUE
))
}

res$ww <- cbind(res$ww, tempww)
res$wwnorm <- cbind(res$wwnorm, tempwwnorm)
res$tt <- cbind(res$tt, temptt)
res$pp <- cbind(res$pp, temppp)

list(
res = res,
break_nt = break_nt,
break_nt_sparse1 = break_nt_sparse1,
should_break = FALSE,
tempwwnorm = tempwwnorm,
temptt = temptt,
temppp = temppp
)
}
