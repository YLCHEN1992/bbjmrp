#TwoSampleMR code
mr_median=function (dat, parameters = default_parameters())
{
if ("mr_keep" %in% names(dat))
dat <- subset(dat, mr_keep)
if (nrow(dat) < 3) {
warning("Need at least 3 SNPs")
return(NULL)
}
b_exp <- dat$beta.exposure
b_out <- dat$beta.outcome
se_exp <- dat$se.exposure
se_out <- dat$se.outcome
sm <- mr_simple_median(b_exp, b_out, se_exp, se_out, parameters)
wm <- mr_weighted_median(b_exp, b_out, se_exp, se_out, parameters)
pm <- mr_penalised_weighted_median(b_exp, b_out, se_exp,
se_out, parameters)
res <- data.frame(id.exposure = dat$id.exposure[1], id.outcome = dat$id.outcome[1],
method = c("Simple median", "Weighted median", "Penalised median"),
nsnp = length(b_exp), b = c(sm$b, wm$b, pm$b), se = c(sm$se,
wm$se, pm$se), stringsAsFactors = FALSE)
res$ci_low <- res$b - stats::qnorm(1 - parameters$alpha/2) *
res$se
res$ci_upp <- res$b + stats::qnorm(1 - parameters$alpha/2) *
res$se
res$pval <- c(sm$pval, wm$pval, pm$pval)
return(res)
}
mr_weighted_median=function (b_exp, b_out, se_exp, se_out, parameters = default_parameters())
{
if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) &
!is.na(se_out)) < 3)
return(list(b = NA, se = NA, pval = NA, nsnp = NA))
b_iv <- b_out/b_exp
VBj <- ((se_out)^2)/(b_exp)^2 + (b_out^2) * ((se_exp^2))/(b_exp)^4
b <- weighted_median(b_iv, 1/VBj)
se <- weighted_median_bootstrap(b_exp, b_out, se_exp, se_out,
1/VBj, parameters$nboot)
pval <- 2 * stats::pnorm(abs(b/se), lower.tail = FALSE)
return(list(b = b, se = se, pval = pval, Q = NA, Q_df = NA,
Q_pval = NA, nsnp = length(b_exp)))
}
weighted_median_bootstrap=function (b_exp, b_out, se_exp, se_out, weights, nboot)
{
med <- rep(0, nboot)
for (i in 1:nboot) {
b_exp.boot = stats::rnorm(length(b_exp), mean = b_exp,
sd = se_exp)
b_out.boot = stats::rnorm(length(b_out), mean = b_out,
sd = se_out)
betaIV.boot = b_out.boot/b_exp.boot
med[i] = weighted_median(betaIV.boot, weights)
}
return(stats::sd(med))
}
weighted_median=function (b_iv, weights)
{
betaIV.order <- b_iv[order(b_iv)]
weights.order <- weights[order(b_iv)]
weights.sum <- cumsum(weights.order) - 0.5 * weights.order
weights.sum <- weights.sum/sum(weights.order)
below <- max(which(weights.sum < 0.5))
b = betaIV.order[below] + (betaIV.order[below + 1] - betaIV.order[below]) *
(0.5 - weights.sum[below])/(weights.sum[below + 1] -
weights.sum[below])
return(b)
}
default_parameters=function ()
{
list(test_dist = "z", nboot = 1000, Cov = 0, penk = 20, phi = 1,
alpha = 0.05, Qthresh = 0.05, over.dispersion = TRUE,
loss.function = "huber", shrinkage = FALSE)
}
mr_simple_median=function (b_exp, b_out, se_exp, se_out, parameters = default_parameters())
{
if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) &
!is.na(se_out)) < 3)
return(list(b = NA, se = NA, pval = NA, nsnp = NA))
b_iv <- b_out/b_exp
b <- weighted_median(b_iv, rep(1/length(b_exp), length(b_exp)))
se <- weighted_median_bootstrap(b_exp, b_out, se_exp, se_out,
rep(1/length(b_exp), length(b_exp)), parameters$nboot)
pval <- 2 * stats::pnorm(abs(b/se), lower.tail = FALSE)
return(list(b = b, se = se, pval = pval, nsnp = length(b_exp)))
}
mr_weighted_mode=function (b_exp, b_out, se_exp, se_out, parameters = default_parameters())
{
index <- !is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) &
!is.na(se_out)
if (sum(index) < 3)
return(list(b = NA, se = NA, pval = NA, nsnp = NA))
b_exp <- b_exp[index]
b_out <- b_out[index]
se_exp <- se_exp[index]
se_out <- se_out[index]
return(mr_mode(data.frame(beta.exposure = b_exp, beta.outcome = b_out,
se.exposure = se_exp, se.outcome = se_out), parameters = parameters,
mode_method = "Weighted mode"))
}
mr_mode=function (dat, parameters = default_parameters(), mode_method = "all")
{
if ("mr_keep" %in% names(dat))
dat <- subset(dat, mr_keep)
if (nrow(dat) < 3) {
warning("Need at least 3 SNPs")
return(NULL)
}
b_exp <- dat$beta.exposure
b_out <- dat$beta.outcome
se_exp <- dat$se.exposure
se_out <- dat$se.outcome
beta <- function(BetaIV.in, seBetaIV.in, phi) {
s <- 0.9 * (min(stats::sd(BetaIV.in), stats::mad(BetaIV.in)))/length(BetaIV.in)^(1/5)
weights <- seBetaIV.in^-2/sum(seBetaIV.in^-2)
beta <- NULL
for (cur_phi in phi) {
h <- max(1e-08, s * cur_phi)
densityIV <- stats::density(BetaIV.in, weights = weights,
bw = h)
beta[length(beta) + 1] <- densityIV$x[densityIV$y ==
max(densityIV$y)]
}
return(beta)
}
boot <- function(BetaIV.in, seBetaIV.in, beta_Mode.in, nboot) {
beta.boot <- matrix(nrow = nboot, ncol = length(beta_Mode.in))
for (i in 1:nboot) {
BetaIV.boot <- stats::rnorm(length(BetaIV.in), mean = BetaIV.in,
sd = seBetaIV.in[, 1])
BetaIV.boot_NOME <- stats::rnorm(length(BetaIV.in),
mean = BetaIV.in, sd = seBetaIV.in[, 2])
beta.boot[i, 1:length(phi)] <- beta(BetaIV.in = BetaIV.boot,
seBetaIV.in = rep(1, length(BetaIV)), phi = phi)
beta.boot[i, (length(phi) + 1):(2 * length(phi))] <- beta(BetaIV.in = BetaIV.boot,
seBetaIV.in = seBetaIV.in[, 1], phi = phi)
weights <- 1/seBetaIV.in[, 1]^2
penalty <- stats::pchisq(weights * (BetaIV.boot -
beta.boot[i, (length(phi) + 1):(2 * length(phi))])^2,
df = 1, lower.tail = FALSE)
pen.weights <- weights * pmin(1, penalty * parameters$penk)
beta.boot[i, (2 * length(phi) + 1):(3 * length(phi))] <- beta(BetaIV.in = BetaIV.boot,
seBetaIV.in = sqrt(1/pen.weights), phi = phi)
beta.boot[i, (3 * length(phi) + 1):(4 * length(phi))] <- beta(BetaIV.in = BetaIV.boot_NOME,
seBetaIV.in = rep(1, length(BetaIV)), phi = phi)
beta.boot[i, (4 * length(phi) + 1):(5 * length(phi))] <- beta(BetaIV.in = BetaIV.boot_NOME,
seBetaIV.in = seBetaIV.in[, 2], phi = phi)
}
return(beta.boot)
}
phi <- parameters$phi
nboot <- parameters$nboot
alpha <- parameters$alpha
BetaIV <- b_out/b_exp
seBetaIV <- cbind(sqrt((se_out^2)/(b_exp^2) + ((b_out^2) *
(se_exp^2))/(b_exp^4)), se_out/abs(b_exp))
beta_SimpleMode <- beta(BetaIV.in = BetaIV, seBetaIV.in = rep(1,
length(BetaIV)), phi = phi)
beta_WeightedMode <- beta(BetaIV.in = BetaIV, seBetaIV.in = seBetaIV[,
1], phi = phi)
weights <- 1/seBetaIV[, 1]^2
penalty <- stats::pchisq(weights * (BetaIV - beta_WeightedMode)^2,
df = 1, lower.tail = FALSE)
pen.weights <- weights * pmin(1, penalty * parameters$penk)
beta_PenalisedMode <- beta(BetaIV.in = BetaIV, seBetaIV.in = sqrt(1/pen.weights),
phi = phi)
beta_WeightedMode_NOME <- beta(BetaIV.in = BetaIV, seBetaIV.in = seBetaIV[,
2], phi = phi)
beta_Mode <- rep(c(beta_SimpleMode, beta_WeightedMode, beta_PenalisedMode,
beta_SimpleMode, beta_WeightedMode_NOME))
beta_Mode.boot <- boot(BetaIV.in = BetaIV, seBetaIV.in = seBetaIV,
beta_Mode.in = beta_Mode, nboot = nboot)
se_Mode <- apply(beta_Mode.boot, 2, stats::mad)
CIlow_Mode <- beta_Mode - stats::qnorm(1 - alpha/2) * se_Mode
CIupp_Mode <- beta_Mode + stats::qnorm(1 - alpha/2) * se_Mode
P_Mode <- stats::pt(abs(beta_Mode/se_Mode), df = length(b_exp) -
1, lower.tail = FALSE) * 2
Method <- rep(c("Simple mode", "Weighted mode", "Penalised mode",
"Simple mode (NOME)", "Weighted mode (NOME)"), each = length(phi))
id.exposure <- ifelse("id.exposure" %in% names(dat), dat$id.exposure[1],
"")
id.outcome <- ifelse("id.outcome" %in% names(dat), dat$id.outcome[1],
"")
Results <- data.frame(id.exposure = id.exposure, id.outcome = id.outcome,
method = Method, nsnp = length(b_exp), b = beta_Mode,
se = se_Mode, ci_low = CIlow_Mode, ci_upp = CIupp_Mode,
pval = P_Mode, stringsAsFactors = FALSE)
if (mode_method == "all") {
return(Results)
}
else {
stopifnot(all(mode_method %in% Results$method))
i <- which(Results$method == mode_method)
return(list(b = Results$b[i], se = Results$se[i], pval = Results$pval[i],
nsnp = length(b_exp)))
}
return(Results)
}