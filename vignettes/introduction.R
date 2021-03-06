## ----shape--------------------------------------------------------------------
# load package
library(eva)

# A naive implementation of the GEV cumulative density function
pgev_naive <- function(q, loc = 0, scale = 1, shape = 1) {
    exp(-(1 + (shape * (q - loc))/scale)^(-1/shape))
}


curve(pgev_naive(1, 0, 1, x), 1e-20, .01, log = "x", n = 1025, 
      xlab = "Shape", ylab = "GEV CDF", col = "red", lty = 1, lwd = 1)
curve(pgev(1, 0, 1, x), 1e-20, .01, log = "x", n = 1025, col = "blue", lty = 6, lwd = 3, add = TRUE)

# Similarly for the GPD cdf
pgpd_naive <- function(q, loc = 0, scale = 1, shape = 1) {
    (1 - (1 + (shape * (q - loc))/scale)^(-1/shape))
}

curve(pgpd_naive(1, 0, 1, x), 1e-20, .01, log = "x", n = 1025, 
      xlab = "Shape", ylab = "GPD CDF", col = "red", lty = 1, lwd = 1)
curve(pgpd(1, 0, 1, x),  1e-20, .01, log = "x", n = 1025, col = "blue", lty = 6, lwd = 3, add = TRUE)



## ----seqtesting---------------------------------------------------------------
data(lowestoft)
gevrSeqTests(lowestoft, method = "ed")

## ----returnlevel, fig.height = 12, fig.width = 8------------------------------

# Make 250 year return level plot using gevr for r = 1 to 10 with the LoweStoft data

data(lowestoft)
result <- matrix(0, 20, 4)
period <- 250

for(i in 1:10) {
  z <- gevrFit(as.matrix(lowestoft[, 1:i]))
  y1 <- gevrRl(z, period, conf = 0.95, method = "delta")
  y2 <- gevrRl(z, period, conf = 0.95, method = "profile", plot = FALSE)
  result[i, 1] <- i
  result[i, 2] <- y1$Estimate
  result[i, 3:4] <- y1$CI
  result[(i + 10), 1] <- i
  result[(i + 10), 2] <- y2$Estimate
  result[(i + 10), 3:4] <- y2$CI
}

result <- cbind.data.frame(result, c(rep("Delta", 10), rep("Profile", 10)))
colnames(result) <- c("r", "Est", "Lower", "Upper", "Method")
result <- as.data.frame(result)

prof <- subset(result, Method == "Profile")
del <- subset(result, Method == "Delta")

par(mfrow = c(2, 1))

plot(prof$r, prof$Est, main = "Profile Likelihood", 
     xlab = "r", ylab = "250 Year Return Level",
     xlim = c(1, 10), ylim = c(4, 7))
polygon(c(rev(prof$r), prof$r), c(rev(prof$Lower), prof$Upper), col = 'grey80', border = NA)
points(prof$r, prof$Est, pch = 19, col = 'black')
lines(prof$r, prof$Est, lty = 'solid', col = 'black')
lines(prof$r, prof$Lower, lty = 'dashed', col = 'red')
lines(prof$r, prof$Upper, lty = 'dashed', col = 'red')


plot(del$r, del$Est, main = "Delta Method", 
     xlab = "r", ylab = "250 Year Return Level",
     xlim = c(1, 10), ylim = c(4, 7))
polygon(c(rev(del$r), del$r), c(rev(del$Lower), del$Upper), col = 'grey80', border = NA)
points(del$r, del$Est, pch = 19, col = 'black')
lines(del$r, del$Est, lty = 'solid', col = 'black')
lines(del$r, del$Lower, lty = 'dashed', col = 'red')
lines(del$r, del$Upper, lty = 'dashed', col = 'red')

par(mfrow = c(1, 1))


## ----nonstatfit1--------------------------------------------------------------

set.seed(7)
n <- 100
r <- 10
x <- rgevr(n, r, loc = 100 + 1:n / 50,  scale = 1 + 1:n / 100, shape = 0)

## Plot the largest order statistic
plot(x[, 1])

## Creating covariates (linear trend first)
covs <- as.data.frame(seq(1, n, 1))
names(covs) <- c("Trend1")
## Create some unrelated covariates
covs$Trend2 <- rnorm(n)
covs$Trend3 <- 30 * runif(n)

## Use full data
fit_full <- gevrFit(data = x, method = "mle", locvars = covs, locform = ~ Trend1 + Trend2*Trend3,
scalevars = covs, scaleform = ~ Trend1)

## Only use r = 1
fit_top_only <- gevrFit(data = x[, 1], method = "mle", locvars = covs, locform = ~ Trend1 + Trend2*Trend3,
scalevars = covs, scaleform = ~ Trend1)


## ----nonstatfit2--------------------------------------------------------------

## Show summary of estimates
fit_full
fit_top_only


## ----nonstatfit3--------------------------------------------------------------

## Compare AIC of three models
fit_reduced1 <- gevrFit(data = x, method = "mle", locvars = covs, locform = ~ Trend1,
scalevars = covs, scaleform = ~ Trend1)

fit_reduced2 <- gevrFit(data = x, method = "mle", locvars = covs, locform = ~ Trend1,
scalevars = covs, scaleform = ~ Trend1, gumbel = TRUE)

AIC(fit_full)
AIC(fit_reduced1)
AIC(fit_reduced2)

LRT <- as.numeric(2 * (logLik(fit_reduced1) - logLik(fit_reduced2)))

pval <- pchisq(LRT, df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)

round(pval, digits = 3)


## ----rfa_example1-------------------------------------------------------------

if (!requireNamespace("SpatialExtremes", quietly=TRUE))
  install.packages("SpatialExtremes", repos = c("http://R-Forge.R-project.org"), deps=TRUE)

set.seed(7)
n.site <- 4
n.obs <- 15
## Simulate a max-stable random field
locations <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
colnames(locations) <- c("lon", "lat")
## Smith model
U <- SpatialExtremes::rmaxstab(n.obs, locations, "gauss", cov11 = 16, cov12 = 0, cov22 = 16)
cor(U, method = "spearman")


## ----rfa_example3-------------------------------------------------------------

## Transform to GEV margins
locations <- c(8, 10, 12, 9)
out <- SpatialExtremes::frech2gev(U, loc = 0, scale = 1, shape = 0.2)
out <- out + t(matrix(rep(locations, nrow(out)), ncol = nrow(out)))
out <- as.vector(out)

## Create design matrix for the location parameters
loc <- cbind.data.frame(as.factor(sort(rep(seq(1, n.site, 1), n.obs))))
colnames(loc) <- c("Site")

z <- gevrFit(out, locvars = loc, locform = ~ Site)
z


