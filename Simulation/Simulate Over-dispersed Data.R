library(data.table)
library(brms)
set.seed(72)
# Generate poisson-generated counts with low and high levels of noise/dispersion
dt <- data.table(predictor = rpois(60,lambda = 15))
hist(dt[,predictor],breaks = 10)
dt[,response:= predictor + round(rnorm(nrow(dt), 0, sd = 1))]
hist(dt[,response], breaks = 10)
dt[,response_overdisp:= abs(predictor+round(rnorm(nrow(dt), 0, sd = 10)))]
hist(dt[,response_overdisp], breaks = 10)
bm0.nb <- brm(response ~ predictor, dt, family = "negbinomial")
bm0.over.nb <- brm(response_overdisp ~ predictor, dt, family = "negbinomial")

model0 <- bm0.nb # different models can be tested
shape_post <- log(posterior_samples(model0)$shape)
means_post <- rowMeans(posterior_predict(model0))
dispersion_post <- 1 + means_post / shape_post
hist(dispersion_post,xlim = c(0.9, max(dispersion_post)))
abline(v = 1, col = "red")
hypothesis(model0, paste("1+", mean(means_post), "/shape > 1", sep=""), class=NULL)
