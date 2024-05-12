set.seed(123)                       # 乱数シードの固定
n <- 5000                           # サンプルサイズ
e <- rnorm(n, mean = 0, sd = 0.6)   # 誤差項
x <- rnorm(n, mean = 0, sd = 1)     # 説明変数                        
y <- 2 + 3 * x  + e                 # 被説明変数

mc <- 10000     # MCMCの長さ
bn <- 2000      # burn-inの長さ

library(brms)

data <- list(N = n, y = y, x = x)
fit.stan <- brm(y ~ x, data = data, iter = mc, warmup = bn, chain=4)

png(filename = "R/plot.png")
plot(fit.stan)
dev.off()

