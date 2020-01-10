library(factiv)
set.seed(2140)
cc.prob <- 0.4
cprobs <- c(cc.prob, rep((1-cc.prob)/8, 8))
N <- 1000

c1 <- t(rmultinom(N, 1, cprobs))
colnames(c1) <- c("cc", "ca", "cn", "ac", "aa", "an", "nc", "na", "nn")

z1 <- sample(rep(c(0,1), N/2))
z2 <- sample(rep(c(0,1), N/2))

d1 <- rowSums(c1[,c("cc", "ca", "cn")])*z1 + rowSums(c1[,c("aa", "ac", "an")])
d2 <- rowSums(c1[,c("cc", "ac", "nc")])*z2 + rowSums(c1[,c("aa", "ca", "na")])


beta.0 <- rnorm(N, rowSums(t(c(1, 0.2, 0.1, -0.1, 0.2, 0.2, -0.4, 0.2, 0.7)*t(c1))), 1)
beta.d1 <- rnorm(N, rowSums(t(c(1, 0.2, 0.1, -0.1, 0.2, 0.2, -0.4, 0.2, 0.7)*t(c1))), 1)
beta.d2 <- rnorm(N, rowSums(t(c(1, 0.2, 0.1, -0.1, 0.2, 0.2, -0.4, 0.2, 0.7)*t(c1))), 1)
beta.int <- rnorm(N, rowSums(t(c(1, 0.2, 0.1, -0.1, 0.2, 0.2, -0.4, 0.2, 0.7)*t(c1))), 1)
epsilon <- rnorm(N,0,.5)

y00 <- beta.0 + epsilon
y10 <- beta.0 + beta.d1 + epsilon
y01 <- beta.0 + beta.d2 + epsilon
y11 <- beta.0 + beta.d1 + beta.d2 + beta.int + epsilon

y <- y00 + d1 * (y10-y00) + d2 * (y01-y00) + d1 * d2 * ((y11 - y01) - (y10 - y00))

mydata <- data.frame(y, d1, d2, z1, z2)

out <- iv_factorial(y ~ d1 + d2 | z1 + z2, data = mydata)
