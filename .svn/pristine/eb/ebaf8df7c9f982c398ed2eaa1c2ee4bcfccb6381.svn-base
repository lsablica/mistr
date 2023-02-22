# Introduction example
library(mistr)
v = 4
op <- par(mfrow = c(1, 2), mar = c(v,v/2,v/2,2*v))

plot(density(stocks$SAP), xlim = c(-0.07, 0.07), xlab = "", ylab = "", main = "", cex.axis = 0.75)
title(main= "Density of SAP log-returns 
     and normal distribution ", cex.main=0.75)
x <- seq(-0.07, 0.07, 0.001)
lines(x, dnorm(x, mean(stocks$SAP), sd(stocks$SAP)), col = "red")
op <- par(mfrow = c(1, 2), mar = c(v,2*v,v/2,v/2))
qqnorm(stocks$SAP, xlab = "", ylab = "", main = "", cex.axis = 0.75)
qqline(stocks$SAP)
title(main= "Normal Q-Q Plot", cex.main=0.75)

par(op)


# Computational framework
N <- normdist(mean = 1, sd = 3)
N

d(N, c(1, 2, 3))
p(N, c(1, 2, 3))
q(N, c(0.1, 0.2, 0.3))
r(N, 3)

B <- binomdist(size = 12, prob = 0.3)
plim(B, c(-3, 0, 3, 12))
qlim(B, plim(B, c(0, 3, 7, 12)))

# Adding transformation
E <- expdist(2)
E * 2
E^2

E2 <- E * -2
E3 <- E2 * 5
E3

Norm_trafo <- (N - 1)^(1/3)
Norm_trafo

Binom_trafo <- -3 * log(B + 4)
q(Binom_trafo, c(0.05, 0.5, 0.95))
plim(Binom_trafo, c(-6, -5, 0))
sudo_support(Binom_trafo)

par(mai = c(0.4, 0.4, 0.2, 0.2))
plot(Norm_trafo, xlim1 = c(-2.5, 2.5), ylab1 = "", cex.axis = 0.75)

library(ggplot2) # requires ggplot2
autoplot(Norm_trafo, xlim1 = c(-2.5, 2.5))

QQplotgg(Norm_trafo, r(Norm_trafo, 1000), conf = 0.99, ylab = NULL, xlab = NULL)

# Mixtures
mixdist(c("norm", "unif"), list(c(2, 2), c(1, 5)), weights = c(0.5, 0.5))

M <- mixdist(Norm_trafo, Binom_trafo, expdist(0.5), weights = c(0.4, 0.2, 0.4))

DM <- mixdist(3 * binomdist(12, 0.4), -2*poisdist(2) + 12, weights=c(0.5, 0.5))
y <- c(0.05, 0.4, p(-DM, c(-5, -10, -15)), 0.95)
x <- q(-DM, y)
autoplot(-DM, which = "cdf", only_mix = TRUE, xlim1 = c(-37, 0)) +
  annotate("point", x, y, col = "white")

sudo_support(M)

M_trans <- -2 * (M)^(1/3)
r(M_trans, 4)

autoplot(M_trans)

# Composite distributions
C <- compdist(-paretodist(1, 1), normdist(0, 2), geomdist(0.3) + 2, 
              weights = c(0.15, 0.7, 0.15), breakpoints = c(-3, 3),
              break.spec = c("L", "R"))
C

C2 <- compdist(-expdist(2), poisdist(), expdist(2),
               weights = c(0.25, 0.5, 0.25), breakpoints = c(0, 0))
C2

par(mai = c(0.4, 0.4, 0.2, 0.2))
plot(C, xlim1 = c(-15, 15), ylab1 = "", cex.axis = 0.75, mtext_cex = 0.75)

autoplot(C2, text_ylim = 0.01)

C_trans <- -0.5 * (C + 7)

q(C_trans, c(0.075, 0.5, 0.7, 0.9))
r(C_trans, 4)

autoplot(C_trans, xlim1 = c(-10,5))

# Combining mixture and composite distributions
C3 <- compdist(M_trans - 3, C_trans, weights = c(0.5, 0.5), breakpoints = -4.5)
C3_trans <- -2 * C3 + 2

plim(C3_trans, c(6, 10, 12))
qlim(C3_trans, c(0.3, 0.5, 0.7))

autoplot(C3_trans, xlim1 = c(0,20), text_ylim = 0.01, grey = TRUE)

autoplot(mixdist( C3_trans, C2 + 5, weights = c(0.7, 0.3)), xlim1 = c(0, 15))

# Data modeling
PNP_model <- PNP_fit(stocks$SAP)
PNP_model

plot(PNP_model, ylab1 = "", ylab2 = "")

GNG_model <- GNG_fit(stocks$SAP)
GNG_model

autoplot(GNG_model)

risk(GNG_model, c(0.02, 0.05, 0.07, 0.1, 0.2, 0.3))