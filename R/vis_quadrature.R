# Visualise the quadrature rules
library(tidyverse)

sigma <- 0.4
n_int <- 7
# gaussian quadrature
gq_classes <- gaussquad_classes(sigma = 0.4,
                                n_classes = n_int)
even_classes <- quantile2_classes_mean(sigma = 0.4,
                                       n_classes = n_int)

xmax <- max(
  qlnorm(0.99, 0, sigma),
  max(gq_classes$activity),
  max(even_classes$activity)
)

par(mfrow = c(3, 1))

plot(function(x) {dlnorm(x, 0, sigma)},
     xlim = c(0, xmax),
     ylab = "p(v)",
     xlab = "v",
     main = "Target distribution",
     axes = FALSE)
box()
axis(side = 1)
axis(side = 2,
     las = 2)

plot(even_classes$fraction ~ even_classes$activity,
     pch = 16,
     xlim = c(0, xmax),
     ylim = c(0, max(even_classes$fraction) * 1.1),
     ylab = "fraction",
     xlab = "v",
     main = "Even classes",
     axes = FALSE)
box()
axis(side = 1)
axis(side = 2,
     at = c(0, 1/n_int),
     labels = c("0", sprintf("1/%s", n_int)),
     las = 2)
arrows(x0 = even_classes$activity,
       x1 = even_classes$activity,
       y0 = 0,
       y1 = even_classes$fraction,
       length = 0)

plot(gq_classes$fraction ~ gq_classes$activity,
     pch = 16,
     xlim = c(0, xmax),
     ylim = c(0, max(gq_classes$fraction) * 1.1),
     ylab = "fraction",
     xlab = "v",
     main = "Gaussian quadrature",
     axes = FALSE)
box()
axis(side = 1)
axis(side = 2,
     las = 2)
arrows(x0 = gq_classes$activity,
       x1 = gq_classes$activity,
       y0 = 0,
       y1 = gq_classes$fraction,
       length = 0)

