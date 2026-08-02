---
title: "Mid-Atlantic Regression"
keep-md: true
format:
  pdf:
    toc: true
    number-sections: true
    keep-tex: true
    header-includes: |
      \usepackage{hyperref}
      \hypersetup{
        pdfpagemode=UseOutlines, % open bookmarks panel
        bookmarksopen=true,      % unfold top-level bookmarks
        bookmarksopenlevel=10,   % unfold all levels
        bookmarksnumbered=true   % number bookmarks
      }
---



---

__Code version: 0.0.18__









# Data





# Angular Distance Between Observations

::: {.cell}

```{.r .cell-code}
Y <- fulldf %>%
  select(starts_with("Y")) %>%
  as.matrix()
ang_dist <- acos(pmin(Y %*% t(Y), 1)) # pmin clamps dot products to <=1, preventing NaN from acos due to floating-point rounding
tibble::enframe(ang_dist[lower.tri(ang_dist)], value = "ang_dist") %>%
  ggplot() +
  geom_histogram(aes(x = ang_dist), bins = 30) +
  geom_rug(aes(x = ang_dist)) +
  scale_x_continuous(name = "Angular Distance (radians)")
```

::: {.cell-output-display}
![](reproduce_midatlantic_files/figure-pdf/unnamed-chunk-1-1.pdf){fig-pos='H'}
:::
:::


# Euc Covariates

::: {.cell}

```{.r .cell-code}
xe <- fulldf %>% 
  select(westedge)
xestd <- xe %>% scale() %>% as_tibble()
```
:::


# Our SvMF Model




# Pretty Data Plot

::: {.cell}

```{.r .cell-code}
ymean <- colMeans(fulldf %>% select(starts_with("Y")))
ymean <- ymean/sqrt(sum(ymean^2))
#project to equator
ymeanproj <- c(ymean[1:2], 0)
ymeanproj/vnorm(ymeanproj)
```

::: {.cell-output .cell-output-stdout}

```
        Y1         Y2            
 0.8791608 -0.4765251  0.0000000 
```


:::

```{.r .cell-code}
target <- c(1,0,0)
prettyrot <- parallel_transport_mat(ymeanproj, target)
prettyY <- (fulldf %>%
  select(Y1, Y2, Y3) %>%
  as.matrix()) %*% t(prettyrot)
colnames(prettyY) <- paste0("Y", 1:3)
prettyX <- (fulldf %>%
  select(X1, X2, X3) %>%
  as.matrix()) %*% t(prettyrot)
colnames(prettyX) <- paste0("X", 1:3)
prettyXmean <- colMeans(prettyX)
prettyXmean <- prettyXmean/vnorm(prettyXmean)
# denomloc <- drop(prettyrot %*% as_mobius_link_cann(mod$mean)$Qs[,1])
plotdata <- bind_cols(prettyY, prettyX, westedge = fulldf$westedge) %>%
  as_tibble() %>%
  ggplot() +
  geom_segment(aes(x=X2, y=X3, xend=Y2, yend=Y3), lty = "dashed", col = "grey") +
  geom_point(aes(x=Y2, y=Y3), shape = 4) +
  geom_point(aes(x=X2, y=X3, fill = westedge), show.legend = FALSE, shape = 21) +
  geom_path(data = circle_df, aes(x = x, y = y), inherit.aes = FALSE, color = "grey") +
  scale_fill_manual(values = c("black", "white")) +
  theme_void() +
  coord_fixed() 
  # mark mean of spherical covariate (obtained earlier)
  # annotate("point", x = prettyXmean[2], y = prettyXmean[3], shape = 24, size = 2)
  # annotate("point", x = denomloc[2], y = denomloc[3], shape = 3, size = 4, col = "blue")
plotdata
```

::: {.cell-output-display}
![](reproduce_midatlantic_files/figure-pdf/unnamed-chunk-3-1.pdf){fig-pos='H'}
:::
:::




# Predictions
The below 'ours' model is using LinEuc's link with an extra covariate for an intercept and G01 free.
The results here are from a local optimisation using gradient from default starting values.
I've checked that global search for the best vMF regression does not find anything better.


::: {.cell}

```{.r .cell-code}
cann <- as_mobius_link_cann(mod$mean)
cann
```

::: {.cell-output .cell-output-stdout}

```
$P
         [,1]       [,2]        [,3]
Y1  0.4913098  0.5626852 -0.66483080
Y2 -0.4803628 -0.4616721 -0.74572813
Y3 -0.7265440  0.6857436  0.04346904

$Bs
          [,1]      [,2]
[1,] 0.9029649 0.0000000
[2,] 0.0000000 0.7702704

$Qs
         [,1]       [,2]       [,3]
X1  0.2657933  0.5868954 -0.7647926
X2 -0.5034120 -0.5920730 -0.6293059
X3 -0.8221498  0.5522710  0.1380812

$Be
           [,1]      [,2]
[1,] 0.07230452 0.0000000
[2,] 0.00000000 0.2525752

$Qe
          [,1]       [,2]      [,3]
dummyzero    1  0.0000000 0.0000000
westedge     0  0.1254316 0.9921023
ones         0 -0.9921023 0.1254316

$ce
[1] 1

attr(,"class")
[1] "mobius_link_cann" "list"            
```


:::
:::


Lets try to interpret the fitted link.


::: {.cell}

```{.r .cell-code}
cann$Bs %*% t(cann$Qs[,-1])
```

::: {.cell-output .cell-output-stdout}

```
             X1         X2        X3
[1,]  0.5299460 -0.5346212 0.4986814
[2,] -0.5890971 -0.4847357 0.1063598
```


:::

```{.r .cell-code}
cann$Be %*% t(cann$Qe[,-1])
```

::: {.cell-output .cell-output-stdout}

```
     dummyzero    westedge        ones
[1,]         0 0.009069271 -0.07173348
[2,]         0 0.250580391  0.03168090
```


:::
:::


The first direction away from $B_{01}$ (first row in above) is roughly equally influenced by X1, X2 and X3 (all are about 0.5) with west edge having very little influence (given the values of standardised westedge).
The second direction away from $B_{01}$ (second row in above) is roughly equally influenced by X1 and X2 but much less by X3 (which is the N-S direction) and westedge plays a role.


::: {.cell}

```{.r .cell-code}
cann$Qs[,1]
```

::: {.cell-output .cell-output-stdout}

```
        X1         X2         X3 
 0.2657933 -0.5034120 -0.8221498 
```


:::
:::


Some general scaling occurs with greater influence from X3 (N-S direction), then X2 then X1.

## Table of Estimates


```{.r .cell-code}
df <- cbind("diag($B_s$)" = diag(cann$Bs),
      "diag($B_e$)" = diag(cann$Be),
      tQe = t(cann$Qe[-1,-1])) %>%
  as.data.frame()
row.names(df) <- c("$t_2$", "$t_3$")
mykbl <- kbl(df, booktabs = TRUE, position = "!h", escape = FALSE, format = "latex", digits = 2) %>%
  add_header_above(c(" "=3, "$R_e^{\\\\top}$" = 2), escape = FALSE)
mykbl
```


\begin{tabular}[t]{lrrrr}
\toprule
\multicolumn{3}{c}{ } & \multicolumn{2}{c}{$R_e^{\top}$} \\
\cmidrule(l{3pt}r{3pt}){4-5}
  & diag($B_s$) & diag($B_e$) & westedge & ones\\
\midrule
$t_2$ & 0.90 & 0.07 & 0.13 & -0.99\\
$t_3$ & 0.77 & 0.25 & 0.99 & 0.13\\
\bottomrule
\end{tabular}


## Plot











Note that we have not optimised Q for ESAG2 because optimisation did not converge.


::: {.cell}
::: {.cell-output-display}
![](reproduce_midatlantic_files/figure-pdf/unnamed-chunk-7-1.pdf)
:::
:::


# Angular Distance Between Predictions

::: {.cell}

```{.r .cell-code}
pred_ang_dist <- acos(pmin(mod$pred %*% t(mod$pred), 1))
tibble::enframe(pred_ang_dist[lower.tri(pred_ang_dist)], value = "ang_dist") %>%
  ggplot() +
  geom_histogram(aes(x = ang_dist), bins = 30) +
  geom_rug(aes(x = ang_dist)) +
  scale_x_continuous(name = "Angular Distance (radians)")
```

::: {.cell-output-display}
![](reproduce_midatlantic_files/figure-pdf/unnamed-chunk-8-1.pdf){fig-pos='H'}
:::
:::


# LOOCV MSE

::: {.cell}

```{.r .cell-code}
loocvmseSvMF <- function(mod){
  stopifnot(inherits(mod$mean, "mobius_link_Omega"))
  dists <- pbapply::pblapply(1:nrow(mod$y), function(idx){
    newmod <- mobius_SvMF(mod$y[-idx,],
                xs = mod$xs[-idx,],
                xe = mod$xe[-idx,c(-1,-ncol(mod$xe)), drop = FALSE],
                fix_qs1 = FALSE,
                type = "LinEuc",
                G01behaviour = "free",
                mean = mod$mean,
                k = mod$k,
                a = mod$a,
                G0 = mod$G0)
    pred <- mobius_link(xs = mod$xs[idx,, drop = FALSE],
                   xe = mod$xe[idx,, drop = FALSE],
                   param = newmod$mean)
    obs <- mod$y[idx,]
    Euc <- vnorm(drop(obs - pred))
    angle <- acos(rowSums(obs * pred))
    return(c(
      Euc = Euc,
      angle = angle
    ))
  })
  dists <- dists %>%
    simplify2array() %>%
    t() %>%
    as_tibble()
  dists %>%
    summarise(across(everything(), ~sum(.x^2)/nrow(mod$y)))
}
loocvmseSvMF(mod)
```

::: {.cell-output .cell-output-stdout}

```
# A tibble: 1 x 2
     Euc  angle
   <dbl>  <dbl>
1 0.0107 0.0107
```


:::
:::


Both of these are smaller than the LOOCV MSE that Rosenthal's PLT acheived of 0.074, which corresponds to the Euc metric MSE.

# Residual Size
Geodesic distance between predicted mean and observation

::: {.cell messages='false'}
::: {.cell-output-display}
![](reproduce_midatlantic_files/figure-pdf/unnamed-chunk-10-1.pdf)
:::
:::


# DoF

::: {.cell}
::: {.cell-output .cell-output-stdout}

```
 Ours  IAG1 ESAG1 ESAG2 
   16    17    22    25 
```


:::
:::


# AIC

::: {.cell}
::: {.cell-output .cell-output-stdout}

```
     Ours      IAG1     ESAG1     ESAG2 
-463.2407 -370.7254 -459.5946 -436.2852 
```


:::
:::


# Main Figure

::: {.cell}

```{.r .cell-code}
(plotdata + ggtitle("(a) data") + theme(plot.title = element_text(hjust = 0.5))) +
(plotours + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ggtitle("(b) ours")) +
  (plotESAG1 + ggtitle("(c) structure 1") + theme(plot.title = element_text(hjust = 0.5))) +
  (plotESAG2 + ggtitle("(d) structure 2") + theme(plot.title = element_text(hjust = 0.5))) +
  plot_layout(ncol = 4, widths = c(3,3,3,3))
```

::: {.cell-output-display}
![](reproduce_midatlantic_files/figure-pdf/unnamed-chunk-13-1.pdf){fig-pos='H'}
:::

```{.r .cell-code}
ggsave("midatlantic_fig.pdf", width = 12, height = 3.5)
```
:::


Caption: Regression for the midatlantic ridge data.
From left: midatlanic ridge (circles) and corresponding locations on the continent (crosses) from Rosenthal el at (2014);
our regression; Paine et al structure 1 regression;
Paine et al structure 2 regression.
The sphere is shown orthogonally projected with north pointing up the page.
Arrows: start at the predicted mean and end at the observed continental location, and thus represent residuals.
Filled symbols: eastern side.
Unfilled symbols: western side.
Triangles: Mean of ridge locations and corresponding predicted mean for western or eastern side.
Plus symbol: estimated value of $r_{s1}$.
Pair of black lines: direction of the estimated second (solid) and third (dashed) columns of $B_0$ located at the estimated first column $b_{01}$ of $B_0$.
Pair of grey lines: estimated directions of $\gamma_{02}$ (solid) and $\gamma_{03}$ (dashed) located at the estimated $\gamma_{01}$.




# Appendix

## Hessian of Likelihood at Optimum
The parameter vector is longer than the DoF because of the constraints in the optimisation.
There should be DoF + `(3-1) * (3 - 2) / 2` positive eigenvalues.
The term `(3-1) * (3 - 2) / 2` if for the commutativity constraint on Omega, which the likelihood computation does not account for, so appears as extra degrees of freedom in the Hessian of the likelihood.


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
 [1] 8.254073e+03 6.672951e+03 2.990692e+03 2.470448e+03 1.236519e+03
 [6] 3.582952e+02 9.508066e+01 2.869560e+01 2.125170e+01 9.630614e+00
[11] 5.113127e+00 4.746740e+00 2.960021e+00 1.951306e+00 1.054189e+00
[16] 4.074127e-01 2.936479e-06
```


:::
:::


These are all positive and non-zero, which confirms that the optimisation routine has found a local maximum of the likelihood.

## Check Other Starts for vMF

::: {.cell}

```{.r .cell-code}
restarts <- pbapply::pblapply(1:100, function(seed){
  start <- rand_mobius_link_cann(p = 3, qs = 3, qe = ncol(xe) + 2, preseed = seed)
  # convert to LinEuc form:
  set.seed(seed+1)
  Qe <- mclust::randomOrthogonalMatrix(ncol(xe)+1, 3-1)
  bigQe <- cbind(0, rbind(0, Qe))
  bigQe[, 1] <- 0
  bigQe[1,1] <- 1
  start$Qe <- bigQe
  start$ce <- 1
  modvMF <- mobius_vMF(y = fulldf %>% select(starts_with("Y")) %>% as.matrix(),
      xs = fulldf %>% select(starts_with("X")) %>% as.matrix(),
      xe = xe %>% as.matrix(),
      fix_qs1 = FALSE, type = "LinEuc",
      start = start)
}, cl = 2)
lapply(restarts, "[[", "obj") %>%
  unlist() %>%
  enframe("seed", "objective") %>%
  ggplot()+
  geom_freqpoly(aes(x = objective), bins = 30) +
  geom_vline(xintercept = mod$preest$nlopt$objective, col = "blue") +
  geom_rug(aes(x = objective))
```

::: {.cell-output-display}
![](reproduce_midatlantic_files/figure-pdf/unnamed-chunk-15-1.pdf){fig-pos='H'}
:::
:::


Other initial parameters have not improved on the default initial parameters.

## Automated multistart for vMF

::: {.cell}

```{.r .cell-code}
mod_vMF_multistart <- mobius_vMF_multistart(
  y  = fulldf %>% select(starts_with("Y")) %>% as.matrix(),
  xs = fulldf %>% select(starts_with("X")) %>% as.matrix(),
  xe = xe %>% as.matrix(),
  fix_qs1 = FALSE, type = "LinEuc"
)
cat("Default vMF obj:    ", mod$preest$nlopt$objective, "\n",
    "Multistart vMF obj: ", mod_vMF_multistart$nlopt$objective, "\n", sep = "")
```

::: {.cell-output .cell-output-stdout}

```
Default vMF obj:    -0.9960845
Multistart vMF obj: -0.9960844
```


:::
:::


## Check Other Starts for SvMF

::: {.cell}

```{.r .cell-code}
restarts <- pbapply::pblapply(1:100, function(seed){
  # randomly generates a SpEuc-form link
  start <- rand_mobius_link_cann(p = 3, qs = 3, qe = ncol(xestd) + 2, preseed = seed)
  # convert to LinEuc form:
  set.seed(seed+1)
  Qe <- mclust::randomOrthogonalMatrix(ncol(xestd)+1, 5-1)
  bigQe <- cbind(0, rbind(0, Qe))
  bigQe[, 1] <- 0
  bigQe[1,1] <- 1
  start$Qe <- bigQe
  start$ce <- 1
  mobius_SvMF(y = fulldf %>% select(starts_with("Y")) %>% as.matrix(),
              xs = fulldf %>% select(starts_with("X")) %>% as.matrix(), 
              xe = xestd %>% as.matrix(),
              type = "LinEuc",
              G01behaviour = "free",
              mean = start)
}, cl = 2)
badrestarts <- unlist(lapply(restarts, inherits, "try-error"))
restarts <- restarts[!badrestarts]
lapply(restarts, "[[", "AIC") %>%
  unlist() %>%
  tibble::enframe("seed", "AIC") %>%
  ggplot()+
  geom_freqpoly(aes(x = AIC), bins = 30) +
  geom_vline(xintercept = mod$AIC, col = "blue") +
  geom_rug(aes(x = AIC))
```

::: {.cell-output-display}
![](reproduce_midatlantic_files/figure-pdf/unnamed-chunk-17-1.pdf){fig-pos='H'}
:::
:::


Other initial parameters have not improved on the default initial parameters.

## Automated multistart for SvMF

::: {.cell}

```{.r .cell-code}
mod_SvMF_multistart <- mobius_SvMF_multistart(
  y  = fulldf %>% select(starts_with("Y")) %>% as.matrix(),
  xs = fulldf %>% select(starts_with("X")) %>% as.matrix(),
  xe = xestd %>% as.matrix(),
  type = "LinEuc",
  G01behaviour = "free"
)
cat("Default SvMF AIC:    ", mod$AIC, "\n",
    "Multistart SvMF AIC: ", mod_SvMF_multistart$AIC, "\n", sep = "")
```

::: {.cell-output .cell-output-stdout}

```
Default SvMF AIC:    -463.2407
Multistart SvMF AIC: -463.2407
```


:::
:::


## Likelihood Ratio Test of vMF vs SvMF

First get best vMF model

::: {.cell}

```{.r .cell-code}
mod_vMF <- mobius_vMF(y = mod$y,
                      xs = mod$xs,
                      xe = mod$xe,
                      fix_qs1 = FALSE,
                      type = "LinEuc")
mod_vMF_restarts <- pbapply::pblapply(1:100, function(seed){mobius_vMF_refit(mod_vMF, seed)})
```

::: {.cell-output .cell-output-stderr}

```
Warning in mobius_vMF(y = mod_vMF$y, xs = mod_vMF$xs, xe = mod_vMF$xe, fix_qs1
= mod_vMF$linktype$fix_qs1, : NLOPT_MAXEVAL_REACHED: Optimization stopped
because maxeval (above) was reached.
```


:::

::: {.cell-output .cell-output-stdout}

```
Error in mobius_link_Omega_check(projectedom) : 
  The following checks failed. Omega_comm: 0.0043
```


:::

```{.r .cell-code}
lapply(mod_vMF_restarts, "[[", "AIC") %>%
  unlist() %>%
  tibble::enframe("seed", "AIC") %>%
  ggplot()+
  geom_freqpoly(aes(x = AIC), bins = 30) +
  geom_vline(xintercept = mod_vMF$AIC, col = "blue") +
  geom_rug(aes(x = AIC))
```

::: {.cell-output-display}
![](reproduce_midatlantic_files/figure-pdf/unnamed-chunk-19-1.pdf){fig-pos='H'}
:::

```{.r .cell-code}
idx <- which.min(lapply(mod_vMF_restarts, "[[", "AIC") %>% unlist())
mod_vMF <- mod_vMF_restarts[[idx]]
mod_vMF$k
```

::: {.cell-output .cell-output-stdout}

```
[1] 255.3967
```


:::

```{.r .cell-code}
mod_vMF$AIC
```

::: {.cell-output .cell-output-stdout}

```
[1] -354.6918
```


:::
:::


Below is plot for comparing likelihoods for give simulated response `y`.


::: {.cell}

```{.r .cell-code}
getlikR <- function(y){
  mod1 <- mobius_vMF(y,
                     xs = mod$xs,
                     xe = mod_vMF$xe,
                     type = "LinEuc",
                     fix_qs1 = FALSE,
                     start = mod_vMF$est)
  mod2 <- mobius_SvMF(y,
                      xs = mod$xs,
                      xe = mod_vMF$xe,
                      type = "LinEuc",
                      fix_qs1 = FALSE,
                      G01behaviour = "free",
                      mean = mod$mean,
                      k = mod$k,
                      a = mod$a,
                      G0 = mod$G0)
  if ((mod1$nlopt$status != 4) || (mod2$nlopt$status != 4)){
    return(list(likR = NA_real_,
         vMF = mod1,
         SvMF = mod2))
  }
  lLik1 <- mobius_SvMF_log_lik(mod1$y, xs = mod1$xs, xe = mod1$xe,
              mean = mod1$est,
              k = mod1$k,
              a = rep(1, 3),
              G0 = mod2$G0) %>%
    colSums()
  lLik2 <- mobius_SvMF_log_lik(mod2$y, xs = mod2$xs, xe = mod2$xe,
              mean = mod2$mean,
              k = mod2$k,
              a = mod2$a,
              G0 = mod2$G0) %>%
    colSums()

  list(likR = -2* (lLik1[["R"]] - lLik2[["R"]]),
       vMF = mod1,
       SvMF = mod2)
}
```
:::


Now we simulate the response 1000 times from the best vMF model and compare the likelihood of the vMF model to the likelihood of the SvMF model on each simulated data.


::: {.cell}

```{.r .cell-code}
null_likRs <- pbapply::pblapply(1:1000, function(seed){
  set.seed(seed)
  y_ld <- rmobius_SvMF(xs = mod$xs,
                   xe = mod_vMF$xe,
                   mean = mod_vMF$mean,
                   k = mod_vMF$k,
                   a = rep(1, 3),
                   G0 = diag(1, 3))
  y <- y_ld[,-4]
  getlikR(y)
}, cl = 2)
obs <- getlikR(mod$y)
```
:::

::: {.cell}

```{.r .cell-code}
null_likRs_vec <- lapply(null_likRs, "[[", "likR") %>% unlist()
sum(is.na(null_likRs_vec))
```

::: {.cell-output .cell-output-stdout}

```
[1] 448
```


:::

```{.r .cell-code}
sum(!is.na(null_likRs_vec))
```

::: {.cell-output .cell-output-stdout}

```
[1] 552
```


:::
:::


A substantial fraction of the simulated fits stop at the evaluation limit rather than the convergence tolerance, and these have been discarded.


::: {.cell}

```{.r .cell-code}
obs$likR
```

::: {.cell-output .cell-output-stdout}

```
[1] 116.5489
```


:::

```{.r .cell-code}
range(null_likRs_vec, na.rm = TRUE)
```

::: {.cell-output .cell-output-stdout}

```
[1]  0.2812075 23.1321965
```


:::
:::

::: {.cell}

```{.r .cell-code}
tibble::enframe(null_likRs_vec, "seed", "likR") %>%
  ggplot() +
  geom_freqpoly(aes(x = likR), bins = 30) +
  geom_vline(xintercept = obs$likR, col = "blue") +
  geom_rug(aes(x = likR))
```

::: {.cell-output .cell-output-stderr}

```
Warning: Removed 448 rows containing non-finite outside the scale range
(`stat_bin()`).
```


:::

::: {.cell-output .cell-output-stderr}

```
Warning: Removed 448 rows containing missing values or values outside the scale range
(`geom_rug()`).
```


:::

::: {.cell-output-display}
![](reproduce_midatlantic_files/figure-pdf/unnamed-chunk-24-1.pdf){fig-pos='H'}
:::
:::


p-value is the probability, under the null, of getting a likelihood-ratio that is at least as large as the observed ratio:

::: {.cell}

```{.r .cell-code}
mean(null_likRs_vec > obs$likR, na.rm = TRUE)
```

::: {.cell-output .cell-output-stdout}

```
[1] 0
```


:::
:::


The observed likelihood ratio is far larger than any value obtained under the null, so this small $p$-value suggests that the data was not drawn from a vMF regression.

