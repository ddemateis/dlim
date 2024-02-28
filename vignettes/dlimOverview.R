## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(dlim)

## -----------------------------------------------------------------------------
data("ex_data")
str(ex_data)

## -----------------------------------------------------------------------------
dlim_fit <- dlim(y = ex_data$y, 
                 x = ex_data$exposure, 
                 modifier = ex_data$modifier, 
                 z = ex_data$z, 
                 df_m = 10, 
                 df_l = 10, 
                 method = "REML")

## -----------------------------------------------------------------------------
dlim_fit

## -----------------------------------------------------------------------------
dlim_pred <- predict(dlim_fit, 
                     newdata = 0.5, 
                     type="CE")
data.frame(cumul_betas = c(dlim_pred$est_dlim$betas_cumul),
           LB = c(dlim_pred$est_dlim$cumul_LB),
           UB = c(dlim_pred$est_dlim$cumul_UB))

## -----------------------------------------------------------------------------
dlim_pred <- predict(dlim_fit, 
                     newdata = 0.5, 
                     type="DLF")
data.frame(betas = c(dlim_pred$est_dlim$betas),
           LB = c(dlim_pred$est_dlim$LB),
           UB = c(dlim_pred$est_dlim$UB))

## -----------------------------------------------------------------------------
plot_cumulative(new_modifiers = seq(0.1,0.9,0.1), 
                mod_fit = dlim_fit, 
                mod_name = "modifier")

## -----------------------------------------------------------------------------
plot_DLF(new_modifiers = seq(0.1,0.9,0.1), 
         mod_fit = dlim_fit, 
         mod_name = "modifier", 
         plot_by = "time", 
         time_pts = c(10,20,30))

## -----------------------------------------------------------------------------
plot_DLF(new_modifiers = c(0.25, 0.5, 0.75),
         mod_fit = dlim_fit, 
         mod_name = "modifier", 
         plot_by = "modifier")

## -----------------------------------------------------------------------------
model_comparison(fit = dlim_fit, 
                 null = "DLM",
                 x = exposure, 
                 B = 5)

