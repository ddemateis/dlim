## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(ggplot2)

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
                 df_l = 10)

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
plot_DLF(new_modifiers = seq(0.1,0.9,0.1), 
         mod_fit = dlim_fit, 
         mod_name = "modifier", 
         plot_by = "time", 
         exposure_time = 10:46,
         time_pts = c(20, 30, 40))

## -----------------------------------------------------------------------------
plot_DLF(new_modifiers = c(0.25, 0.5, 0.75),
         mod_fit = dlim_fit, 
         mod_name = "modifier", 
         plot_by = "modifier",
         exposure_time = 10:46) +
  xlab("months after parturition")

## -----------------------------------------------------------------------------
#predict
dlim_pred <- predict(dlim_fit,
                     newdata = seq(0.1, 0.9, 0.1))

#create data frame for plotting
dlim_cumul_df <- data.frame(Estimate = c(dlim_pred$est_dlim$betas_cumul),
                            LB = c(dlim_pred$est_dlim$cumul_LB),
                            UB = c(dlim_pred$est_dlim$cumul_UB),
                            Modifier = c(dlim_pred$est_dlim$modifiers))

#plotting
ggplot(dlim_cumul_df, aes(x = Modifier, y = Estimate)) +
    geom_point(color = "blue") +
    geom_errorbar(aes(ymin=LB, ymax=UB)) +
    geom_hline(yintercept = 0,  color = "black", size=1) +
    xlab("Modifier") + 
    ylab("Change in response per unit of exposure") +
    ggtitle("Cumulative Effect Esimates") +
    theme_bw()

## -----------------------------------------------------------------------------
#predict
new_mods <- c(0.25, 0.5, 0.75)
dlim_pred <- predict(dlim_fit, 
                     newdata = c(0.25, 0.5, 0.75), 
                     type = "DLF")

#create data frame for plotting
dlim_pred_df <- data.frame(Estimate = c(t(dlim_pred$est_dlim$betas)),
                           LB = c(t(dlim_pred$est_dlim$LB)),
                           UB = c(t(dlim_pred$est_dlim$UB)),
                           Week = rep(1:37,length(new_mods)),
                           Modifier = rep(new_mods, each = 37))

#plotting
ggplot(dlim_pred_df, aes(x = Week, y = Estimate)) +
    geom_point(color = "blue") +
    geom_errorbar(aes(ymin=LB, ymax=UB)) +
    geom_hline(yintercept = 0,  color = "black", size=1) +
    facet_grid(cols = vars(Modifier), labeller = "label_both") + 
    xlab("Exposure week") + 
    ylab("Change in response per unit of exposure") +
    theme_bw()

## -----------------------------------------------------------------------------
model_comparison(fit = dlim_fit, 
                 null = "none",
                 x = exposure, 
                 B = 5)

