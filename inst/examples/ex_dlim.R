library(dlim)
data("ex_data")
dlim_fit <- dlim(y = ex_data$y, 
                 x = ex_data$exposure, 
                 modifier = ex_data$modifier, 
                 z = ex_data$z, 
                 df_m = 10, 
                 df_l = 10)
dlim_pred <- predict(dlim_fit, 
                     newdata = 0.5, 
                     type="CE")