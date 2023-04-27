#' Plot Distributed Lag Function
#' @description Plot estimated distributed lag function values from a DLIM object, can also compare those of a DLM
#' @export
#' @import ggplot2
#' @param new_modifiers a vector of new modifier values for prediction (numeric)
#' @param mod_fit an object of class dlim (dlim)
#' @param dlm_fit a list containing a \code{crossbasis} object from the \pkg{dlnm} package as the first element and a DLM model object as the second element (list)
#' @param trans_fn if modifiers are transformed, specify back transformation function (character)
#' @param mod_name modifier name (character)
#' @return This function returns ggplot

plot_cumulative <- function(new_modifiers, mod_fit, dlm_fit=NULL, trans_fn = NULL, mod_name = NULL){
  library(ggplot2)
  library(viridis)

  #predict DLIM
  model_pred <- predict(object=mod_fit, newdata = new_modifiers, type="CE")

  if(!is.null(dlm_fit)){
    #predict DLM
    cb_dlm <- dlm_fit[[1]]
    model_dlm <- dlm_fit[[2]]
    dlm_crosspred <- crosspred(cb_dlm,model_dlm,at=rep(1,37),cen = F)
    cumul_betas <- dlm_crosspred$allfit
    cumul_lb <- dlm_crosspred$alllow
    cumul_ub <- dlm_crosspred$allhigh
  }

  #back transform modifiers if specified
  if(!is.null(trans_fn)){
    new_modifiers <- do.call(trans_fn,list(new_modifiers))
  }

  #plot
  if(is.null(dlm_fit)){
    df_cumul <- data.frame(Modifiers = c(new_modifiers),
                           Cumul_Effect = c(model_pred$est_dlim$betas_cumul),
                           LB = c(model_pred$est_dlim$cumul_LB),
                           UB = c(model_pred$est_dlim$cumul_UB))

    ggplot(df_cumul, aes(x=Modifiers,y=Cumul_Effect)) +
      geom_hline(yintercept = 0) +
      geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.5 , fill = "grey70")+
      geom_line()+
      xlab(ifelse(is.null(mod_name), "Modifier", mod_name)) +
      ylab("Cumulative Effect") +
      theme_classic() +
      scale_fill_viridis(discrete=T) +
      scale_color_viridis(discrete=T)
  }else{
    model_name <- paste0("DLIM(", mod_fit$cb$df_m,",",mod_fit$cb$df_l,")")
    df_cumul <- data.frame(Modifiers = c(new_modifiers, new_modifiers),
                           Cumul_Effect = c(model_pred$est_dlim$betas_cumul, rep(cumul_betas, length(new_modifiers))),
                           LB = c(model_pred$est_dlim$cumul_LB, rep(cumul_lb, length(new_modifiers))),
                           UB = c(model_pred$est_dlim$cumul_UB, rep(cumul_ub, length(new_modifiers))),
                           Model = factor(c(rep(model_name, length(new_modifiers)), rep("DLM", length(new_modifiers))), levels = c(model_name,"DLM"))
    )

    ggplot(df_cumul, aes(x=Modifiers,y=Cumul_Effect, color=Model, fill=Model)) +
      geom_hline(yintercept = 0) +
      geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.2, color=NA)+
      geom_line()+
      xlab(ifelse(is.null(mod_name), "Modifier", mod_name)) +
      ylab("Cumulative Effect") +
      theme_classic() +
      ylim(-0.1,0.1) +
      scale_fill_viridis(discrete=T) +
      scale_color_viridis(discrete=T)
  }

}
