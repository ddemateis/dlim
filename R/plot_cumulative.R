#' Plot Distributed Lag Function
#' @description Plot estimated distributed lag function values from a DLIM object, can also compare those of a DLM
#' @export
#' @import ggplot2
#' @import viridis
#' @param new_modifiers a vector of new modifier values for prediction (class "\code{numeric}")
#' @param mod_fit DLIM model object (class "\code{dlim}")
#' @param dlm_fit a list containing a \code{crossbasis} object from the \pkg{dlnm} package as the first element and a DLM model object as the second element (class "\code{list}")
#' @param mod_name modifier name (character)
#' @param mod_trans if modifiers are transformed, specify back transformation function (class "\code{character}")
#' @param link_trans if family for \code{glm} is not Gaussian, specify back transformation to undo link function (class "\code{character}")
#' @return This function returns a ggplot for cumulative effects, including for a DLM if specified

plot_cumulative <- function(new_modifiers, mod_fit, dlm_fit=NULL, mod_name = NULL, mod_trans = NULL, link_trans = NULL){
  library(ggplot2)
  library(viridis)

  #predict DLIM
  model_pred <- predict(object=mod_fit, newdata = new_modifiers, type="CE")

  if(!is.null(dlm_fit)){
    #predict DLM
    cb_dlm <- dlm_fit[[1]]
    model_dlm <- dlm_fit[[2]]
    nlag <- attr(cb_dlm, "lag")[2]+1
    dlm_crosspred <- crosspred(cb_dlm,model_dlm,at=rep(1,nlag),cen = F)
    cumul_betas <- dlm_crosspred$allfit
    z <- qnorm(1 - (1 - 0.95)/2)
    cumul_lb <- dlm_crosspred$allfit - z * dlm_crosspred$allse #for some reason $alllow is gone
    cumul_ub <- dlm_crosspred$allfit + z * dlm_crosspred$allse #for some reason $allhigh is gone
  }

  #back transform modifiers if specified
  if(!is.null(mod_trans)){
    new_modifiers <- do.call(mod_trans,list(new_modifiers))
  }

  #plot
  ref_line <- 0
  if(is.null(dlm_fit)){

    df_cumul <- data.frame(Modifiers = c(new_modifiers),
                           Cumul_Effect = c(model_pred$est_dlim$betas_cumul),
                           LB = c(model_pred$est_dlim$cumul_LB),
                           UB = c(model_pred$est_dlim$cumul_UB))

    if(is.null(link_trans)){
      if(mod_fit$fit$family$family!="gaussian"){
        warning("Family is not Gaussian. Use the link_trans argument to transform estimate.")
      }
    }else{
      ref_line <- do.call(link_trans, list(0))
      df_cumul[,2:4] <- do.call(link_trans,list(df_cumul[,2:4]))
    }

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

    if(is.null(link_trans)){
      if(mod_fit$fit$family$family!="gaussian"){
        warning("Family is not Gaussian. Use the link_trans argument to transform estimate.")
      }
    }else{
      ref_line <- do.call(link_trans, list(0))
      df_cumul[,2:4] <- do.call(link_trans,list(df_cumul[,2:4]))
    }

    ggplot(df_cumul, aes(x=Modifiers,y=Cumul_Effect, color=Model, fill=Model)) +
      geom_hline(yintercept = ref_line) +
      geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.2, color=NA)+
      geom_line()+
      xlab(ifelse(is.null(mod_name), "Modifier", mod_name)) +
      ylab("Cumulative Effect") +
      theme_classic() +
      scale_fill_viridis(discrete=T) +
      scale_color_viridis(discrete=T)

  }

}
