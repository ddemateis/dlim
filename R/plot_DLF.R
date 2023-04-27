#' Plot Cumulative Effects
#' @description Plot estimated cumulative effects from a DLIM object, can also compare estimated cumulative effects between a DLM and DLIM
#' @export
#' @import ggplot2
#' @import reshape2
#' @param new_modifiers a vector of new modifier values for prediction (numeric)
#' @param mod_fit an object of class dlim (dlim)
#' @param mod_name modifier name that follows variable name nomenclature (character)
#' @param dlm_fit a list containing a \code{crossbasis} object from the \pkg{dlnm} package as the first element and a DLM model object as the second element (list)
#' @param plot_by choose to create plots for particular modifier values, "modifier", or particular time points, "time", (character)
#' @param time_pts a set of time points if plotting by time (numeric)
#' @param trans_fn if modifiers are transformed, specify back transformation function (character)
#' @return This function returns ggplot

plot_DLF <- function(new_modifiers, mod_fit, mod_name, dlm_fit=NULL, plot_by, time_pts=NULL, trans_fn = NULL){
  library(ggplot2)
  library(reshape2)
  library(viridis)

  #predict DLIM
  model_pred <- predict(object=mod_fit, newdata = new_modifiers, type="DLF")

  if(!is.null(dlm_fit)){
    #predict DLM
    cb_dlm <- dlm_fit[[1]]
    model_dlm <- dlm_fit[[2]]
    dlm_crosspred <- crosspred(cb_dlm,model_dlm,at=rep(1,37),cen = F)
    beta_by_lag <- dlm_crosspred$matfit
    dlm_lb <- dlm_crosspred$matlow
    dlm_ub <- dlm_crosspred$mathigh
  }

  #back transform if specified
  if(!is.null(trans_fn)){
    new_modifiers <- do.call(trans_fn,list(new_modifiers))
  }

  if(is.null(time_pts)){
    time_pts <- 1:ncol(model_pred$est_dlim$betas)
  }

  est_dlf <- model_pred$est_dlim$betas[,time_pts] #modifiers x time_pts
  colnames(est_dlf) <- time_pts
  rownames(est_dlf) <- new_modifiers

  est_lb <- model_pred$est_dlim$LB[,time_pts] #modifiers x time_pts
  colnames(est_lb) <- time_pts
  rownames(est_lb) <- new_modifiers

  est_ub <- model_pred$est_dlim$UB[,time_pts] #modifiers x time_pts
  colnames(est_ub) <- time_pts
  rownames(est_ub) <- new_modifiers

  dlf_df <- melt(est_dlf)
  colnames(dlf_df) <-c("Modifiers", "Week", "Effect")

  lb_df <- melt(est_lb)
  colnames(lb_df) <-c("Modifiers", "Week", "Effect")

  ub_df <- melt(est_ub)
  colnames(ub_df) <-c("Modifiers", "Week", "Effect")

  df <- data.frame(dlf_df, LB=lb_df$Effect, UB=ub_df$Effect)

  if(plot_by=="modifier"){
    df$Modifiers <- signif(df$Modifiers,3)
  }else if(plot_by=="week"){
    df$Week <- paste("Week", df$Week)
  }


  m <- length(new_modifiers)

  if(!is.null(dlm_fit)){ #DLM and DLIM
    model_name <- paste0("DLIM(", mod_fit$cb$df_m,",",mod_fit$cb$df_l,")")
    df_time_pts <- data.frame(Modifiers = c(df$Modifiers,df$Modifiers),
                           Week = c(df$Week, df$Week),
                           Effect = c(df$Effect, rep(beta_by_lag[time_pts],each=m)),
                           LB = c(df$LB, rep(dlm_lb[time_pts],each=m)),
                           UB = c(df$UB, rep(dlm_ub[time_pts], each=m)),
                           Model = factor(c(rep(model_name, length(df$Modifiers)), rep("DLM", length(df$Modifiers))))
    )

    if(plot_by=="modifier"){
      colnames(df_time_pts)[which(colnames(df_time_pts)=="Modifiers")] <- mod_name
      ggplot(df_time_pts, aes(x=Week,y=Effect, color=Model, fill = Model)) +
        geom_hline(yintercept = 0)+
        geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.2, color=NA)+
        geom_line()+
        facet_wrap(mod_name, ncol = 3, labeller = label_both) +
        xlab("Time") +
        ylab("Effect") +
        theme_classic() +
        scale_fill_viridis(discrete=T) +
        scale_color_viridis(discrete=T)

    }else if(plot_by=="time"){
      ggplot(df_time_pts, aes(x=Modifiers,y=Effect, color=Model, fill=Model)) +
        geom_hline(yintercept = 0) +
        geom_ribbon(aes(ymin=LB, ymax=UB),alpha=0.2, color=F)+
        geom_line()+
        facet_wrap(vars(Week), ncol = 3, labeller = label_both) +
        xlab(ifelse(is.null(mod_name), "Modifier", mod_name)) +
        ylab("Effect") +
        theme_classic() +
        scale_fill_viridis(discrete=T) +
        scale_color_viridis(discrete=T)

    }
  }else{ #just DLIM
    df_time_pts <- data.frame(Modifiers = c(df$Modifiers),
                           Week = c(df$Week),
                           Effect = c(df$Effect),
                           LB = c(df$LB),
                           UB = c(df$UB)
    )

    if(plot_by=="modifier"){
      colnames(df_time_pts)[which(colnames(df_time_pts)=="Modifiers")] <- mod_name
      ggplot(df_time_pts, aes(x=Week,y=Effect)) +
        geom_hline(yintercept = 0)+
        geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.5, color=F)+
        geom_line()+
        facet_wrap(mod_name, ncol = 3, labeller = label_both) +
        xlab("Time") +
        ylab("Effect") +
        theme_classic() +
        scale_fill_viridis(discrete=T) +
        scale_color_viridis(discrete=T)
    }else if(plot_by=="time"){
      ggplot(df_time_pts, aes(x=Modifiers,y=Effect)) +
        geom_hline(yintercept = 0)+
        geom_ribbon(aes(ymin=LB, ymax=UB),alpha=0.5, color=F)+
        geom_line()+
        facet_wrap(vars(Week), ncol = 3, labeller = label_both) +
        xlab(ifelse(is.null(mod_name), "Modifier", mod_name)) +
        ylab("Effect") +
        theme_classic() +
        scale_fill_viridis(discrete=T) +
        scale_color_viridis(discrete=T)
    }
  }


}
