#' Plot Cumulative Effects
#' @description Plot estimated cumulative effects from a DLIM object, can also compare estimated cumulative effects between a DLM and DLIM
#' @seealso \link[dlim]{dlim}
#' @seealso Type \code{vignette('dlimOverview')} for a detailed description.
#' @export
#' @import ggplot2 
#' @import reshape2 
#' @import viridis 
#' @import dlnm 
#' @importFrom stats qnorm
#' @importFrom stats predict
#' @importFrom rlang .data
#' @param new_modifiers a vector of new modifier values for prediction (class "\code{numeric}")
#' @param mod_fit DLIM model object (class "\code{dlim}")
#' @param mod_name modifier name that follows variable name nomenclature (class "\code{character}")
#' @param dlm_fit a list containing a \code{crossbasis} object from the \pkg{dlnm} package as the first element and a DLM model object as the second element (class "\code{list}")
#' @param plot_by choose to create plots for particular modifier values, "modifier", or particular time points, "time", (class "\code{character}")
#' @param exposure_time optional vector of exposure-time points if the first time point does not correspond to exposure-time 1. Must have the same length as the number of exposure-time points (class "\code{numeric}")
#' @param exp_time_unit option to provide the unit for the exposure time points, e.g., "month" or "week". Only used with \code{plot_by = "time"} for labeling cross-sections (class "\code{character}")
#' @param time_pts a subset of exposure-time points if \code{plot_by = "time"}. Must be a subset of exposure_time points (class "\code{numeric}")
#' @param mod_trans if modifiers are transformed, specify back transformation function (class "\code{character}")
#' @param link_trans if family for \code{glm} is not Gaussian, specify back transformation to undo link function (class "\code{character}")
#' @return This function returns a ggplot for point-wise effects isolated by either time points or modifier, including a DLM if specified

plot_DLF <- function(new_modifiers, mod_fit, mod_name, 
                     dlm_fit=NULL, plot_by, exposure_time = NULL, 
                     exp_time_unit = "Time", time_pts=NULL,
                     mod_trans = NULL, link_trans = NULL){

  #predict DLIM
  model_pred <- predict(object=mod_fit, newdata = new_modifiers, type="DLF")

  if(!is.null(exposure_time) & length(exposure_time) != ncol(model_pred$est_dlim$betas)){
    stop("The number of exposure-time points provided does not match the number of exposure-time points in the data.")
  }
  
  if(!is.null(dlm_fit)){
    #predict DLM
    cb_dlm <- dlm_fit[[1]]
    model_dlm <- dlm_fit[[2]]
    nlag <- attr(cb_dlm, "lag")[2]+1
    dlm_crosspred <- crosspred(cb_dlm,model_dlm,at=rep(1,nlag),cen = FALSE)
    beta_by_lag <- dlm_crosspred$matfit
    z <- qnorm(1 - (1 - 0.95)/2)
    dlm_lb <- dlm_crosspred$matfit - z * dlm_crosspred$matse #for some reason $matlow is gone
    dlm_ub <- dlm_crosspred$matfit + z * dlm_crosspred$matse #for some reason $mathigh is gone
  }

  #back transform if specified
  if(!is.null(mod_trans)){
    new_modifiers <- do.call(mod_trans,list(new_modifiers))
  }
  
  #exposure-time points default if not specified
  if(is.null(exposure_time)){
    exposure_time <- 1:ncol(model_pred$est_dlim$betas)
  }

  #use all time points if not specified
  if(is.null(time_pts)){
    time_pts <- exposure_time
  }
  
  #make time_pts the index of exposure_time corresponding to subset
  time_pts_idx <- which(exposure_time %in% time_pts)

  #set up dataframe with DLIM estimates
  est_dlf <- model_pred$est_dlim$betas[,time_pts_idx] #modifiers x time_pts
  colnames(est_dlf) <- time_pts
  rownames(est_dlf) <- new_modifiers

  est_lb <- model_pred$est_dlim$LB[,time_pts_idx] #modifiers x time_pts
  colnames(est_lb) <- time_pts
  rownames(est_lb) <- new_modifiers

  est_ub <- model_pred$est_dlim$UB[,time_pts_idx] #modifiers x time_pts
  colnames(est_ub) <- time_pts
  rownames(est_ub) <- new_modifiers

  dlf_df <- melt(est_dlf)
  colnames(dlf_df) <-c("Modifiers", "Time", "Effect")

  lb_df <- melt(est_lb)
  colnames(lb_df) <-c("Modifiers", "Time", "Effect")

  ub_df <- melt(est_ub)
  colnames(ub_df) <-c("Modifiers", "Time", "Effect")

  df <- data.frame(dlf_df, LB=lb_df$Effect, UB=ub_df$Effect)


  #Add the modifiers or Times for plotting by
  if(plot_by=="modifier"){
    df$Modifiers <- signif(df$Modifiers,3)
  }#else if(plot_by=="Time"){
  #   df$Time <- paste("Time", df$Time)
  # }

  #get number of modifier values
  m <- length(new_modifiers)

  #reference line
  ref_line <- 0


  if(!is.null(dlm_fit)){ #DLM and DLIM
    model_type <- attr(mod_fit, "model_type")
    if(model_type=="nonlinear"){
      model_name <- paste0("DLIM(", mod_fit$cb$df_m,",",mod_fit$cb$df_l,")")
    }else{
      model_name <- paste0("DLIM-", model_type, "(", mod_fit$cb$df_l,")")
    }
    df_time_pts <- data.frame(Modifiers = c(df$Modifiers,df$Modifiers),
                           Time = c(df$Time, df$Time),
                           Effect = c(df$Effect, rep(beta_by_lag[time_pts_idx],each=m)),
                           LB = c(df$LB, rep(dlm_lb[time_pts_idx],each=m)),
                           UB = c(df$UB, rep(dlm_ub[time_pts_idx], each=m)),
                           Model = factor(c(rep(model_name, length(df$Modifiers)), rep("DLM", length(df$Modifiers))))
    )

    if(is.null(link_trans)){
      if(mod_fit$fit$family$family!="gaussian"){
        warning("Family is not Gaussian. Use the link_trans argument to transform estimate.")
      }
    }else{
      ref_line <- do.call(link_trans, list(0))
      df_time_pts[,3:5] <- do.call(link_trans,list(df_time_pts[,3:5]))
    }

    if(plot_by=="modifier"){

      colnames(df_time_pts)[which(colnames(df_time_pts)=="Modifiers")] <- mod_name
      ggplot(df_time_pts, aes_string(x="Time",y="Effect", color="Model", fill = "Model")) +
        facet_wrap(mod_name, ncol = 3, labeller = label_both) +
        geom_hline(yintercept = ref_line)+
        geom_ribbon(aes_string(ymin="LB", ymax="UB"), alpha=0.2, color=NA)+
        geom_line()+
        xlab("Time") +
        ylab("Effect") +
        theme_classic() +
        scale_fill_viridis(discrete=TRUE, begin = 0.6, end = 0) +
        scale_color_viridis(discrete=TRUE, begin = 0.6, end = 0)

    }else if(plot_by=="time"){

      ggplot(df_time_pts, aes_string(x="Modifiers",y="Effect", color="Model", fill="Model")) +
        facet_wrap(vars(.data$Time), ncol = 3, labeller = label_both) +
        geom_hline(yintercept = ref_line) +
        geom_ribbon(aes_string(ymin="LB", ymax="UB"),alpha=0.2, color=FALSE)+
        geom_line()+
        xlab(ifelse(is.null(mod_name), "Modifier", mod_name)) +
        ylab("Effect") +
        theme_classic() +
        scale_fill_viridis(discrete=TRUE, begin = 0.6, end = 0) +
        scale_color_viridis(discrete=TRUE, begin = 0.6, end = 0)

    }
  }else{ #just DLIM

    df_time_pts <- data.frame(Modifiers = c(df$Modifiers),
                           Time = c(df$Time),
                           Effect = c(df$Effect),
                           LB = c(df$LB),
                           UB = c(df$UB))

    if(is.null(link_trans)){
      if(mod_fit$fit$family$family!="gaussian"){
        warning("Family is not Gaussian. Use the link_trans argument to transform estimate.")
      }
    }else{
      ref_line <- do.call(link_trans, list(0))
      df_time_pts[,3:5] <- do.call(link_trans,list(df_time_pts[,3:5]))
    }

    if(plot_by=="modifier"){
      colnames(df_time_pts)[which(colnames(df_time_pts)=="Modifiers")] <- mod_name
      ggplot(df_time_pts, aes_string(x="Time",y="Effect")) +
        geom_hline(yintercept = ref_line)+
        geom_ribbon(aes_string(ymin="LB", ymax="UB"), alpha=0.5, color=FALSE)+
        geom_line()+
        facet_wrap(mod_name, ncol = 3, labeller = label_both) +
        xlab("Time") +
        ylab("Effect") +
        theme_classic() +
        scale_fill_viridis(discrete=TRUE, begin = 0.6, end = 0) +
        scale_color_viridis(discrete=TRUE, begin = 0.6, end = 0)
    }else if(plot_by=="time"){
      ggplot(df_time_pts, aes_string(x="Modifiers",y="Effect")) +
        geom_hline(yintercept = 0)+
        geom_ribbon(aes_string(ymin="LB", ymax="UB"),alpha=0.5, color=FALSE)+
        geom_line()+
        facet_wrap(vars(.data$Time), ncol = 3, labeller = label_both) +
        xlab(ifelse(is.null(mod_name), "Modifier", mod_name)) +
        ylab("Effect") +
        theme_classic() +
        scale_fill_viridis(discrete=TRUE, begin = 0.6, end = 0) +
        scale_color_viridis(discrete=TRUE, begin = 0.6, end = 0)
    }
  }


}
