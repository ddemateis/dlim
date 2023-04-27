#' Simulate Distributed Lag Functions
#' @description generate true distributed lag function values for a given type of simulation
#' @seealso \link[dlim]{sim_data}
#' @export
#' @param L number of lags minus 1
#' @param M vector of modifiers
#' @param type simulation type: 1 is no modification, 2 is linear modification, 3 is complex modification, 4 is stronger complex modification

sim_dlf <- function(L,modifiers,type){
  method<-"normal" #can remove later
  width <- 5
  lags <- 0:L

  curve <- function(l,center,method="quad"){
    if(method=="quad"){
      0.2*pmax(1- ((l-center)^2)/100,0)
    }else if(method=="normal"){
      2.5*dnorm(l,center,width)
    }

  }

  #betas_l <- dnorm(lags,20,width)
  betas_l <- curve(lags,20, method=method) #matching new type 3

  if(type==1){ #no modification
    betas <- matrix(rep(betas_l,length(modifiers)),nrow=L+1,byrow = F) #lag x obs
  }else if(type==2){ #linear modification
    betas <- sweep(matrix(rep(((modifiers+1)/2),L+1),nrow=L+1,byrow = T), MARGIN=1, betas_l, `*`)
  # }else if(type==3){ #complex modification
  #   centers <- M*20 + 10
  #   if(length(M)==1){
  #     betas <- as.matrix(sapply(lags,dnorm,mean=centers,sd=width))
  #   }else{
  #     betas <- t(sapply(lags,dnorm,mean=centers,sd=width))
  #   }
  #   idx <- betas < 0.005
  #   betas[idx] <- 0
  }else if(type==3){ #complex modification second attempt
    centers <- modifiers*20 + 10
    if(length(modifiers)==1){
      betas <- as.matrix(sapply(lags,curve,center=centers, method=method))
    }else{
      betas <- t(sapply(lags,curve,center=centers, method=method))
    }
  }else if(type==4){ #stronger complex modification
    centers <- 37/(1+exp(-20*(modifiers-0.5)))
    if(length(modifiers)==1){
      betas <- as.matrix(sapply(lags,curve,center=centers, method=method))
    }else{
      betas <- t(sapply(lags,curve,center=centers, method=method))
    }
  }else if(type==5){ #type 2 and type 4
    centers <- 37/(1+exp(-20*(modifiers-0.5)))
    if(length(modifiers)==1){
      betas_l <- as.matrix(sapply(lags,curve,center=centers, method=method))
      betas <- matrix(rep(((modifiers+1)/2),L+1),nrow=L+1,byrow = T) * betas_l
    }else{
      betas_l <- t(sapply(lags,curve,center=centers, method=method))
      betas <- matrix(rep(((modifiers+1)/2),L+1),nrow=L+1,byrow = T) * betas_l
    }
  }else{
  # }else if(type==3){ #linear modification function of lag
  #   M_i_l <- sweep(matrix(rep(M,L+1),nrow=L+1,byrow = T), MARGIN=1, c(0:L)/L, `*`) #multiply each row of modifier matrix by lags
  #   betas <- sweep(M_i_l, MARGIN=1, betas_l, `*`) #multiply each col of mod/lag matrix by betas
  # }else if(type==4){ #nonlinear modification
  #   M_i_l <- sweep(matrix(rep(M,L+1),nrow=L+1,byrow = T), MARGIN=1, c(0:L)/L, `*`) #multiply each row of modifier matrix by lags
  #   betas <- sweep(M_i_l, MARGIN=1, betas_l, `^`) #raise each col of mod/lag matrix by betas
  # }else{
    #stop("Incorrect simulation type entered. Choose a type between 1 and 3. See ?dlim::sim_data for more information.")
  }


  return(betas)
}
