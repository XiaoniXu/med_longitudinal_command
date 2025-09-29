
# Main function to conduct mediation analysis in presence of time-varying mediators, a survival outcome and competing risks in difference in hazards scale

# L: list of models for dropout or censoring
# M: list of models for time-varying mediators
# m: character vector of mediator variable names
# Y: a model object for the survival outcome
# treat: name of the treatment variable
# control.value, treat.value: values of the treatment # peryr: scaling factor, default to 100000 person-years

med_longitudinal=function(L=NULL, M, m, Y, treat='logcocr', control.value=a, treat.value=a_star, data, time_points, peryr=100000){
  
  # Gets the number of rows (subjects) in the dataset
  N=dim(data)[1]
  # How many dropout/censoring models are in L
  NL=length(L)
  # How many mediator models are in M (e.g., 2 if there are M_1, M_2).
  NM=length(M)
  
  # bootstrap loop, 10000 iterations
  boot <- foreach(i=1:1000, .combine='rbind') %dopar% {
    ind <- sample(1:N, replace=TRUE)
   
    # parameters for mediator models
    MModel = list()
    for (i in 1:NM){
      MModel[[i]] <- rmvnorm(1, mean = coef(M[[i]]), sigma = vcov(M[[i]]))
    }
    # parameters for outcome model
    YModel = rmvnorm(1, mean = Y$gamma, sigma = Y$robvar.gamma)
    

    
    # initialize empty matrices to store simulated values
    # store binary predictions of whether each subject drops out at each time point under the counterfactual
    # a: control treatmentf
    # astar: treated
    # predictL_a: dropout under control
    # predictL_astar: dropout under treatment
    # predictL_astar_a: dropout under treatment with a mediator
    # predictM_a: mediator under control
    # predictM_astar: mediator under treatment
    # predictM_a_astar: mediator under control with treatment a*
    PredictL_a <- PredictL_astar <- PredictL_astar_a <- matrix(NA,nrow=N,ncol=NL)
    PredictM_a <- PredictM_astar <- PredictM_a_astar <- matrix(NA,nrow=N,ncol=NM)
    
    # simulate first time-point mediator value under two different counterfactual treatment scenarios for each individual in the bootstrap sample
    # Predict M1
    # create counterfactual datasets for mediator 1
    pred.data.astar.m1 <- pred.data.a.m1 <- model.frame(M[[1]])[ind,]
    # set treatment values
    pred.data.astar.m1[, treat] <- treat.value
    pred.data.a.m1[, treat] <- control.value
    # generate model matrices
    m1mat.astar <- model.matrix(terms(M[[1]]), data = pred.data.astar.m1)
    m1mat.a <- model.matrix(terms(M[[1]]), data = pred.data.a.m1)
    # simulate predicted mediator values
    PredictM_astar[,1] <- tcrossprod(MModel[[1]], m1mat.astar)
    PredictM_a[,1] <- tcrossprod(MModel[[1]], m1mat.a)
    
    # simulate dropout/censoring at the first time point under three counterfactual scenarios
    # Predict L1
    if(NL > 0){
      pred.data.astar.l1 <- pred.data.a.l1 <- pred.data.astar.a.l1 <- model.frame(L[[1]])[ind,]
      pred.data.astar.l1[, treat] <- pred.data.astar.a.l1[, treat] <- treat.value
      pred.data.a.l1[, treat] <- control.value
      pred.data.astar.l1[, m[1]] <- PredictM_astar[,1]
      pred.data.a.l1[, m[1]] <- pred.data.astar.a.l1[, m[1]] <- PredictM_a[,1]
      # simulate dropout indicators via logistic regression
      # PredictL_a: direct effect
      # PredictL_astar: indirect effect through M
      # PredictL_astar_a: indirect effect through dropout
      PredictL_astar_a[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.astar.a.l1, type='response'))
      PredictL_a[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.a.l1, type='response'))
      PredictL_astar[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.astar.l1, type='response'))
    }
    
    # simulate mediator and dropout beyond time 1
    # Predict Li (only if more than one time point)
    if (NM > 1){
      for (i in 2:NM){
        pred.data.a.m <- pred.data.astar.m <- pred.data.a.astar.m <- as.data.frame(matrix(nrow=N, ncol=(dim(model.frame(M[[i]]))[2]-1)))
        colnames(pred.data.a.m) <- colnames(pred.data.astar.m) <- colnames(pred.data.a.astar.m) <- attr(terms(M[[i]]),"term.labels")
        names <- colnames(pred.data.a.m)[which(colnames(pred.data.a.m) %in% attr(terms(M[[1]]),"term.labels"))]
        # set predictor values, fill baseline covariates using bootstrap sample
        pred.data.a.m[, names] <- pred.data.a.astar.m[, names] <- pred.data.astar.m[, names] <- model.frame(M[[1]])[ind,names]
        # assign treatment values for each scenario
        pred.data.a.m[, treat] <- pred.data.a.astar.m[, treat] <- control.value
        pred.data.astar.m[, treat] <- treat.value
        # assign treatment values
        pred.data.a.m[, m[i-1]] <- PredictM_a[,i-1]
        pred.data.astar.m[, m[i-1]] <- PredictM_astar[,i-1]
        pred.data.a.astar.m[, m[i-1]] <- PredictM_a_astar[,i-1]
        if(i==2){
          pred.data.a.astar.m[, m[1]] <- PredictM_a[,1]
        }
        
        # predict M_i
        # case when dropout model exists: only simulate M_i for those not dropped out at tim i-1
        if(NL > 1){
          m1mat.a.m <- model.matrix(~.,data=pred.data.a.m[which(PredictL_a[,i-1]==0),])
          m1mat.astar.m <- model.matrix(~., data = pred.data.astar.m[which(PredictL_astar[,i-1]==0),])
          m1mat.a.astar.m <- model.matrix(~., data = pred.data.a.astar.m[which(PredictL_astar_a[,i-1]==0),])
          
          PredictM_a[which(PredictL_a[,i-1]==0),i] <- tcrossprod(MModel[[i]], m1mat.a.m)
          PredictM_astar[which(PredictL_astar[,i-1]==0),i] <- tcrossprod(MModel[[i]], m1mat.astar.m)
          PredictM_a_astar[which(PredictL_astar_a[,i-1]==0),i] <- tcrossprod(MModel[[i]], m1mat.a.astar.m)
        } else{
          # case when no dropout model
          m1mat.a.m <- model.matrix(~.,data=pred.data.a.m)
          m1mat.astar.m <- model.matrix(~., data = pred.data.astar.m)
          m1mat.a.astar.m <- model.matrix(~., data = pred.data.a.astar.m)
          
          PredictM_a[,i] <- tcrossprod(MModel[[i]], m1mat.a.m)
          PredictM_astar[,i] <- tcrossprod(MModel[[i]], m1mat.astar.m)
          PredictM_a_astar[,i] <- tcrossprod(MModel[[i]], m1mat.a.astar.m)
        }
        # predict dropout at time i
        if(NL > 1 & i<=NL){
          pred.data.a.l <- pred.data.astar.l <- pred.data.astar.a.l <- as.data.frame(matrix(nrow=N, ncol=(dim(model.frame(L[[i]]))[2]-1)))
          colnames(pred.data.a.l) <- colnames(pred.data.astar.l) <- colnames(pred.data.astar.a.l) <- attr(terms(L[[i]]),"term.labels")
          names <- colnames(pred.data.a.l)[which(colnames(pred.data.a.l) %in% attr(terms(L[[1]]),"term.labels"))]
          
          pred.data.a.l[, names] <- pred.data.astar.a.l[, names] <- pred.data.astar.l[, names] <- model.frame(L[[1]])[ind,names]
          pred.data.a.l[, treat] <- control.value
          pred.data.astar.l[, treat] <- pred.data.astar.a.l[, treat] <- treat.value
          
          pred.data.a.l[, m[i]] <- PredictM_a[,i]
          pred.data.astar.l[, m[i]] <- PredictM_astar[,i]
          pred.data.astar.a.l[, m[i]] <- PredictM_a_astar[,i]
          
          PredictL_a[which(PredictL_a[,i-1]==0),i] <- rbinom(length(which(PredictL_a[,i-1]==0)), size=1, prob=predict(L[[i]], pred.data.a.l[which(PredictL_a[,i-1]==0),], type='response'))
          PredictL_astar[which(PredictL_astar[,i-1]==0),i] <- rbinom(length(which(PredictL_astar[,i-1]==0)), size=1, prob=predict(L[[i]], pred.data.astar.l[which(PredictL_astar[,i-1]==0),], type='response'))
          PredictL_astar_a[which(PredictL_astar_a[,i-1]==0),i] <- rbinom(length(which(PredictL_astar_a[,i-1]==0)), size=1, prob=predict(L[[i]], pred.data.astar.a.l[which(PredictL_astar_a[,i-1]==0),], type='response'))
        }
      }
    }
    # by the end, would have PredictM_*: matrix of predicted mediator values for all time points
    # PredictL_*: matrix of predicted dropout indicators for each time point
    
    
    # create counterfactual datasets for predicting the outcome Y (hazard), using the mediator values simulated earlier
    
    # Predict Y
    # Data augmentation method for person-time database
    # PredictY_DEIEM: a*, D1_a, M1_aD1a, D2_aD1aM1aD1a, M2_aD1aM1aD2a
      # direct effect
    # PredictY_TEDE_2: a, D1_a, M1_aD1a, D2_aD1a M1aD1a, M2_aD1aM1aD2a
      # total effect baseline
    # PredictY_IEMIED: a*, D1_a, M1_a*D1a, D2_aD1a M1a*D1a, M2_a*D1aM1a*D2a
      # indirect effect through M
    # PredictY_IEDTE_1: a*, D1_a*, M1_a*D1a*, D2_a*D1a*M1a*D1a*, M2_a*D1a*M1a*D2a*
      # total effect full
    
    pred.data.a.y <- pred.data.astar.y <- pred.data.astar.a.y <- pred.data.astar.a.astar.a.y <- data[ind,c('idno',getvarnames(Y$call)$xvar[-2],m,colnames(data)[grep('time.since.first.exam', colnames(data))])]
    pred.data.a.y[, treat] <- control.value
    pred.data.astar.y[, treat] <- pred.data.astar.a.y[, treat] <- pred.data.astar.a.astar.a.y[, treat] <- treat.value
    pred.data.a.y[, m] <- pred.data.astar.a.y[, m] <- PredictM_a
    pred.data.astar.y[, m] <- PredictM_astar
    pred.data.astar.a.astar.a.y[, m] <- PredictM_a_astar
    pred.data.astar.a.astar.a.y[, m[1]] <- PredictM_a[,1]
    


    
    
    ########################
    # prepare counterfactual datasets for predicting outcomes from the model
    # where mediators and time-varying covariates change across time

    # Build varying lists correctly: one entry per time-varying variable
    M_cols  <- m
    TS_cols <- paste0("time.since.first.exam_", time_points)
    stopifnot(length(M_cols) == length(TS_cols))

    reshape_to_long <- function(df_wide){
      df_wide$id_boot <- seq_len(nrow(df_wide))
      df_long <- reshape(
        df_wide,
        direction = "long",
        varying   = list(M_cols, TS_cols),
        v.names   = c("M","time_since"),
        times     = time_points,
        timevar   = "tidx",
        idvar     = "id_boot"
      )
      df_long[order(df_long$id_boot, df_long$tidx), , drop = FALSE]
    }

    # pred.data.a.y
    df_tv <- reshape_to_long(pred.data.a.y)
    df_pred.data.a.y <- df_tv[, match(getvarnames(Y$call)$xvar, colnames(df_tv)), drop = FALSE]
    df_pred.data.a.y <- model.matrix(~., data = df_pred.data.a.y)[, -1, drop = FALSE]

    # pred.data.astar.y
    df_tv <- reshape_to_long(pred.data.astar.y)
    df_pred.data.astar.y <- df_tv[, match(getvarnames(Y$call)$xvar, colnames(df_tv)), drop = FALSE]
    df_pred.data.astar.y <- model.matrix(~., data = df_pred.data.astar.y)[, -1, drop = FALSE]

    # pred.data.astar.a.astar.a.y
    df_tv <- reshape_to_long(pred.data.astar.a.astar.a.y)
    df_pred.data.astar.a.astar.a.y <- df_tv[, match(getvarnames(Y$call)$xvar, colnames(df_tv)), drop = FALSE]
    df_pred.data.astar.a.astar.a.y <- model.matrix(~., data = df_pred.data.astar.a.astar.a.y)[, -1, drop = FALSE]

    # pred.data.astar.a.y
    df_tv <- reshape_to_long(pred.data.astar.a.y)
    df_pred.data.astar.a.y <- df_tv[, match(getvarnames(Y$call)$xvar, colnames(df_tv)), drop = FALSE]
    df_pred.data.astar.a.y <- model.matrix(~., data = df_pred.data.astar.a.y)[, -1, drop = FALSE]

    
# predict counterfactual outcomes
    # each line computes a predicted hazard under a specific counterfactual scenario
    PredictY_DEIEM <- mean(tcrossprod(YModel, df_pred.data.astar.a.y))
    PredictY_TEDE_2 <- mean(tcrossprod(YModel, df_pred.data.a.y))
    PredictY_IEMIED <- mean(tcrossprod(YModel, df_pred.data.astar.a.astar.a.y))
    PredictY_IEDTE_1 <- mean(tcrossprod(YModel, df_pred.data.astar.y))
    
    # calculate mediation effects
    DE <- mean(PredictY_DEIEM - PredictY_TEDE_2)*100000
    IEM <- mean(PredictY_IEDTE_1 - PredictY_IEMIED)*100000
    IED <- mean(PredictY_IEMIED - PredictY_DEIEM)*100000
    TE <- mean(PredictY_IEDTE_1 - PredictY_TEDE_2)*100000
    # return effects for this sample
    effects <- cbind(DE, IEM, IED, TE)
    effects
  }
  
  # summary step of the function; summarizes the results into point estimates and confidence intervals
  
  # Calculate the effects
  # Direct effect
  DE <- quantile(boot[,1], 0.5)
  DE_low <- quantile(boot[,1], 0.025)
  DE_up <- quantile(boot[,1], 0.975)
  DE_result <- paste0(round(DE,2), ' (', round(DE_low, 2), ', ', round(DE_up, 2), ')')
  
  # Indirect effect through M
  IEM <- quantile(boot[,2], 0.5)
  IEM_low <- quantile(boot[,2], 0.025)
  IEM_up <- quantile(boot[,2], 0.975)
  IEM_result <- paste0(round(IEM,2), ' (', round(IEM_low, 2), ', ', round(IEM_up, 2), ')')
  
  # Indirect effect through D
  IED <- quantile(boot[,3], 0.5)
  IED_low <- quantile(boot[,3], 0.025)
  IED_up <- quantile(boot[,3], 0.975)
  IED_result <- paste0(round(IED,2), ' (', round(IED_low, 2), ', ', round(IED_up, 2), ')')
  
  # Total effect
  TE <- quantile(boot[,4], 0.5)
  TE_low <- quantile(boot[,4], 0.025)
  TE_up <- quantile(boot[,4], 0.975)
  TE_result <- paste0(round(TE,2), ' (', round(TE_low, 2), ', ', round(TE_up, 2), ')')
  
  # Relative indirect effect
  Q <- round(IEM/TE*100,2)
  
  res <- list(DE=DE_result, IEM=IEM_result, IED=IED_result, TE=TE_result, Q=Q)
  return(res)
}

