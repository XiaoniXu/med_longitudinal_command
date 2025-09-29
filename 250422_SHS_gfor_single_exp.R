
library(SimDesign)
library(timereg)
library(tidyr)
library(dplyr)
library(mice)
library(splitstackshape)
library(foreach)

library(doParallel)
# cl <- makeCluster(detectCores(logical = TRUE))
registerDoParallel(31)

# 8 instead of 31 on my computer


# Load SHS data
# setwd("/home/ardore/ftp_ardore/Epigenetica/MESA")
# setwd("/Volumes/emergentese/datos/jobscrew/ftp/ardore/Epigenetica/dat/CVD")
setwd("/Users/ad3531/Library/CloudStorage/OneDrive-ColumbiaUniversityIrvingMedicalCenter/Epigenetica/dat/CVD")
shs_cvd <- read.csv("shs_ccvd2017.csv")
names(shs_cvd) <- casefold(names(shs_cvd))
shs_cvd <- shs_cvd[,c(1,5:63, 65:72, 74:164)]
# setwd("/Volumes/emergentese/datos/jobscrew/ftp/ardore/Epigenetica/dat/")
setwd("/Users/ad3531/Library/CloudStorage/OneDrive-ColumbiaUniversityIrvingMedicalCenter/Epigenetica/dat")
shs.metal <- read.csv("SHS_P123_100421.csv")
shs.metal <- shs.metal[!duplicated(shs.metal$idno),]
shs.metal <- shs.metal[,c(1:661, 664:783)]
shs.metal <- merge(shs.metal, shs_cvd, by='idno', all=FALSE)


# Log-transformed metal variables
# cd.cr.ln, zn.cr.ln, u.cr.ln, w.cr.ln, sb.cr.ln, se.cr.ln


# Remove missing values
shs.metal$s1ldl[!is.na(shs.metal$s1ldlest)] <- shs.metal$s1ldlest[!is.na(shs.metal$s1ldlest)]
shs.metal$s1ldl[is.na(shs.metal$s1ldlest)]  <- shs.metal$s1ldlbq[is.na(shs.metal$s1ldlest)]
shs.metal <- shs.metal[!is.na(shs.metal$s1bmi) & !is.na(shs.metal$s1ldl)&                                       #    - LDL cholesterol ** derived from mix of LDLEST & LDLBQ
                         !is.na(shs.metal$s1sbp) & !is.na(shs.metal$s1dmwhb) & !is.na(shs.metal$s1acr) &
                         !is.na(shs.metal$s1hdl) & !is.na(shs.metal$cd.cr.ln) &
                         !is.na(shs.metal$sum.ias.cr.ln) & !is.na(shs.metal$s1ckdepi.new),]
shs.metal <- shs.metal[-which(shs.metal$s1smoke==''),]

# Convert variables to date format
shs.metal$closingdate <- as.Date(as.character(shs.metal$closingdate), format="%m/%d/%Y")
shs.metal$s1examdt <- as.Date(as.character(shs.metal$s1examdt), format="%m/%d/%Y")
shs.metal$s2examdt <- as.Date(as.character(shs.metal$s2examdt), format="%m/%d/%Y")
shs.metal$s3examdt <- as.Date(as.character(shs.metal$s3examdt), format="%m/%d/%Y")
shs.metal$dod <- as.Date(as.character(shs.metal$dod), format="%m/%d/%Y")
shs.metal$anycvddts1 <- as.Date(as.character(shs.metal$anycvddts1), format="%m/%d/%Y")

# Restrict follow-up to 2009
shs.metal$year <- as.numeric(substr(shs.metal$anycvddts1,1,4))
shs.metal$anycvds1[which(!is.na(shs.metal$anycvddts1) & shs.metal$anycvddts1 > as.Date("12/31/2009", format="%m/%d/%Y"))] <-  0
shs.metal$anycvddts1[which(!is.na(shs.metal$anycvddts1) & shs.metal$anycvddts1 > as.Date("12/31/2009", format="%m/%d/%Y"))] <-  as.Date("12/31/2009", format="%m/%d/%Y")
shs.metal[which(!is.na(shs.metal$dod) & shs.metal$dod>as.Date("12/31/2009", format="%m/%d/%Y")), 'dod'] <- NA
shs.metal$closingdate <- rep(as.Date("12/31/2009", format="%m/%d/%Y"), dim(shs.metal)[1])

# CVD follow-up
shs.metal$peryr.exam.anycvd <- NULL
shs.metal$peryr.exam.anycvd <- ifelse(shs.metal$closingdate > shs.metal$anycvddts1,(shs.metal$anycvddts1 - shs.metal$s1examdt)/365.25,
                                      (shs.metal$closingdate - shs.metal$s1examdt)/365.25)
shs.metal$peryr.exam.anycvd[is.na(shs.metal$anycvddts1)&!is.na(shs.metal$dod)] <- (shs.metal$dod[is.na(shs.metal$anycvddts1)&!is.na(shs.metal$dod)] -
                                                                                     shs.metal$s1examdt[is.na(shs.metal$anycvddts1)&!is.na(shs.metal$dod)])/365.25
shs.metal$peryr.exam.anycvd[is.na(shs.metal$anycvddts1)&is.na(shs.metal$peryr.exam.anycvd)] <-
  (shs.metal$closingdate[is.na(shs.metal$anycvddts1)&is.na(shs.metal$peryr.exam.anycvd)] - shs.metal$s1examdt[is.na(shs.metal$anycvddts1)&is.na(shs.metal$peryr.exam.anycvd)])/365.25
shs.metal$cvda <-  ifelse(shs.metal$anycvds1==1&!is.na(shs.metal$anycvds1),1,0) 

# Mortality follow-up
shs.metal$peryr.exam.anycvd.death <- NULL
shs.metal$peryr.exam.anycvd.death <- ifelse((shs.metal$closingdate > shs.metal$dod), (shs.metal$dod - shs.metal$s1examdt)/365.25,
                                            (shs.metal$closingdate - shs.metal$s1examdt)/365.25)
shs.metal$peryr.exam.anycvd.death[is.na(shs.metal$dod)] <- (shs.metal$closingdate[is.na(shs.metal$dod)] - shs.metal$s1examdt[is.na(shs.metal$dod)])/365.25
shs.metal$death <- ifelse(!is.na(shs.metal$dod),1,0)


# Creation of a death indicator for each visit
shs.metal$death_2 <- ifelse(is.na(shs.metal$s2examdt) & (!is.na(shs.metal$dod) & (shs.metal$dod < "1995-12-18")) & shs.metal$cvda!=1,1,0)
shs.metal$death_3 <- ifelse((is.na(shs.metal$s3examdt) & (!is.na(shs.metal$dod) & (shs.metal$dod < "1999-12-18")) & shs.metal$cvda!=1) | shs.metal$death_2==1,1,0)
shs.metal$death_4 <- ifelse((!is.na(shs.metal$dod) & (shs.metal$dod < "2009-12-18") & shs.metal$cvda!=1) | shs.metal$death_3==1,1,0)

# Creation of the variable "time since first exam"
shs.metal$time.since.first.exam_1 <- rep(0, dim(shs.metal)[1])
shs.metal$time.since.first.exam_2 <- ifelse(!is.na(shs.metal$s2examdt), (shs.metal$s2examdt - shs.metal$s1examdt)/365.25, NA)
shs.metal$time.since.first.exam_3 <- ifelse(!is.na(shs.metal$s3examdt), (shs.metal$s3examdt - shs.metal$s1examdt)/365.25, NA)


# Selection of relevant variables from the database
shs.metal <- shs.metal[,c('idno', 's1sbp', 's2sbp', 's3sbp',
                          's1g0', 's2g0', 's3g0',
                          's1age', 'sex', 'cd.cr.ln', 'zn.cr.ln', 'u.cr.ln', 'sum.ias.cr.ln',
                          'peryr.exam.anycvd', 'cvda', 'death_2', 'death_3', 'death_4',
                          'center', 's1smoke', 's1bmi', 's1ldl', 's1dmwhb', 's1acr', 's1hdl', 's1ckdepi.new', 's1htnrx2',
                          'time.since.first.exam_1', 'time.since.first.exam_2', 'time.since.first.exam_3')]


# Imputation of missing values in sbp using MICE
init = mice(shs.metal, maxit=0) 
meth = init$method
meth[5:30] <- ""
predM = init$predictorMatrix
predM['s2sbp', ] <- 0
predM['s2sbp', c('s1sbp', 's1age', 'sex', 'center', 's1bmi', 's1ckdepi.new',
                 's1smoke', 's1ldl', 's1hdl', 's1dmwhb', 's1htnrx2', 's1acr')] <- 1
predM['s3sbp', ] <- 0
predM['s3sbp', c('s1sbp', 's2sbp', 's1age', 'sex', 'center', 's1bmi', 's1ckdepi.new',
                 's1smoke', 's1ldl', 's1hdl', 's1dmwhb', 's1htnrx2', 's1acr')] <- 1

set.seed(111)
imputed_1 = mice(shs.metal, method=meth, predictorMatrix=predM, m=5)
# Change the action parameter to get each of the imputed datasets
imputed <- complete(imputed_1, action = 4)


# Comprobation
# shs.metal[1:15,c('s1sbp',  's2sbp',  's3sbp')]
# imputed[1:15,c('s1sbp',  's2sbp',  's3sbp')]


# Assign NA to sbp in visits in which the participant had already died
for (i in 1:dim(imputed)[1]){
  if(imputed[i,'death_2']==1){
    imputed[i,'s2sbp'] <- imputed[i,'s3sbp'] <- NA
  } else if(imputed[i,'death_3']==1){
    imputed[i,'s3sbp'] <- NA
  }}


# Convert categorical variables to factors
imputed$death_2 <- as.factor(imputed$death_2)
imputed$death_3 <- as.factor(imputed$death_3)
imputed$death_4 <- as.factor(imputed$death_4)
imputed$sex <- as.factor(imputed$sex)
imputed$center <- as.factor(imputed$center)
imputed$s1smoke <- as.factor(imputed$s1smoke)
imputed$s1acr <- as.factor(imputed$s1acr)
imputed$s1dmwhb <- as.factor(imputed$s1dmwhb)
imputed$s1acr <- as.factor(imputed$s1acr)


# Impute time since first exam to those that have it missing
for (i in 1:(dim(imputed)[1])){
  if(is.na(imputed$time.since.first.exam_2[i])){
    imputed$time.since.first.exam_2[i] <- mean(imputed$time.since.first.exam_2, na.rm=T)
  }
  if(is.na(imputed$time.since.first.exam_3[i])){
    imputed$time.since.first.exam_3[i] <- mean(imputed$time.since.first.exam_3, na.rm=T)
}}

colnames(imputed)[which(colnames(imputed)=='s1sbp')] <- 'sbp_1'
colnames(imputed)[which(colnames(imputed)=='s2sbp')] <- 'sbp_2'
colnames(imputed)[which(colnames(imputed)=='s3sbp')] <- 'sbp_3'



# Bootstrap
start_time <- Sys.time()

boot <- foreach(i=1:1000, .combine='rbind') %dopar% {
  ind <- sample(1:dim(imputed)[1], replace=TRUE)
  data <- imputed[ind,]
  data$id_unique <- 1:dim(data)[1]
  
  M1 = lm(sbp_1 ~ cd.cr.ln + s1age + sex + center + s1smoke + s1bmi + s1ckdepi.new, data = data)
  L1 = glm(death_2 ~ cd.cr.ln + sbp_1 + s1age + sex + center + s1smoke + s1bmi + s1ckdepi.new, data = data, family=binomial(link='logit'))
  M2 = lm(sbp_2 ~ cd.cr.ln + sbp_1 + s1age + sex + center + s1smoke + s1bmi + s1ckdepi.new, data = data)
  L2 = glm(death_3 ~ cd.cr.ln + sbp_2 + s1age + sex + center + s1smoke + s1bmi + s1ckdepi.new, data = data, family=binomial(link='logit'))
  # no M3/L3
  M3 = lm(sbp_3 ~ cd.cr.ln + sbp_2 + s1age + sex + center + s1smoke + s1bmi + s1ckdepi.new, data = data)
  L3 = glm(death_4 ~ cd.cr.ln + sbp_3 + s1age + sex + center + s1smoke + s1bmi + s1ckdepi.new, data = data, family=binomial(link='logit'))
  
  # Data augmentation method for counting process format for the outcome model
  df <- tmerge(data, data, id=id_unique, endpt=event(peryr.exam.anycvd,cvda))
  # df[250:300,c('idno', 'peryr.exam.anycvd', 'tstart', 'tstop', 'cvda', 'endpt')]
  
  # Convert the already data dataset with the death indicator from wide to long format
  df_tv <- reshape(data, direction = "long", varying = c("sbp_1", "time.since.first.exam_1", "sbp_2", "time.since.first.exam_2", "sbp_3",     
                                                            "time.since.first.exam_3"), sep = "_", times=c('1','2','3'), idvar='id_unique')
  df_tv <- df_tv[order(df_tv$id_unique),]
  df <- tmerge(df, df_tv, id=id_unique, sbp=tdc(time.since.first.exam,sbp))
  
  # In the Aalen model, factor variables need to be introduced as dummy variables for the getvarnames() fuction to work
  Y = aalen(Surv(tstart, tstop, endpt) ~ const(cd.cr.ln) + const(sbp) +
              const(s1age) + const(sex) + const(center) + const(s1smoke) + const(s1bmi) +
              const(s1ckdepi.new) , data = df, resample.iid=1)
  # + const(s1ldl) + const(s1hdl) + const(s1dmwhb) + const(s1acr)
  
  # Mediation function parameters specification
  treat <- 'cd.cr.ln'
  L = list(L1=L1, L2=L2, L3=L3)
  M = list(M1=M1, M2=M2, M3=M3)
  m <- c('sbp_1', 'sbp_2', 'sbp_3')
  time_points <- c("time.since.first.exam_1", "time.since.first.exam_2", "time.since.first.exam_3") # adapt to yuchen's data (take one out)
  a <- quantile(data$cd.cr.ln, probs=0.25)
  a_star <- quantile(data$cd.cr.ln, probs=0.75)
  control.value=a
  treat.value=a_star
  peryr=100000
  
  # Mediational g formula call
  mesa_mediation <- med_longitudinal(L=L, M, m, Y, treat=treat, control.value=a, treat.value=a_star, data, time_points, peryr=100000)
  # mesa_mediation <- med_sp(L=L, M, m, Y, treat=treat, control.value=a, treat.value=a_star, data, time_points, t=20)
  mesa_mediation
}

end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

# Calculate the effects
# Direct effect
DE <- quantile(boot[,1], 0.5)
DE_low <- quantile(boot[,1], 0.025)
DE_up <- quantile(boot[,1], 0.975)
DE_result <- paste0(round(DE,3), ' (', round(DE_low, 3), ', ', round(DE_up, 3), ')')

# Indirect effect through M
IEM <- quantile(boot[,2], 0.5)
IEM_low <- quantile(boot[,2], 0.025)
IEM_up <- quantile(boot[,2], 0.975)
IEM_result <- paste0(round(IEM,3), ' (', round(IEM_low, 3), ', ', round(IEM_up, 3), ')')

# Indirect effect through D
IED <- quantile(boot[,3], 0.5)
IED_low <- quantile(boot[,3], 0.025)
IED_up <- quantile(boot[,3], 0.975)
IED_result <- paste0(round(IED,3), ' (', round(IED_low, 3), ', ', round(IED_up, 3), ')')

# Total effect
TE <- quantile(boot[,4], 0.5)
TE_up <- quantile(boot[,4], 0.975)
TE_low <- quantile(boot[,4], 0.025)
TE_result <- paste0(round(TE,3), ' (', round(TE_low, 3), ', ', round(TE_up, 3), ')')

res <- list(DE=DE_result, IEM=IEM_result, IED=IED_result, TE=TE_result)
res

# Relative indirect effect
# Q <- round(IEM/TE*100,2)

# Calculate the effects
# Direct effect
DE <- quantile(boot[,1], 0.5) # median, point estimate
DE_low <- quantile(boot[,1], 0.0041) # 2.5% percentile, 0.025
DE_up <- quantile(boot[,1], 0.9958) # 0.975
DE_result <- paste0(round(DE,3), ' (', round(DE_low, 3), ', ', round(DE_up, 3), ')')

# Indirect effect through M
IEM <- quantile(boot[,2], 0.5)
IEM_low <- quantile(boot[,2], 0.0041)
IEM_up <- quantile(boot[,2], 0.9958)
IEM_result <- paste0(round(IEM,3), ' (', round(IEM_low, 3), ', ', round(IEM_up, 3), ')')

# Indirect effect through D
IED <- quantile(boot[,3], 0.5)
IED_low <- quantile(boot[,3], 0.0041)
IED_up <- quantile(boot[,3], 0.9958)
IED_result <- paste0(round(IED,3), ' (', round(IED_low, 3), ', ', round(IED_up, 3), ')')

# Total effect
TE <- quantile(boot[,4], 0.5)
TE_low <- quantile(boot[,4], 0.0041)
TE_up <- quantile(boot[,4], 0.9958)
TE_result <- paste0(round(TE,3), ' (', round(TE_low, 3), ', ', round(TE_up, 3), ')')

res <- list(DE=DE_result, IEM=IEM_result, IED=IED_result, TE=TE_result)
res



# Functions

# Get variable names of a regression model (suitable for Aalen additive hazards models)
getvarnames <- function(formula, data = NULL)
{
  if (is.character(formula))
    return(list(varnames=formula, xvar=formula, yvar=NULL))
  if (is.null(formula)) return(list(varnames=NULL, xvar=NULL, yvar=NULL))
  
  formula <- formula(formula)
  lyv <- NULL
  lxv <- lvnm <- all.vars(formula[1:2])
  if (length(formula)==3) {
    lyv <- lxv 
    lxv <- all.vars(formula[-2])
    if ("." %in% lxv) {
      if (length(data)==0)
        stop("!getvarnames! '.' in formula and no 'data'")
      lform <- formula(terms(formula, data=data))
      lxv <- all.vars(lform[-2])
    }
    lvnm <- c(lxv, lvnm)
  }
  list(varnames=lvnm, xvar=lxv, yvar=lyv)
}




# Mediation main function
med_longitudinal=function(L=NULL, M, m, Y, treat='logcocr', control.value=a, treat.value=a_star, data, time_points, peryr=100000){
  
  N=dim(data)[1]
  NL=length(L)
  NM=length(M)
  
  MModel = list()
  for (i in 1:NM){
    MModel[[i]] <- rmvnorm(1, mean = coef(M[[i]]), sigma = vcov(M[[i]]))
  }
  
  YModel = rmvnorm(1, mean = Y$gamma, sigma = Y$robvar.gamma)
  
  PredictL_a <- PredictL_astar <- PredictL_astar_a <- matrix(NA,nrow=N,ncol=NL)
  PredictM_a <- PredictM_astar <- PredictM_a_astar <- matrix(NA,nrow=N,ncol=NM)
  
  # Predict M1
  pred.data.astar.m1 <- pred.data.a.m1 <- model.frame(M[[1]])
  pred.data.astar.m1[, treat] <- treat.value
  pred.data.a.m1[, treat] <- control.value
  
  m1mat.astar <- model.matrix(terms(M[[1]]), data = pred.data.astar.m1)
  m1mat.a <- model.matrix(terms(M[[1]]), data = pred.data.a.m1)
  
  PredictM_astar[,1] <- tcrossprod(MModel[[1]], m1mat.astar)
  PredictM_a[,1] <- tcrossprod(MModel[[1]], m1mat.a)
  
  # Predict L1
  if(NL > 0){
    pred.data.astar.l1 <- pred.data.a.l1 <- pred.data.astar.a.l1 <- model.frame(L[[1]])
    pred.data.astar.l1[, treat] <- pred.data.astar.a.l1[, treat] <- treat.value
    pred.data.a.l1[, treat] <- control.value
    pred.data.astar.l1[, m[1]] <- PredictM_astar[,1]
    pred.data.a.l1[, m[1]] <- pred.data.astar.a.l1[, m[1]] <- PredictM_a[,1]
    
    PredictL_astar_a[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.astar.a.l1, type='response'))
    PredictL_a[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.a.l1, type='response'))
    PredictL_astar[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.astar.l1, type='response'))
  }
  
  # Predict Li (only if more than one time point)
  if (NM > 1){
    for (i in 2:NM){
      pred.data.a.m <- pred.data.astar.m <- pred.data.a.astar.m <- as.data.frame(matrix(nrow=N, ncol=(dim(model.frame(M[[i]]))[2]-1)))
      colnames(pred.data.a.m) <- colnames(pred.data.astar.m) <- colnames(pred.data.a.astar.m) <- attr(terms(M[[i]]),"term.labels")
      names <- colnames(pred.data.a.m)[which(colnames(pred.data.a.m) %in% attr(terms(M[[1]]),"term.labels"))]
      
      pred.data.a.m[, names] <- pred.data.a.astar.m[, names] <- pred.data.astar.m[, names] <- model.frame(M[[1]])[,names]
      pred.data.a.m[, treat] <- pred.data.a.astar.m[, treat] <- control.value
      pred.data.astar.m[, treat] <- treat.value
      
      pred.data.a.m[, m[i-1]] <- PredictM_a[,i-1]
      pred.data.astar.m[, m[i-1]] <- PredictM_astar[,i-1]
      pred.data.a.astar.m[, m[i-1]] <- PredictM_a_astar[,i-1]
      if(i==2){
        pred.data.a.astar.m[, m[1]] <- PredictM_a[,1]
      }
      
      if(NL > 1){
        m1mat.a.m <- model.matrix(~.,data=pred.data.a.m[which(PredictL_a[,i-1]==0),])
        m1mat.astar.m <- model.matrix(~., data = pred.data.astar.m[which(PredictL_astar[,i-1]==0),])
        m1mat.a.astar.m <- model.matrix(~., data = pred.data.a.astar.m[which(PredictL_astar_a[,i-1]==0),])
        
        PredictM_a[which(PredictL_a[,i-1]==0),i] <- tcrossprod(MModel[[i]], m1mat.a.m)
        PredictM_astar[which(PredictL_astar[,i-1]==0),i] <- tcrossprod(MModel[[i]], m1mat.astar.m)
        PredictM_a_astar[which(PredictL_astar_a[,i-1]==0),i] <- tcrossprod(MModel[[i]], m1mat.a.astar.m)
      } else{
        m1mat.a.m <- model.matrix(~.,data=pred.data.a.m)
        m1mat.astar.m <- model.matrix(~., data = pred.data.astar.m)
        m1mat.a.astar.m <- model.matrix(~., data = pred.data.a.astar.m)
        
        PredictM_a[,i] <- tcrossprod(MModel[[i]], m1mat.a.m)
        PredictM_astar[,i] <- tcrossprod(MModel[[i]], m1mat.astar.m)
        PredictM_a_astar[,i] <- tcrossprod(MModel[[i]], m1mat.a.astar.m)
      }
      
      if(NL > 1 & i<=NL){
        pred.data.a.l <- pred.data.astar.l <- pred.data.astar.a.l <- as.data.frame(matrix(nrow=N, ncol=(dim(model.frame(L[[i]]))[2]-1)))
        colnames(pred.data.a.l) <- colnames(pred.data.astar.l) <- colnames(pred.data.astar.a.l) <- attr(terms(L[[i]]),"term.labels")
        names <- colnames(pred.data.a.l)[which(colnames(pred.data.a.l) %in% attr(terms(L[[1]]),"term.labels"))]
        
        pred.data.a.l[, names] <- pred.data.astar.a.l[, names] <- pred.data.astar.l[, names] <- model.frame(L[[1]])[,names]
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
  
  
  # Predict Y
  # Data augmentation method for person-time database
  # PredictY_DEIEM: a*, D1_a, M1_aD1a, D2_aD1aM1aD1a, M2_aD1aM1aD2a
  # PredictY_TEDE_2: a, D1_a, M1_aD1a, D2_aD1a M1aD1a, M2_aD1aM1aD2a
  # PredictY_IEMIED: a*, D1_a, M1_a*D1a, D2_aD1a M1a*D1a, M2_a*D1aM1a*D2a
  # PredictY_IEDTE_1: a*, D1_a*, M1_a*D1a*, D2_a*D1a*M1a*D1a*, M2_a*D1a*M1a*D2a*
  
  pred.data.a.y <- pred.data.astar.y <- pred.data.astar.a.y <- pred.data.astar.a.astar.a.y <- data[,c('id_unique',getvarnames(Y$call)$xvar[-2],m,colnames(data)[grep('time.since.first.exam', colnames(data))])]
  pred.data.a.y[, treat] <- control.value
  pred.data.astar.y[, treat] <- pred.data.astar.a.y[, treat] <- pred.data.astar.a.astar.a.y[, treat] <- treat.value
  pred.data.a.y[, m] <- pred.data.astar.a.y[, m] <- PredictM_a
  pred.data.astar.y[, m] <- PredictM_astar
  pred.data.astar.a.astar.a.y[, m] <- PredictM_a_astar
  pred.data.astar.a.astar.a.y[, m[1]] <- PredictM_a[,1]
  
  
  ########################
  
  # Data augmentation method for the counterfactuals
  vector_time_points <- c()
  for (i in 1:length(time_points)){
    vector_time_points <- c(vector_time_points, m[i], time_points[i])
  }
  
  # pred.data.a.y
  df_tv <- reshape(pred.data.a.y, direction = "long", varying = vector_time_points,
                   sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_unique')
  df_tv <- df_tv[order(df_tv$id_unique),]
  df_pred.data.a.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
  df_pred.data.a.y <- model.matrix(~.,data=df_pred.data.a.y)[,-1]
  
  # pred.data.astar.y
  df_tv <- reshape(pred.data.astar.y, direction = "long", varying = vector_time_points,
                   sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_unique')
  df_tv <- df_tv[order(df_tv$id_unique),]
  df_pred.data.astar.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
  df_pred.data.astar.y <- model.matrix(~.,data=df_pred.data.astar.y)[,-1]
  
  # pred.data.astar.a.astar.a.y
  df_tv <- reshape(pred.data.astar.a.astar.a.y, direction = "long", varying = vector_time_points,
                   sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_unique')
  df_tv <- df_tv[order(df_tv$id_unique),]
  df_pred.data.astar.a.astar.a.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
  df_pred.data.astar.a.astar.a.y <- model.matrix(~.,data=df_pred.data.astar.a.astar.a.y)[,-1]
  
  # pred.data.astar.a.y
  df_tv <- reshape(pred.data.astar.a.y, direction = "long", varying = vector_time_points,
                   sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_unique')
  df_tv <- df_tv[order(df_tv$id_unique),]
  df_pred.data.astar.a.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
  df_pred.data.astar.a.y <- model.matrix(~.,data=df_pred.data.astar.a.y)[,-1]
  
  #######################
  
  PredictY_DEIEM <- mean(tcrossprod(YModel, df_pred.data.astar.a.y))
  PredictY_TEDE_2 <- mean(tcrossprod(YModel, df_pred.data.a.y))
  PredictY_IEMIED <- mean(tcrossprod(YModel, df_pred.data.astar.a.astar.a.y))
  PredictY_IEDTE_1 <- mean(tcrossprod(YModel, df_pred.data.astar.y))
  
  DE <- mean(PredictY_DEIEM - PredictY_TEDE_2)*100000
  IEM <- mean(PredictY_IEDTE_1 - PredictY_IEMIED)*100000
  IED <- mean(PredictY_IEMIED - PredictY_DEIEM)*100000
  TE <- mean(PredictY_IEDTE_1 - PredictY_TEDE_2)*100000
  effects <- cbind(DE, IEM, IED, TE)
  return(effects)
}




med_sp=function(L=NULL, M, m, Y, treat='logcocr', control.value=a, treat.value=a_star, data, time_points, t=20){
  
  N=dim(data)[1]
  NL=length(L)
  NM=length(M)
  
  PredictL_a <- PredictL_astar <- PredictL_astar_a <- matrix(NA,nrow=N,ncol=NL)
  PredictM_a <- PredictM_astar <- PredictM_a_astar <- matrix(NA,nrow=N,ncol=NM)
  
  # Predict M1
  pred.data.astar.m1 <- pred.data.a.m1 <- model.frame(M[[1]])
  pred.data.astar.m1[, treat] <- treat.value
  pred.data.a.m1[, treat] <- control.value
  PredictM_astar[,1] <- predict(M[[1]], pred.data.astar.m1)
  PredictM_a[,1] <- predict(M[[1]], pred.data.a.m1)
  
  # Predict L1  
  if(NL > 1){
    pred.data.astar.l1 <- pred.data.a.l1 <- pred.data.astar.a.l1 <- model.frame(L[[1]])
    pred.data.astar.l1[, treat] <- pred.data.astar.a.l1[, treat] <- treat.value
    pred.data.a.l1[, treat] <- control.value
    pred.data.astar.l1[, m[1]] <- PredictM_astar[,1]
    pred.data.a.l1[, m[1]] <- pred.data.astar.a.l1[, m[1]] <- PredictM_a[,1]
    
    PredictL_astar_a[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.astar.a.l1, type='response'))
    PredictL_a[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.a.l1, type='response'))
    PredictL_astar[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.astar.l1, type='response'))
  }
  
  # Predict Li and Mi (only if more than one time point)
  if (NM > 1){
    for (i in 2:NM){
      pred.data.a.m <- pred.data.astar.m <- pred.data.a.astar.m <- as.data.frame(matrix(nrow=N, ncol=(dim(model.frame(M[[i]]))[2]-1)))
      colnames(pred.data.a.m) <- colnames(pred.data.astar.m) <- colnames(pred.data.a.astar.m) <- attr(terms(M[[i]]),"term.labels")
      names <- colnames(pred.data.a.m)[which(colnames(pred.data.a.m) %in% attr(terms(M[[1]]),"term.labels"))]
      
      pred.data.a.m[, names] <- pred.data.a.astar.m[, names] <- pred.data.astar.m[, names] <- model.frame(M[[1]])[,names]
      pred.data.a.m[, treat] <- pred.data.a.astar.m[, treat] <- control.value
      pred.data.astar.m[, treat] <- treat.value
      
      pred.data.a.m[, m[i-1]] <- PredictM_a[,i-1]
      pred.data.astar.m[, m[i-1]] <- PredictM_astar[,i-1]
      pred.data.a.astar.m[, m[i-1]] <- PredictM_a_astar[,i-1]
      if(i==2){
        pred.data.a.astar.m[, m[1]] <- PredictM_a[,1]
      }
      if(NL > 1){
        PredictM_a[which(PredictL_a[,i-1]==0),i] <- predict(M[[i]], pred.data.a.m[which(PredictL_a[,i-1]==0),])
        PredictM_astar[which(PredictL_astar[,i-1]==0),i] <- predict(M[[i]], pred.data.astar.m[which(PredictL_astar[,i-1]==0),])
        PredictM_a_astar[which(PredictL_astar_a[,i-1]==0),i] <- predict(M[[i]], pred.data.a.astar.m[which(PredictL_astar_a[,i-1]==0),])
      } else{
        PredictM_a[,i] <- predict(M[[i]], pred.data.a.m)
        PredictM_astar[,i] <- predict(M[[i]], pred.data.astar.m)
        PredictM_a_astar[,i] <- predict(M[[i]], pred.data.a.astar.m)
      }
      
      if(NL > 1 & i<=NL){
        pred.data.a.l <- pred.data.astar.l <- pred.data.astar.a.l <- as.data.frame(matrix(nrow=N, ncol=(dim(model.frame(L[[i]]))[2]-1)))
        colnames(pred.data.a.l) <- colnames(pred.data.astar.l) <- colnames(pred.data.astar.a.l) <- attr(terms(L[[i]]),"term.labels")
        names <- colnames(pred.data.a.l)[which(colnames(pred.data.a.l) %in% attr(terms(L[[1]]),"term.labels"))]
        
        pred.data.a.l[, names] <- pred.data.astar.a.l[, names] <- pred.data.astar.l[, names] <- model.frame(L[[1]])[,names]
        pred.data.a.l[, treat] <- control.value
        pred.data.astar.l[, treat] <- pred.data.astar.a.l[, treat] <- treat.value
        
        pred.data.a.l[, m[i]] <- PredictM_a[,i]
        pred.data.astar.l[, m[i]] <- PredictM_astar[,i]
        pred.data.astar.a.l[, m[i]] <- PredictM_a_astar[,i]
        
        PredictL_a[which(PredictL_a[,i-1]==0),i] <- rbinom(length(which(PredictL_a[,i-1]==0)), size=1, prob=predict(L[[i]], pred.data.a.l[which(PredictL_a[,i-1]==0),], type='response'))
        PredictL_astar[which(PredictL_astar[,i-1]==0),i] <- rbinom(length(which(PredictL_astar[,i-1]==0)), size=1, prob=predict(L[[i]], pred.data.astar.l[which(PredictL_astar[,i-1]==0),], type='response'))
        PredictL_astar_a[which(PredictL_astar_a[,i-1]==0),i] <- rbinom(length(which(PredictL_astar_a[,i-1]==0)), size=1, prob=predict(L[[i]], pred.data.astar.a.l[which(PredictL_astar_a[,i-1]==0),], type='response'))
      }}}
  
  # Predict Y
  # Data augmentation method for person-time database
  # PredictY_DEIEM: a*, D1_a, M1_aD1a, D2_aD1aM1aD1a, M2_aD1aM1aD2a
  # PredictY_TEDE_2: a, D1_a, M1_aD1a, D2_aD1a M1aD1a, M2_aD1aM1aD2a
  # PredictY_IEMIED: a*, D1_a, M1_a*D1a, D2_aD1a M1a*D1a, M2_a*D1aM1a*D2a
  # PredictY_IEDTE_1: a*, D1_a*, M1_a*D1a*, D2_a*D1a*M1a*D1a*, M2_a*D1a*M1a*D2a*
  
  pred.data.a.y <- pred.data.astar.y <- pred.data.astar.a.y <- pred.data.astar.a.astar.a.y <- data[,c('id_unique',getvarnames(Y$call)$xvar[-2],m,colnames(data)[grep('time.since.first.exam', colnames(data))])]
  pred.data.a.y[, treat] <- control.value
  pred.data.astar.y[, treat] <- pred.data.astar.a.y[, treat] <- pred.data.astar.a.astar.a.y[, treat] <- treat.value
  pred.data.a.y[, m] <- pred.data.astar.a.y[, m] <- PredictM_a
  pred.data.astar.y[, m] <- PredictM_astar
  pred.data.astar.a.astar.a.y[, m] <- PredictM_a_astar
  pred.data.astar.a.astar.a.y[, m[1]] <- PredictM_a[,1]
  
  
  ########################
  
  # Data augmentation method for the counterfactuals
  vector_time_points <- c()
  for (i in 1:length(time_points)){
    vector_time_points <- c(vector_time_points, m[i], time_points[i])
  }
  
  # pred.data.a.y
  df_tv <- reshape(pred.data.a.y, direction = "long", varying = vector_time_points,
                   sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_unique')
  df_tv <- df_tv[order(df_tv$id_unique),]
  df_pred.data.a.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
  df_pred.data.a.y <- model.matrix(~.,data=df_pred.data.a.y)[,-1]
  
  # pred.data.astar.y
  df_tv <- reshape(pred.data.astar.y, direction = "long", varying = vector_time_points,
                   sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_unique')
  df_tv <- df_tv[order(df_tv$id_unique),]
  df_pred.data.astar.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
  df_pred.data.astar.y <- model.matrix(~.,data=df_pred.data.astar.y)[,-1]
  
  # pred.data.astar.a.astar.a.y
  df_tv <- reshape(pred.data.astar.a.astar.a.y, direction = "long", varying = vector_time_points,
                   sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_unique')
  df_tv <- df_tv[order(df_tv$id_unique),]
  df_pred.data.astar.a.astar.a.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
  df_pred.data.astar.a.astar.a.y <- model.matrix(~.,data=df_pred.data.astar.a.astar.a.y)[,-1]
  
  # pred.data.astar.a.y
  df_tv <- reshape(pred.data.astar.a.y, direction = "long", varying = vector_time_points,
                   sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_unique')
  df_tv <- df_tv[order(df_tv$id_unique),]
  df_pred.data.astar.a.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
  df_pred.data.astar.a.y <- model.matrix(~.,data=df_pred.data.astar.a.y)[,-1]
  
  #######################
  
  PredictY_DEIEM <- predict.aalen(object=Y, Z=df_pred.data.astar.a.y, times=t, n.sim=0)
  PredictY_TEDE_2 <- predict.aalen(object=Y, Z=df_pred.data.a.y, times=t, n.sim=0)
  PredictY_IEMIED <- predict.aalen(object=Y, Z=df_pred.data.astar.a.astar.a.y, times=t, n.sim=0)
  PredictY_IEDTE_1 <- predict.aalen(object=Y, Z=df_pred.data.astar.y, times=t, n.sim=0)
  
  DE <- mean(PredictY_DEIEM$S0) - mean(PredictY_TEDE_2$S0)*100
  IEM <- mean(PredictY_IEDTE_1$S0) - mean(PredictY_IEMIED$S0)*100
  IED <- mean(PredictY_IEMIED$S0) - mean(PredictY_DEIEM$S0)*100
  TE <- mean(PredictY_IEDTE_1$S0) - mean(PredictY_TEDE_2$S0)*100
  effects <- cbind(DE, IEM, IED, TE)
  return(effects)
}


