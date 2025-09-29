library(SimDesign)
library(timereg)
library(tidyr)
library(dplyr)
library(mice)
library(splitstackshape)
library(foreach)
# Only if you want to run parallel models
#library(doParallel)

#num_cores <- parallel::detectCores() - 2
#cl <- makeCluster(num_cores)
#registerDoParallel(cl)
# Ensure the cluster is stopped when the function/script exits
#on.exit(stopCluster(cl))



load("simData.mimick.SHS.RDa")


set.seed(2)  
# subset the data into 10000 instead of 2000. 2025/8/26
simData <- simData[sample(1:nrow(simData), size = 10000, replace = FALSE), ]



# ---- Start: data prep for simData ----
df <- as.data.frame(simData)

# 1) Rename & add required columns (wide format expected)
df$idno <- df$id
df$M_1  <- df$M1
df$M_2  <- df$M2

# time index helpers required by your function’s reshape

# simulated time: 4 and 8 years after
df$time.since.first.exam_1 <- 4
df$time.since.first.exam_2 <- 8

# coerce types if needed
to_num <- c("E","C","M_1","M_2","time_to_event")
for (cc in intersect(to_num, names(df))) df[[cc]] <- suppressWarnings(as.numeric(as.character(df[[cc]])))
df$eventHappened <- as.integer(df$eventHappened > 0)
if ("D1" %in% names(df)) df$D1 <- as.integer(df$D1 > 0)
if ("D2" %in% names(df)) df$D2 <- as.integer(df$D2 > 0)

# keep only rows complete for BOTH mediators and predictors they use
keep <- complete.cases(df[, c("E","C","M_1","M_2")])
df2  <- df[keep, , drop = FALSE]   # this will be ~1795 rows 

cat("Rows kept:", nrow(df2), "\n")  #  should print 1795




df2 <- within(df, {
  idno <- id
  M_1  <- M1
  M_2  <- M2
  time.since.first.exam_1 <- 1
  time.since.first.exam_2 <- 2
})
df2 <- df2[complete.cases(df2[, c("M_1","M_2")]), ]

m <- c("M_1","M_2")
time_points <- c("1","2")

control.value <- as.numeric(quantile(df2$E, 0.25, na.rm = TRUE))
treat.value   <- as.numeric(quantile(df2$E, 0.75, na.rm = TRUE))
treat <- "E"

M <- list(
  lm(M_1 ~ E + C,           data = df2),
  lm(M_2 ~ E + M_1 + C,     data = df2)
)

L <- list(
  glm(D1 ~ E + M_1 + C, data = df2, family = binomial()),
  glm(D2 ~ E + M_2 + C, data = df2, family = binomial())
)

base <- df2[, c("idno","E","C","M_1","M_2","time_to_event","eventHappened")]

# Build a long(person-time) data with exactly two rows per subject
longA <- transform(base,
                   tstart = 0, tstop = pmin(1, time_to_event),
                   endpt  = as.integer(eventHappened == 1 & time_to_event <= 1),
                   M      = M_1)
longB <- transform(base,
                   tstart = 1, tstop = time_to_event,
                   endpt  = as.integer(eventHappened == 1 & time_to_event > 1),
                   M      = M_2)

## Keep valid intervals only
longA <- longA[longA$tstop > longA$tstart, ]
longB <- longB[longB$tstop > longB$tstart, ]
df_long <- rbind(longA, longB)

## Fit additive hazards model
Y <- aalen(Surv(tstart, tstop, endpt) ~ const(E) + const(M) + const(C),
           data = df_long, n.sim = 0, resample.iid = 1)

stopifnot(all(c("idno", treat, m,
                "time.since.first.exam_1","time.since.first.exam_2")
              %in% names(df2)))

res <- med_longitudinal(
  L = L,             # or NULL if you don't model dropout
  M = M,
  m = m,
  Y = Y,
  treat = treat,
  control.value = control.value,
  treat.value   = treat.value,
  data = df2,
  time_points = time_points,
  peryr = 100000
)
print(res)