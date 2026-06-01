get_chai_seed <- function(default = 123L) {
  option_seed <- getOption("chai.seed", default)
  seed_value <- suppressWarnings(as.integer(option_seed))
  if (is.na(seed_value)) {
    default
  } else {
    seed_value
  }
}


# Fit Mclust with fallback
fit_mclust_with_fallback <- function(df, K_vec = 2:10, timeout_sec = 30, jitter_sd = 1e-6) {
  fit <- NULL
  fit_error <- NULL

  if (.Platform$OS.type == "unix" && is.finite(timeout_sec) && timeout_sec > 0) {
    job <- parallel::mcparallel(mclust::Mclust(df, G = K_vec), silent = TRUE)
    collected <- parallel::mccollect(job, wait = FALSE, timeout = timeout_sec)

    if (length(collected)) {
      fit_candidate <- collected[[1]]
      if (inherits(fit_candidate, "try-error") || inherits(fit_candidate, "error")) {
        fit_error <- simpleError(as.character(fit_candidate))
      } else {
        fit <- fit_candidate
      }
    } else {
      tools::pskill(job$pid, tools::SIGKILL)
      parallel::mccollect(job, wait = FALSE)
      fit_error <- simpleError(sprintf("Mclust exceeded %s seconds.", timeout_sec))
    }
  } else {
    fit_candidate <- tryCatch(mclust::Mclust(df, G = K_vec), error = function(e) e)
    if (inherits(fit_candidate, "error")) {
      fit_error <- fit_candidate
    } else {
      fit <- fit_candidate
    }
  }

  if (!is.null(fit)) {
    return(fit)
  }

  warning(
    "Primary Mclust fit failed or timed out; retrying with scaled+jittered input. Error: ",
    conditionMessage(fit_error)
  )

  df_stable <- as.data.frame(lapply(df, function(col) {
    scaled <- as.numeric(scale(col))
    scaled[!is.finite(scaled)] <- 0
    scaled
  }))

  set.seed(get_chai_seed())
  for (col_name in names(df_stable)) {
    df_stable[[col_name]] <- df_stable[[col_name]] + stats::rnorm(nrow(df_stable), sd = jitter_sd)
  }

  mclust::Mclust(df_stable, G = K_vec)
}

# This FDP is without cutoff
computeFDP<-function(lFDR,trueIndex,alpha) {
  df<-data.frame(cbind(1:length(lFDR),lFDR))
  df<-df[order(df[,2]),]
  df$cumFDR<-cumsum(df[,2])
  df$cumFDR<-df$cumFDR/(1:nrow(df))
  df<-df[order(df[,1]),]
  select<-which(df$cumFDR<=alpha)
  if(length(select)==0){
    FDP<-0
  }else{
    FDP<-1-length(intersect(select,trueIndex))/length(select)
  }
  FDP
}


computePower<-function(lFDR,trueIndex,alpha) {
  df<-data.frame(cbind(1:length(lFDR),lFDR))
  df<-df[order(df[,2]),]
  df$cumFDR<-cumsum(df[,2])
  df$cumFDR<-df$cumFDR/(1:nrow(df))
  df<-df[order(df[,1]),]
  select<-which(df$cumFDR<=alpha)
  length(intersect(select,trueIndex))/length(trueIndex)
}

computeFDPIndex<-function(index,trueIndex) {
  if(length(index)==0){
    FDP<-0
  }else{
    FDP<-1-length(intersect(index,trueIndex))/length(index)
  }
  FDP
}


computePowerIndex<-function(index,trueIndex) {
  length(intersect(index,trueIndex))/length(trueIndex)
}


ratio_z_given_x_custom<- function(z, post_w, m, v) {
  log_pz <- log(
    sum(post_w * dnorm(z, mean=m, sd=sqrt(v), log=FALSE))
  )
  log_phi_z <- dnorm(z, mean=0, sd=1, log=TRUE)
  # log ratio
  log_ratio <- exp(log_pz - log_phi_z)
  return(log_ratio)  # still can be a big or small #, but less likely Inf
}

minRatioForX_custom <- function(post_w, m, v, lower=-10, upper=10) {
  f <- function(z) ratio_z_given_x_custom(z, post_w, m, v)
  res <- optimize(f, interval=c(lower, upper), maximum=FALSE)
  list(z_star = res$minimum, ratio_star = res$objective)
}

