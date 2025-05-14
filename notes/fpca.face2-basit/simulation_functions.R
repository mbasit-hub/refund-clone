sim_single_rep <- function(seed = 0, K, n, J, sigma,
                           equally, design, level){
  set.seed(seed)
  Data = GeneData(n=n, J=J, sigma=sigma, equally=equally,
                  desig=design, level=level)
  Y =  Data$Y
  scores_TRUE = Data$scores
  evalues_TRUE = Data$evalues
  efunctions_TRUE  = Data$efunctions
  argvals = Data$argvals
  
  # Fit
  start <- Sys.time()
  fit1 <- fpca.face(Y, center = TRUE, argvals=argvals, pve=0.99)
  end <- Sys.time()
  time1 = end - start 
  start <- Sys.time()
  fit2 <- fpca.face2(Y, center = TRUE, argvals=argvals, pve=0.99)
  end <- Sys.time()
  time2 = end - start 
  
  #print(sum((fit2$ret$Yhat - predict.fpca.face2(fit2,Y)$Yhat)^2))
  # MISE of eigenfucntions
  MISE_eigen1 <- sum(unlist(lapply(1:K, function(x){
    min(sum((efunctions_TRUE[,x] - sqrt(J)*fit1$efunctions[,x])^2),
        sum((efunctions_TRUE[,x] + sqrt(J)*fit1$efunctions[,x])^2))})))/(K*J)
  MISE_eigen2 <- sum(unlist(lapply(1:K, function(x){
    min(sum((efunctions_TRUE[,x] - fit2$ret$efunctions[,x])^2),
        sum((efunctions_TRUE[,x] + fit2$ret$efunctions[,x])^2))})))/(K*J)
  
  # MSE for  scores
  MSE_score1 <- sum(unlist(lapply(1:K, function(x){
    min(mean((scores_TRUE[,x] - fit1$scores[,x]/sqrt(J))^2),
        mean((scores_TRUE[,x] + fit1$scores[,x]/sqrt(J))^2))})))/K
  MSE_score2 <- sum(unlist(lapply(1:K, function(x){
    min(mean((scores_TRUE[,x] - fit2$ret$scores[,x])^2),
        mean((scores_TRUE[,x] + fit2$ret$scores[,x])^2))})))/K
  
  
  return(
    list(
      time_old = time1,
      MISEeigen_old = MISE_eigen1,
      MSEscore_old = MSE_score1,
      
      time_new = time2,
      MISEeigen_new = MISE_eigen2,
      MSEscore_new = MSE_score2
    )
  )  
}



sim_multiple_rep <- function(Nel, K, n, J, sigma = 0.5, equally = TRUE,
                             design = "regular", level, parallel = FALSE){
  no.cores <- parallel::detectCores()
  if (parallel) {
    out <- parallel::mclapply(1:Nel, function(iter) {
      res <- tryCatch(
        sim_single_rep(iter, K, n, J, sigma, equally, design, level),
        error = function(e) {
          print(e)
        }
      )
      unlist(res)
    }, mc.cores = no.cores)
  } else {
    out <- lapply(1:Nel, function(iter) {
      res <- tryCatch(
        sim_single_rep(iter, K, n, J, sigma, equally, design, level),
        error = function(e) {
          print(e)
        }
      )
      unlist(res)
    })
  }
  out <- do.call(rbind, out)
  
  return(colMeans(out))
}



sim_final_results <- function(Nel, K, n, J, sigma,
                             equally, design, level, parallel = FALSE){

  row_names <-
    c(
      "J", "n", "design", "time_old", "MISEeigen_old", "MSEscore_old", "time_new",
      "MISEeigen_new", "MSEscore_new"
    )
  
  file_name <- paste0("sim/simulation_results_equal_", equally, ".csv")
  
  cat(row_names, "\n", file = file_name, sep = ",")
  
  out <-
    Map(f = function(y, x, z) {
      result <- sim_multiple_rep(
        Nel = Nel, K = K, n = x, J = y, sigma = sigma, equally = equally,
        design = z, level = level, parallel = parallel
      )
      
      table_output <- c(x, y, z, result)
      cat(table_output, "\n", file = file_name, sep = ",", append = T)
    }, x = n, y = J, z = design)
  
}  
