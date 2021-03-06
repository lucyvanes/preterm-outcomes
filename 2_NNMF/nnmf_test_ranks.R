
nnmf_test_ranks <- function(rank, data, WH_true){
  
  # input:
  # rank = rank at which to be run
  # data = input data matrix
  # WH_true = should wold-holdouts be applied (i.e. mask out 20% of matrix elements)?
  
  library(NNLM)
  library(Metrics)
  
  # Run NMF
  #===================================
  print(paste("running single nnmf with rank", rank), quote=F)
  print("=====================================================", quote=F)
  
  res <- nnmf(data, rank, method="lee")
  
  w <- res$W
  h <- res$H
  V.hat <- w %*% h
  
  ## RECONSTRUCTION ERROR
  #=============================
  # frobenius norm
  recon_error_frob <- norm(data, type="F") - norm(V.hat, type="F")
  
  # write data
  write.table(recon_error_frob, paste("mclapply_recon_error_frob_",rank,".txt", sep=""), row.names=F, quote=F)
  
  if (WH_true==1){
    # Wold holdouts
    #================================
    wold_holdouts <- data.frame("rep"=1:50, "recon_error_frob"=NA, "rmse_train"=NA, "rmse_test"=NA,
                                "mse_train"=NA, "mse_test"=NA)
    for (r in 1:50){
      data_holdout <- data
      index <- sample(length(data), length(data)*0.2);
      data_holdout[index] <- NA;
      
      print(paste("running CV for rank",rank, "with 20% Wold holdouts - fold", r, sep=" "))
      res_nnlm <- nnmf(data_holdout, rank, method="lee")
      
      V.hat <- with(res_nnlm, W %*% H)
      
      # reconstruction error - whole matrix
      wold_holdouts$recon_error_frob[wold_holdouts$rep==r] <- norm(data, type="F") - norm(V.hat, type="F")
      
      # reconstruction error - included values
      not_index <- setdiff(1:length(data), index)
      wold_holdouts$rmse_train[wold_holdouts$rep==r] <- rmse(V.hat[not_index], data[not_index])
      
      # reconstruction error - imputed values 
      wold_holdouts$rmse_test[wold_holdouts$rep==r] <- rmse(V.hat[index], data[index])
    }
    
    mean_WH_recon_error_frob <- mean(wold_holdouts$recon_error_frob)
    RMSE_train <- mean(wold_holdouts$rmse_train)
    RMSE_test <- mean(wold_holdouts$rmse_test)
    
    write.table(mean_WH_recon_error_frob, paste("mean_WH_recon_error_frob_",rank,".txt", sep=""), row.names=F, quote=F)
    write.table(RMSE_train, paste("RMSE_train_",rank,".txt", sep=""), row.names=F, quote=F)
    write.table(RMSE_test, paste("RMSE_test_",rank,".txt", sep=""), row.names=F, quote=F)
  
    list(res=res, 
         recon_error_frob=recon_error_frob, 
         wold_holdouts=wold_holdouts, mean_WH_recon_error_frob=mean_WH_recon_error_frob,
         RMSE_train=RMSE_train, RMSE_test=RMSE_test)
    
  } else if (WH_true==0){
    list(res=res, recon_error_frob=recon_error_frob)
  }
}

