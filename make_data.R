set.seed(123)

for(N in SampleSize)
{
  yMat = matrix(quantile.function(runif(N*Nreps),true.par), nrow = N)
  
  saveRDS(yMat, file =  paste("sample_",mc.name,"_N",N, ".RDS",sep=""))
}