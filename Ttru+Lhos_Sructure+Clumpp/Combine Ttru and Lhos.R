combine.Ttru.Lhos.data <- function(Ttru.msats,Ttru.strata,Lhos.msats,Lhos.strata){
  
  #merge strata
  diff.cols <- setdiff(colnames(Ttru.strata),colnames(Lhos.strata))
  na.mat <- matrix(nrow=nrow(Lhos.strata),ncol=length(diff.cols))
  colnames(na.mat) <- diff.cols
  Lhos.strata <- cbind(Lhos.strata,na.mat)
  strat.cols <- intersect(colnames(Ttru.strata),colnames(Lhos.strata))
  strata <- rbind(Ttru.strata[,strat.cols],Lhos.strata[,strat.cols])
  colnames(strata)[1] <- "LabID"
#  write.csv(strata,"AS183 AS179 strata.csv", row.names=FALSE)
  
  #merge.msats
  msats <- rbind(Ttru.msats[,colnames(Ttru.msats)],Lhos.msats[,colnames(Ttru.msats)])
  colnames(msats)[1] <- "LabID"
#  msat.merge <- merge(strata, msats, by = "LabID", all.y = TRUE)
#  msat.merge <- msat.merge[order(as.numeric(msat.merge$LABID)), ]
#  locus.col <- ncol(strata)+1
#  write.csv(msats,"AS183 AS179 msats.csv",row.names=FALSE)

  return(list(strata=strata, msats=msats))
}

