plot_ts_fideos <- function(data=syr,data_ev=NULL,ylabel="Spawning biomass (t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=500,breaks_x =c(2020,2025,2030),yref=Srms,nview=6,seed=NULL){
  if(is.null(seed)) seed <- 356
  set.seed(seed)
  dfp <- data
  Fstr <- rep(Str,nsim)
  dfp <- as.data.frame(dfp)
  colnames(dfp) <- yrs
  dfp <- cbind(dfp,Str=Fstr)
  dfe <- as.data.frame(data_ev)
  colnames(dfe) <- yr
  df <- NULL
  for(i in 1:nstr){
    tmp_dfp <- dfp[dfp$Str==Str[i],,] 
    tmp_df <- cbind(dfe,tmp_dfp)
    df <- rbind(df,tmp_df)
  }
  selec_sim <- NULL
  for(i in 1:nstr){
    ind1 <- nsim*i-nsim+1
    ind2 <- nsim*i
    tmp_sel <- floor(runif(1,min=ind1,max=ind2))
    selec_sim <- c(selec_sim,tmp_sel)
  }
  #selec_sim <- floor(runif(nview,min=1,max=nsim))
  #selec_sim <- c(3,505,1015,)
  df2 <- df[selec_sim,]
  df2 <- melt(df2)
  colnames(df2) <- c("u","Year","Var")
  
  p <- ggplot(data=df2)+
    #geom_line(aes(x=Year,y=Var,group=u,linetype=u))+
    geom_line(aes(x=Year,y=Var,group=u,col=u))+
    xlab("Year")+
    ylab(ylabel)+
    geom_hline(yintercept=yref,linetype="dashed")+
    #theme(axis.title.y = element_text(size = rel(0.1), angle = 90))+
    #theme(axis.title.x = element_text(size = rel(0.1), angle = 00))+
    #theme(legend.position = "none")+
    scale_x_discrete(breaks = breaks_x)+mi.tema()+
    theme(legend.position = "none")
  p
}
