plot_ts <- function(data=syr,data_ev=NULL,data_es=NULL,ylabel="Spawning biomass (t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=2,ytick=0.4,yref=Target){
  dfp <- data
  Fstr <- rep(Str,nsim)
  dfp <- as.data.frame(dfp)
  dfes <- data_es
  dfes <- as.data.frame(dfes)
  colnames(dfp) <- colnames(dfes) <- yrs
  dfp <- cbind(dfp,Str=Fstr)
  dfe <- as.data.frame(data_ev)
  dfes <- cbind(dfes,Str=Fstr)
  colnames(dfe) <- yr
  df <- NULL
  for(i in 1:nstr){
    tmp_dfp <- dfp[dfp$Str==Str[i],,] 
    tmp_df <- cbind(dfe,tmp_dfp)
    tmp_df <- t(apply(tmp_df[,1:(dim(tmp_df)[2]-1)],2,quantile,c(0.025,0.1,0.5,0.9,0.975),na.rm=TRUE))
    df <- rbind(df,tmp_df)
  }
  colnames(df) <- c("Li","L1","Mediana","L2","Ls")
  iyr <- c(yr,yrs)
  #hcr <- c(rep(Str[1],length(iyr)),rep(Str[2],length(iyr)),rep(Str[3],length(iyr)),rep(Str[4],length(iyr)),rep(Str[5],length(iyr)),rep(Str[6],length(iyr)))
  iyr2 <- rep(iyr,nstr)
  hcr <- NULL
  for(i in 1:nstr){
    tmp_hcr <- rep(Str[i],length(iyr))
    hcr <- c(hcr,tmp_hcr)
  }
  df <- data.frame(df,Year=iyr2,Str=hcr)
  df2 <- NULL
  for(i in 1:nstr){
    tmp_dfp <- dfes[dfes$Str==Str[i],,] 
    tmp_df <- cbind(dfe,tmp_dfp)
    tmp_df <- t(apply(tmp_df[,1:(dim(tmp_df)[2]-1)],2,quantile,c(0.025,0.1,0.5,0.9,0.975),na.rm=TRUE))
    df2 <- rbind(df2,tmp_df)
  }
  colnames(df2) <- c("Li","L1","Mediana","L2","Ls")
  hcr <- NULL
  for(i in 1:nstr){
    tmp_hcr <- rep(Str[i],length(iyr))
    hcr <- c(hcr,tmp_hcr)
  }
  df2 <- data.frame(df2,Year=iyr2,Str=hcr)
  
    p <- ggplot(data=df,aes(x=Year,y=Mediana,group=Str))+
      geom_ribbon(aes(ymin=Li,ymax=Ls),fill="grey70")+
      geom_ribbon(aes(ymin=L1,ymax=L2),fill="grey30")+
      geom_point(data=df2,aes(x=Year,y=Mediana),col="black",size=0.6)+
      geom_rect(xmin = yr[1],
                xmax = yr[length(yr)],
                ymin = ymin, ymax = ymax,
                fill = "grey", alpha = 0.01)+
      geom_line()+facet_wrap(~Str,ncol=nstr)+
      geom_hline(yintercept=yref,linetype="dashed")+
      #geom_hline(yintercept=yref/2,linetype="dashed")+
      scale_y_continuous(name=ylabel,breaks=seq(ymin,ymax,by=ytick),limits = c(ymin,ymax))+mi.tema()
    p
}
