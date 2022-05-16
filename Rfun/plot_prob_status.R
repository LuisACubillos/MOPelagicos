plot_prob_status <- function(data=nrpr,depl=0.2,targt=0.4,Str=Str,nstr=nstr,nsim=nsim,breaks_x=c(2020,2025,2030)){
  df_tmp <- data
  Fstr <- rep(Str,nsim)
  dfp_tmp <- as.data.frame(df_tmp)
  colnames(dfp_tmp) <- yrs
  dfp <- cbind(dfp_tmp,Str=Fstr)
  dfd <- NULL
  dft <- NULL
  nstr=nstr
  for(i in 1:nstr){
    tmp_dfd <- dfp[dfp$Str==Str[i],,]
    tmp_dfd <- tmp_dfd[,1:length(yrs)]
    tmp_dfd <- t(apply(tmp_dfd, 2,function(x) sum(ifelse(x<depl,1,0))/nsim))
    dfd <- rbind(dfd,tmp_dfd)
    tmp_dft <- dfp[dfp$Str==Str[i],,]
    tmp_dft <- tmp_dft[,1:length(yrs)]
    tmp_dft <- t(apply(tmp_dft, 2,function(x) sum(ifelse(x>targt,1,0))/nsim))
    dft <- rbind(dft,tmp_dft)
  }
  dat_prob <- data.frame(dfd)
  colnames(dat_prob) <- yrs
  dat_prob$u <- Str
  dat_prob_depletion <- melt(dat_prob)
  p1 <- ggplot(dat_prob_depletion)+xlab("Year")+
    #geom_line(aes(x=variable,y=value,group=u,linetype=u))+
    geom_line(aes(x=variable,y=value,group=u,col=u))+
    mi.tema()+theme(legend.position = c(0.8,0.8))+
    scale_x_discrete(breaks = breaks_x)+
    scale_y_continuous(name="Probability",breaks=seq(0,1,by=0.2),limits = c(0,1))
  
  dat_prob <- data.frame(dft)
  colnames(dat_prob) <- yrs
  dat_prob$u <- Str
  dat_prob_target <- melt(dat_prob)
  p2 <- ggplot(dat_prob_target)+xlab("Year")+
    #geom_line(aes(x=variable,y=value,group=u,linetype=u))+
    geom_line(aes(x=variable,y=value,group=u,col=u))+
    mi.tema()+theme(legend.position = "none")+
    scale_x_discrete(breaks = breaks_x)+
    scale_y_continuous(name="Probability",breaks=seq(0,1,by=0.2),limits = c(0,1))
  
  p3 <- ggarrange(p1,p2,labels = c("A","B"),nrow=1,ncol = 2)
}
