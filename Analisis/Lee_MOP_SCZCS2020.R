# Lee Modelos Operativos de Anchoveta Zona Centro Norte

md <- read.admb(paste0("Sczcs2020/OP4/","OP4sczcs")) #Lee el reporte y resultados de ajuste
names(md)
yr <- md$years
md$Years <- md$years
# Lee los datos de los MOP y los prepara para graficar --------------------

nstr=1
nsim= 200
nr=nstr*nsim #numero de filas
So= md$BD_virgen_LP/1000 #Biomasa virginal
Srms <- 0.55*So


#### MODELO B/B0 dinámico
##### MODELO OPERATIVO 4
carpeta <- 'Sczcs2020/OP1/'
ceval = matrix(md$desembarqueobs,nrow=200,ncol=length(yr),byrow = TRUE)/1000
#Biomasa
byr=as.matrix(read.table(paste0(carpeta,"02BiomasaTotal_op.txt"),nrow=nr,fill=T))/1000
beval=as.matrix(read.table(paste0(carpeta,"09BiomTotal_hist.txt"),nrow=nr,fill=T))/1000
best =as.matrix(read.table(paste0(carpeta,"02BiomasaTotal_est.txt"),nrow=nr,fill=T))/1000
#Biomasa desovante
syr=as.matrix(read.table(paste0(carpeta,"03Desovante_op.txt"),nrow=nr,fill=T))/1000
seval=as.matrix(read.table(paste0(carpeta,"08Desovante_hist.txt"),nrow=nr,fill=T))/1000
sest =as.matrix(read.table(paste0(carpeta,"03Desovante_est.txt"),nrow=nr,fill=T))/1000
#Reclutamiento
ryr=as.matrix(read.table(paste0(carpeta,"04Reclutamiento_op.txt"),nrow=nr,fill=T))/1000
reval=as.matrix(read.table(paste0(carpeta,"10Reclutas_hist.txt"),nrow=nr,fill=T))/1000
rest =as.matrix(read.table(paste0(carpeta,"04Reclutamiento_est.txt"),nrow=nr,fill=T))/1000
#Fmort
fyr=as.matrix(read.table(paste0(carpeta,"05FMort_op.txt"),nrow=nr,fill=T))
feval=as.matrix(read.table(paste0(carpeta,"07FMort_hist.txt"),nrow=nr,fill=T))
fest =as.matrix(read.table(paste0(carpeta,"05FMort_est.txt"),nrow=nr,fill=T))
#Agotamiento
dpyr  =as.matrix(read.table(paste0(carpeta,"12RPR_op.txt"),nrow=nr,fill=T))
dpeval=as.matrix(read.table(paste0(carpeta,"12RPR_hist.txt"),nrow=nr,fill=T))
dpest =as.matrix(read.table(paste0(carpeta,"11RPR_est.txt"),nrow=nr,fill=T))
#Capturas
cyr=as.matrix(read.table(paste0(carpeta,"01Capturas_proyectadas.txt"),nrow=nr,fill=T))/1000


Target = 0.55
Rprom <- mean(md$Reclutas)/1000
range(ryr)
Str <- "MOP1"
yrs <- seq(range(md$Years)[2]+1,range(md$Years)[2]+20,1)
op1.f1 <- plot_ts(data=ryr,data_ev=reval,data_es=rest,ylabel="Reclut. (millones)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=1200,ytick=100,yref =Rprom)
op1.f1
range(byr)
op1.f2 <- plot_ts(data = byr,data_ev = beval,data_es = best,ylabel="Biomasa (miles t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=7000,ytick=500,yref=NULL)
op1.f2

op1.f3 <- plot_ts(data = syr,data_ev = seval,data_es = sest,ylabel="Biom. Desov.(miles t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=4000,ytick=500,yref=Srms)
op1.f3

op1.f4 <- plot_ts(data=cyr,data_ev=ceval,data_es=cyr,ylabel="Captura (miles t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=1200,ytick=200,yref = NULL)
op1.f4

op1.f5 <- plot_ts(data = dpyr,data_ev = dpeval,data_es = dpest,ylabel="BD/Bo",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=2.5,ytick=0.5,yref=Target)
op1.f5

Target=0.3
op1.f6 <- plot_ts(data = fyr,data_ev = feval,data_es = fest,ylabel="Ft",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=1.5,ytick=0.5,yref=Target)
op1.f6



##### MODELO OPERATIVO 2
carpeta <- 'Sczcs2020/OP2/'
#Biomasa total
byr=as.matrix(read.table(paste0(carpeta,"02BiomasaTotal_op.txt"),nrow=nr,fill=T))/1000
beval=as.matrix(read.table(paste0(carpeta,"09BiomTotal_hist.txt"),nrow=nr,fill=T))/1000
best =as.matrix(read.table(paste0(carpeta,"02BiomasaTotal_est.txt"),nrow=nr,fill=T))/1000
#Biomasa desovante
syr=as.matrix(read.table(paste0(carpeta,"03Desovante_op.txt"),nrow=nr,fill=T))/1000
seval=as.matrix(read.table(paste0(carpeta,"08Desovante_hist.txt"),nrow=nr,fill=T))/1000
sest =as.matrix(read.table(paste0(carpeta,"03Desovante_est.txt"),nrow=nr,fill=T))/1000
#Reclutamiento
ryr=as.matrix(read.table(paste0(carpeta,"04Reclutamiento_op.txt"),nrow=nr,fill=T))/1000
reval=as.matrix(read.table(paste0(carpeta,"10Reclutas_hist.txt"),nrow=nr,fill=T))/1000
rest =as.matrix(read.table(paste0(carpeta,"04Reclutamiento_est.txt"),nrow=nr,fill=T))/1000
#Fmort
fyr=as.matrix(read.table(paste0(carpeta,"05FMort_op.txt"),nrow=nr,fill=T))
feval=as.matrix(read.table(paste0(carpeta,"07FMort_hist.txt"),nrow=nr,fill=T))
fest =as.matrix(read.table(paste0(carpeta,"05FMort_est.txt"),nrow=nr,fill=T))
#Agotamiento
dpyr  =as.matrix(read.table(paste0(carpeta,"12RPR_op.txt"),nrow=nr,fill=T))
dpeval=as.matrix(read.table(paste0(carpeta,"12RPR_hist.txt"),nrow=nr,fill=T))
dpest =as.matrix(read.table(paste0(carpeta,"11RPR_est.txt"),nrow=nr,fill=T))
#Capturas
cyr=as.matrix(read.table(paste0(carpeta,"01Capturas_proyectadas.txt"),nrow=nr,fill=T))/1000


Target = 0.55
Rprom <- mean(md$Reclutas)/1000
range(ryr)
Str <- "MOP2"
yrs <- seq(range(md$Years)[2]+1,range(md$Years)[2]+20,1)
op2.f1 <- plot_ts(data=ryr,data_ev=reval,data_es=rest,ylabel="Reclut. (millones)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=1200,ytick=100,yref =Rprom)
op2.f1
range(byr)
op2.f2 <- plot_ts(data = byr,data_ev = beval,data_es = best,ylabel="Biomasa (miles t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=7000,ytick=500,yref=NULL)
op2.f2

op2.f3 <- plot_ts(data = syr,data_ev = seval,data_es = sest,ylabel="Biom. Desov.(miles t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=4000,ytick=500,yref=Srms)
op2.f3

op2.f4 <- plot_ts(data=cyr,data_ev=ceval,data_es=cyr,ylabel="Captura (miles t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=1200,ytick=200,yref = NULL)
op2.f4

op2.f5 <- plot_ts(data = dpyr,data_ev = dpeval,data_es = dpest,ylabel="BD/Bo",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=2.5,ytick=0.5,yref=Target)
op2.f5

Target=0.3
op2.f6 <- plot_ts(data = fyr,data_ev = feval,data_es = fest,ylabel="Ft",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=1.5,ytick=0.5,yref=Target)
op2.f6

##### MODELO OPERATIVO 3
carpeta <- 'Sczcs2020/OP3/'
#Biomasa total
byr=as.matrix(read.table(paste0(carpeta,"02BiomasaTotal_op.txt"),nrow=nr,fill=T))/1000
beval=as.matrix(read.table(paste0(carpeta,"09BiomTotal_hist.txt"),nrow=nr,fill=T))/1000
best =as.matrix(read.table(paste0(carpeta,"02BiomasaTotal_est.txt"),nrow=nr,fill=T))/1000
#Biomasa desovante
syr=as.matrix(read.table(paste0(carpeta,"03Desovante_op.txt"),nrow=nr,fill=T))/1000
seval=as.matrix(read.table(paste0(carpeta,"08Desovante_hist.txt"),nrow=nr,fill=T))/1000
sest =as.matrix(read.table(paste0(carpeta,"03Desovante_est.txt"),nrow=nr,fill=T))/1000
#Reclutamiento
ryr=as.matrix(read.table(paste0(carpeta,"04Reclutamiento_op.txt"),nrow=nr,fill=T))/1000
reval=as.matrix(read.table(paste0(carpeta,"10Reclutas_hist.txt"),nrow=nr,fill=T))/1000
rest =as.matrix(read.table(paste0(carpeta,"04Reclutamiento_est.txt"),nrow=nr,fill=T))/1000
#Fmort
fyr=as.matrix(read.table(paste0(carpeta,"05FMort_op.txt"),nrow=nr,fill=T))
feval=as.matrix(read.table(paste0(carpeta,"07FMort_hist.txt"),nrow=nr,fill=T))
fest =as.matrix(read.table(paste0(carpeta,"05FMort_est.txt"),nrow=nr,fill=T))
#Agotamiento
dpyr  =as.matrix(read.table(paste0(carpeta,"12RPR_op.txt"),nrow=nr,fill=T))
dpeval=as.matrix(read.table(paste0(carpeta,"12RPR_hist.txt"),nrow=nr,fill=T))
dpest =as.matrix(read.table(paste0(carpeta,"11RPR_est.txt"),nrow=nr,fill=T))
#Capturas
cyr=as.matrix(read.table(paste0(carpeta,"01Capturas_proyectadas.txt"),nrow=nr,fill=T))/1000


Target = 0.55
Rprom <- mean(md$Reclutas)/1000
range(ryr)
Str <- "MOP3"
yrs <- seq(range(md$Years)[2]+1,range(md$Years)[2]+20,1)
op3.f1 <- plot_ts(data=ryr,data_ev=reval,data_es=rest,ylabel="Reclut. (millones)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=1200,ytick=100,yref =Rprom)
op3.f1
range(byr)
op3.f2 <- plot_ts(data = byr,data_ev = beval,data_es = best,ylabel="Biomasa (miles t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=7000,ytick=500,yref=NULL)
op3.f2

op3.f3 <- plot_ts(data = syr,data_ev = seval,data_es = sest,ylabel="Biom. Desov.(miles t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=4000,ytick=500,yref=Srms)
op3.f3

op3.f4 <- plot_ts(data=cyr,data_ev=ceval,data_es=cyr,ylabel="Captura (miles t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=1200,ytick=200,yref = NULL)
op3.f4

op3.f5 <- plot_ts(data = dpyr,data_ev = dpeval,data_es = dpest,ylabel="BD/Bo",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=2.5,ytick=0.5,yref=Target)
op3.f5

Target=0.3
op3.f6 <- plot_ts(data = fyr,data_ev = feval,data_es = fest,ylabel="Ft",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=1.5,ytick=0.5,yref=Target)
op3.f6


f01 <- ggarrange(op1.f1,op1.f3,op2.f1,op2.f3,op3.f1,op3.f3,labels=c("A","B","C","D","E","F"),nrow=3,ncol=2)
f01
f02 <- ggarrange(op1.f4,op1.f6,op2.f4,op2.f6,op3.f4,op3.f6,labels=c("A","B","C","D","E","F"),nrow=3,ncol=2)
f02
f03 <- ggarrange(op1.f5,op2.f5,op3.f5,labels=c("A","B","C"),nrow=3,ncol=1)
f03

ggsave("Figs/fig1_rec_ssb_sczcs.jpg",plot=f01,dpi=300,width = 20,height = 28,units = "cm")
ggsave("Figs/fig2_catch_Ft_sczcs.jpg",plot=f02,dpi=300,width = 20,height = 28,units = "cm")
ggsave("Figs/fig3_BBo_sczcs.jpg",plot=f03,dpi=300,width = 14,height = 18,units = "cm")

#### MODELO B/B0 dinámico
##### MODELO OPERATIVO 4
carpeta <- 'Sczcs2020/OP4/'
ceval = matrix(md$desembarqueobs,nrow=200,ncol=length(yr),byrow = TRUE)/1000
#Biomasa total
byr=as.matrix(read.table(paste0(carpeta,"02BiomasaTotal_op.txt"),nrow=nr,fill=T))/1000
beval=as.matrix(read.table(paste0(carpeta,"08BiomTotal_hist.txt"),nrow=nr,fill=T))/1000
best =as.matrix(read.table(paste0(carpeta,"02BiomasaTotal_est.txt"),nrow=nr,fill=T))/1000
#Biomasa desovante
syr=as.matrix(read.table(paste0(carpeta,"03Desovante_op.txt"),nrow=nr,fill=T))/1000
seval=as.matrix(read.table(paste0(carpeta,"09Desovante_hist.txt"),nrow=nr,fill=T))/1000
sest =as.matrix(read.table(paste0(carpeta,"03Desovante_est.txt"),nrow=nr,fill=T))/1000
#Reclutamiento
ryr=as.matrix(read.table(paste0(carpeta,"04Reclutamiento_op.txt"),nrow=nr,fill=T))/1000
reval=as.matrix(read.table(paste0(carpeta,"10Reclutas_hist.txt"),nrow=nr,fill=T))/1000
rest =as.matrix(read.table(paste0(carpeta,"04Reclutamiento_est.txt"),nrow=nr,fill=T))/1000
#Fmort
fyr=as.matrix(read.table(paste0(carpeta,"05FMort_op.txt"),nrow=nr,fill=T))
feval=as.matrix(read.table(paste0(carpeta,"10FMort_hist.txt"),nrow=nr,fill=T))
fest =as.matrix(read.table(paste0(carpeta,"05FMort_est.txt"),nrow=nr,fill=T))
#Agotamiento
dpyr  =as.matrix(read.table(paste0(carpeta,"06RPR_op.txt"),nrow=nr,fill=T))
dpeval=as.matrix(read.table(paste0(carpeta,"11RPR_hist.txt"),nrow=nr,fill=T))
dpest =as.matrix(read.table(paste0(carpeta,"06RPR_est.txt"),nrow=nr,fill=T))
#Agotamiento dinámco
dinpyr  =as.matrix(read.table(paste0(carpeta,"07RPRdin_op.txt"),nrow=nr,fill=T))
dinpeval=as.matrix(read.table(paste0(carpeta,"12RPRdin_hist.txt"),nrow=nr,fill=T))
dinpest =as.matrix(read.table(paste0(carpeta,"07RPRdin_est.txt"),nrow=nr,fill=T))
#Capturas
cyr=as.matrix(read.table(paste0(carpeta,"01Capturas_proyectadas.txt"),nrow=nr,fill=T))/1000

Target = 0.55
Rprom <- mean(md$Reclutas)/1000
range(ryr)
Str <- "MOP4"
yrs <- seq(range(md$Years)[2]+1,range(md$Years)[2]+20,1)
op4.f1s <- plot_ts(data=ryr,data_ev=reval,data_es=rest,ylabel="Reclut. (millones)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=1300,ytick=100,yref =Rprom)
op4.f1s
range(byr)
op4.f2s <- plot_ts(data = byr,data_ev = beval,data_es = best,ylabel="Biomasa (miles t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=9000,ytick=500,yref=NULL)
op4.f2s

op4.f3s <- plot_ts(data = syr,data_ev = seval,data_es = sest,ylabel="Biom. Desov.(miles t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=5000,ytick=500,yref=Srms)
op4.f3s

op4.f4s <- plot_ts(data=cyr,data_ev=ceval,data_es=cyr,ylabel="Captura (miles t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=1200,ytick=200,yref = NULL)
op4.f4

Target=0.3
op4.f5 <- plot_ts(data = fyr,data_ev = feval,data_es = fest,ylabel="Ft",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=1.5,ytick=0.5,yref=Target)
op4.f5

Target=0.55
op4.f6 <- plot_ts(data = dpyr,data_ev = dpeval,data_es = dpest,ylabel="BD/Bo estático",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=2.5,ytick=0.5,yref=Target)
op4.f6

op4.f7 <- plot_ts(data = dinpyr,data_ev = dinpeval,data_es = dinpest,ylabel="BD/Bo dinámico",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ymax=2.5,ytick=0.5,yref=Target)
op4.f7



