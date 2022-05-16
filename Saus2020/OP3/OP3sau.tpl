//========================================================================

GLOBALS_SECTION

//========================================================================
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc.csv");
 #undef rep
 #define rep(object) report << #object "\n" << object << endl; // "\n" = FIN DE LINEA;
 
 #undef R_Report
 #define R_Report(object) R_report << #object "\n" << object << endl;

 adstring simname;
 //Nombre de archivo gradiente
 ofstream rep_grad("00rep_convergencia.txt",ios::app);
 //Nombre de archivos para las variables de estado del simulador
 ofstream rep1("01Capturas_proyectadas.txt",ios::app);
 ofstream rep2("02BiomasaTotal_op.txt",ios::app);
 ofstream rep3("03Desovante_op.txt",ios::app);
 ofstream rep4("04Reclutamiento_op.txt",ios::app);
 ofstream rep5("05FMort_op.txt",ios::app);
 //Nombre de archivos para las variables de estado del estimador
 ofstream rep6("02BiomasaTotal_est.txt",ios::app);
 ofstream rep7("03Desovante_est.txt",ios::app);
 ofstream rep8("04Reclutamiento_est.txt",ios::app);
 ofstream rep9("05FMort_est.txt",ios::app);
 //Nombre de archivos para fase histórica del simulador agregado por LCubillos
 ofstream rep10("07Fmort_hist.txt",ios::app);
 ofstream rep11("08Desovante_hist.txt",ios::app);
 ofstream rep12("09BiomTotal_hist.txt",ios::app);
 ofstream rep13("10Reclutas_hist.txt",ios::app);   
 ofstream rep14("11RPR_est.txt",ios::app);   
 ofstream rep15("12RPR_op.txt",ios::app);  
 ofstream rep16("12RPR_hist.txt",ios::app);  
  
 
//========================================================================

TOP_OF_MAIN_SECTION

//========================================================================
 time(&start);
 arrmblsize = 90000000; 
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7); 
 gradient_structure::set_MAX_NVAR_OFFSET(5000); 
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000); 

//========================================================================

DATA_SECTION
  //Modificaciones por LA Cubillos y MJ Cuevas
  int iseed
  !!long int lseed=iseed;
  !!CLASS random_number_generator rng(iseed);
 
//========================================================================
 init_int ntime  
 init_int nedades
 init_int ntallas
 init_matrix mdatos(1,ntime,1,13)
 init_vector Tallas(1,ntallas)
 init_matrix Ctot(1,ntime,1,ntallas)
 init_matrix Ncru(1,ntime,1,ntallas)

 
 !!ad_comm::change_datafile_name("OPsau_control.txt");
 init_vector msex(1,ntallas)
 init_vector Wmed(1,ntallas)
 init_number sigmaR
 init_vector dt(1,3)
 init_vector Par_bio(1,7)
 init_vector cv_Par_bio(1,7)
 init_int    minedad
 init_number bprior //Hiperestabilidad

  number log_Lr_prior
  number log_sr_prior
  number log_b_prior
  number log_beta_prior
  number log_Linf_prior
  number log_k_prior
   
  !! log_Linf_prior = log(Par_bio(1));
  !! log_k_prior = log(Par_bio(2));
  !! log_Lr_prior = log(Par_bio(3));
  !! log_sr_prior = log(Par_bio(4));
  !! log_beta_prior= log(Par_bio(5));
  !! log_b_prior = log(bprior);

 init_number L50prior
 init_number s1prior
 init_number s2prior
 init_int opt_sel2 // opcion domo en flota
 init_int opt_sel3 // opcion domo en crucero

 number log_L50prior
 number log_s1prior
 number log_s2prior

 !! log_L50prior = log(L50prior);
 !! log_s1prior = log(s1prior);
 !! log_s2prior = log(s2prior);

 init_int    nbloques1
 init_vector ybloques1(1,nbloques1)
 init_int    nbloques2
 init_vector ybloques2(1,nbloques2)
 init_int    nqbloques
 init_vector yqbloques(1,nqbloques)
 init_int    nqbloquesc
 init_vector yqbloquesc(1,nqbloquesc)
 init_int    nqbloquesmph
 init_vector yqbloquesmph(1,nqbloquesmph)

 init_number      prior_qf 
 init_number      prior_qc 
 init_number      prior_qmph 

 number      log_prior_qf 
 number      log_prior_qc 
 number      log_prior_qmph 

 !!log_prior_qf=log(prior_qf);
 !!log_prior_qc=log(prior_qc); 
 !!log_prior_qmph =log(prior_qmph);

 init_vector cv_q(1,3) 

 init_int    opt_qf
 init_int    opt_qc
 init_int    opt_qmph
 
 init_int    opt1_fase // selectividad flota
 init_int    opt2_fase // selectividad cruceros

 init_int    opt_Lr
 init_int    opt_sr
 init_int    opt_beta
 init_int    opt_vb1
 init_int    opt_vb2
 
 init_int    opt_Ro       //Opcion para estimar o dejar fijo log_Ro
 init_number log_priorRo  // Valor fijo para log_Ro, perfiles de verosimilitud 
 
 init_int    opt_F
 init_int    opt_devRt
 init_int    opt_devNo//Condicion inicial (Si no estima significa poblaciónen equilibrio)
 init_int    opt_bpow //hiperestabilidad

 init_int    npbr
 init_vector pbr(1,npbr)
 init_int ntime_sim

 init_number  festim
 init_number  oprec
 // 14. Para la fase de proyeccion con relacion stock-recluta
 init_int nsim //numero de años de simulacion
 init_int nhcr //numero de harvest control rule
 //init_number h
 init_int SrType
 //init rho
 init_number rho_r
 
//========================================================================

INITIALIZATION_SECTION

//========================================================================
  log_Rmed       8.8
  log_Lr         log_Lr_prior
  log_sr         log_sr_prior
  log_L50        log_L50prior 
  log_L50c       log_L50prior 
  log_sigma1     log_s1prior 
  log_sigma2     log_s2prior
  log_sigma1c    log_s1prior 
  log_sigma2c    log_s2prior
  log_b          log_b_prior 
  log_beta       log_beta_prior 
  log_k          log_k_prior
  log_Linf       log_Linf_prior
  log_qflo       log_prior_qf
  log_qcru       log_prior_qc
  log_qmph       log_prior_qmph
  //log_fpbr  -1.204
//========================================================================

PARAMETER_SECTION

//========================================================================
// selectividad paramétrica a la talla común
// init_bounded_vector log_L50f(1,nbloques1,-5,8,opt1_fase)  
 init_vector log_L50(1,nbloques1,opt1_fase)  
 init_vector log_sigma1(1,nbloques1,opt1_fase)
 init_vector log_sigma2(1,nbloques1,opt_sel2)

 init_vector log_L50c(1,nbloques2,opt2_fase)  
 init_vector log_sigma1c(1,nbloques2,opt2_fase)
 init_vector log_sigma2c(1,nbloques2,opt_sel3)

// parametros reclutamientos y mortalidades)
 init_number log_Rmed(opt_Ro)
 init_bounded_dev_vector log_desv_Rt(1,ntime,-10,10,opt_devRt)
 init_bounded_vector log_desv_No(1,nedades,-10,10,opt_devNo)
 init_bounded_vector log_F(1,ntime,-20,0.8,opt_F) // log  mortalidad por pesca por flota
 //Agrega Cubillos
 //init_bounded_number log_fpbr(-20,0.8,1) //para calcular f60
// capturabilidades
 init_vector log_qflo(1,nqbloques,opt_qf)
 init_vector log_qcru(1,nqbloquesc,opt_qc)
 init_vector log_qmph(1,nqbloquesmph,opt_qmph)
 //vector log_qmph(1,nqbloquesmph) //Cubillos

 init_number log_b(opt_bpow)

// Crecimiento
 init_number log_Lr(opt_Lr)
 init_number log_sr(opt_sr)
 init_number log_beta(opt_beta)
 init_number log_Linf(opt_vb1)
 init_number log_k(opt_vb2)

//---------------------------------------------------------------------------------
//Defino las variables de estado 
 vector BMflo(1,ntime)
 vector BMcru(1,ntime)
 vector Brec(1,ntime)
 vector BMmph(1,ntime)
 vector pred_CPUE(1,ntime);
 vector pred_Bcru(1,ntime);
 vector pred_Desemb(1,ntime);
 vector likeval(1,9);
 vector Neq(1,ntallas);

 vector Rpred(1,ntime);
 vector Unos_edad(1,nedades);
 vector Unos_year(1,ntime);
 vector Unos_tallas(1,ntallas);
 vector delta(1,ntallas)
 vector Lesp(1,ntallas)
 vector sigmaL(1,ntallas)
 vector pre(1,ntallas)

 vector mu_edad(1,nedades)
 vector sigma_edad(1,nedades)
 vector BDo(1,ntime);
 vector No(1,ntallas)
 vector prior(1,7)
 vector yrs(1,ntime)
 vector Desemb(1,ntime);
 vector CPUE(1,ntime);
 vector Bcru(1,ntime);
 vector dt_C(1,ntime);
 vector Frms_bloque(1,ntime);
 vector mph(1,ntime);
 vector Lmed_obs(1,ntime)
 vector Lmed_pred(1,ntime)
 vector Lmed_obsc(1,ntime)
 vector Lmed_predc(1,ntime)
 vector edades(1,nedades)
 sdreport_vector Reclutas(1,ntime)
 vector nm(1,ntime)
 vector nmc(1,ntime)
 vector penalty(1,7)


 matrix cv_index(1,4,1,ntime)

 matrix S1(1,nbloques1,1,ntallas)
 matrix S2(1,nbloques2,1,ntallas)
 matrix Sel(1,ntime,1,ntallas)
 matrix Selc(1,ntime,1,ntallas)
 matrix F(1,ntime,1,ntallas)
 matrix Z(1,ntime,1,ntallas)
 matrix S(1,ntime,1,ntallas)
 matrix N(1,ntime,1,ntallas)

 matrix NM(1,ntime,1,ntallas)
 matrix NMD(1,ntime,1,ntallas)
 matrix NDv(1,ntime,1,ntallas)
 matrix Nrec(1,ntime,1,ntallas)
 matrix NVflo(1,ntime,1,ntallas)
 matrix NVcru(1,ntime,1,ntallas)
 matrix NVmph(1,ntime,1,ntallas)
 matrix TEMP(1,ntime,1,ntallas)
 matrix pred_Ctot(1,ntime,1,ntallas)
 matrix pobs(1,ntime,1,ntallas)
 matrix ppred(1,ntime,1,ntallas)
 matrix pobsc(1,ntime,1,ntallas)
 matrix ppredc(1,ntime,1,ntallas)
 matrix T(1,ntallas,1,ntallas)
 //agrega Id Cubillos
 matrix Id(1,ntallas,1,ntallas); //Matriz identidad - Cubillos
 matrix Nv(1,ntime,1,nedades)
 matrix NMDv(1,ntime,1,ntallas)
 number suma1
 number suma2
 number suma3
 number suma4
 number suma5
  
 number So
 number alfa
 number beta

 number Linf
 number k
 number Linfh
 number M
 number Lr
 number sr
 number Lm
 number Rm
 number h

 number BDp
 number Npplus
 number Bp_anch 

 number nm1;
 number cuenta1;
 number alfa_sr;
 number beta_sr;
 number pF

 vector Npact(1,ntallas)
 vector Np(1,ntallas)
 vector Zpbr(1,ntallas)
 vector Fpbr(1,ntallas)
 vector Sp(1,ntallas)

 matrix Bp(1,npbr,1,ntime_sim)
 
 //cba actual
 vector NV(1,ntallas)
 vector CTP(1,ntallas)
 vector Yact(1,npbr)
 //cba proyectada
 vector NVp(1,ntallas) 
 vector CTPp(1,ntallas)
 matrix Yp(1,npbr,1,ntime_sim)
 
 matrix Rpp(1,npbr,1,ntime_sim)
 
 objective_function_value f
  
 sdreport_vector BD(1,ntime) // 
 sdreport_vector BT(1,ntime) // 
 vector RPRlp(1,ntime) // 
 vector RPR(1,ntime)
 vector Frpr(1,ntime)
 vector pred_mph(1,ntime);
 sdreport_number SSBo
 
 //Agrega Cubillos NPR
 number SPR0
 number SPRf
 vector Nspr(1,ntallas)
 vector Nsprf(1,ntallas)
 vector zpbr(1,ntallas)
 number fpbr
 number penpbr
 number ratio_spr
 
 //================ Modificacion LA Cubillos MJ Cuevas =====
 //Relacion stock-recluta
 number RO
 number SO
 //number alfa_sr
 //number beta_sr
 //hcr
 number p1
 number p2
 number r1
 number r2
 //Reporte para R
 sdreport_vector fyr(1,ntime)
 //vector wi(1,nedades)
 vector wm(1,ntallas)
 vector sel(1,ntallas)
	 
   //Reporte para proyeccion a largo plazo
 //number Fref
 //number ssb_fut
 //number SObj
 //number FObj
 //vectores por edad
 //vector N_fut(1,ntallas)
 //vector Fmort(1,ntallas)
 //vector Z_fut(1,ntallas)
 //vector S_fut(1,ntallas)
 //vector CBAn(1,ntallas)
 //vector NVpel_fut(1,nedades)
 //vector NVrec_fut(1,nedades)
 //vectores en tiempo
 //vector Bpel_fut(1,nsim)
 //vector Brec_fut(1,nsim)
 //matrices por hcr y tiempo
 //matrix BT_fut(1,nhcr,1,nsim)
 //matrix BD_fut(1,nhcr,1,nsim)
 //matrix R_fut(1,nhcr,1,nsim)
 //matrix F_fut(1,nhcr,1,nsim)
 //matrix Depl(1,nhcr,1,nsim)
 //matrix Y_fut(1,nhcr,1,nsim)

	 
 // VARIABLES PARA EL MODELO OPERATIVO
 // Proyeccion futura
	  vector yrs_fut(ntime+1,ntime+nsim);
      //number SO;
     number SObj;
	  //number RO;
      number ssb_fut;
	  number BTp;
	  number Kobs_tot_catch;
	  matrix C_fut(1,npbr,ntime+1,ntime+nsim);
	  //matrix future_biomass(1,npbr,ntime+1,ntime+nsim);//matrix o vector??
	  vector expl_biom(ntime+1,ntime+nsim);	  
	  //matrix future_ssbiom(1,npbr,ntime+1,ntime+nsim);
	  vector N_fut(1,ntallas);
	  matrix F_fut(ntime+1,ntime+nsim,1,ntallas);
	  //matrix Zpbr(nanos+1,nanos+nproy,1,nedades);
	  //vector Zpbr(2,nedades);
	  //vector Sp(1,nedades);
	  matrix check_convergence(1,npbr,ntime+1,ntime+nsim);
	  matrix Fyr_future(1,npbr,ntime+1,ntime+nsim);
	  vector NVflo_fut(1,ntallas);
	  vector NVcru_fut(1,ntallas);
	  vector NMDp(1,ntallas);
		  
		 
	  matrix F_tmp(ntime+1,ntime+nsim,1,ntallas);
	  matrix M_tmp(ntime+1,ntime+nsim,1,ntallas);
      vector BMflo_fut(ntime+1,ntime+nsim);
	  vector BMcru_fut(ntime+1,ntime+nsim);
	  vector BMpel_fut(ntime+1,ntime+nsim);
	  vector BMmph_fut(ntime+1,ntime+nsim);
	  //vector NMD_fut(1,ntallas);

	 //Datos simulados
	  matrix fut_Ctot(ntime+1,ntime+nsim,1,ntallas);
	  matrix fut_pcru(ntime+1,ntime+nsim,1,ntallas);
	   //vector wm(1,ntallas);
	  //vector Wm_fut(1,ntallas);
	  //vector Win_fut(1,ntallas);
	  vector error_edad_fut(1,ntallas);	
   
	  //vector Win_fut(1,nanos);	
	  number R_fut;
	  vector Sel_f(1,ntallas);
	  vector Sel_c(1,ntallas);
	  //matrix Scru_pela_fut(nanos+1,nanos+nproy,1,nedades);
	  //vector Sel_f_fut(1,nedades);
	  //vector sel(1,nedades);
	  //vector Np(1,nanos);
	  vector rec_epsilon_fut(ntime+1,ntime+nsim);
	  //vector CPUEI_epsilon_fut(ntime+1,ntime+nsim);
	  vector CPUE_epsilon_fut(ntime+1,ntime+nsim);
	  vector Bcru_epsilon_fut(ntime+1,ntime+nsim);
	  vector mph_epsilon_fut(ntime+1,ntime+nsim);
	  
	  //vector sim_CPUEI(ntime+1,ntime+nsim);
	  vector sim_CPUE(ntime+1,ntime+nsim);
	  vector sim_mph(ntime+1,ntime+nsim);
	  vector sim_Bcru(ntime+1,ntime+nsim);
	  
	  matrix rpr_fut(1,npbr,ntime+1,ntime+nsim);
	  matrix des_Btot(1,npbr,ntime+1,ntime+nsim);
	  matrix des_SSBt(1,npbr,ntime+1,ntime+nsim);
	  matrix des_Rt(1,npbr,ntime+1,ntime+nsim);
	  matrix des_F(1,npbr,ntime+1,ntime+nsim);
	  matrix des_Yt(1,npbr,ntime+1,ntime+nsim);
	  matrix des_rpr(1,npbr,ntime+1,ntime+nsim);
	  matrix keep_Btot(1,npbr,ntime+1,ntime+nsim);
	  matrix keep_SSBt(1,npbr,ntime+1,ntime+nsim);
	  matrix keep_Rt(1,npbr,ntime+1,ntime+nsim);
	  matrix keep_F(1,npbr,ntime+1,ntime+nsim);
	  matrix keep_RBB(1,npbr,ntime+1,ntime+nsim);
	  matrix des_RBB(1,npbr,ntime+1,ntime+nsim);
	  
	  vector rec_epsilon(ntime,ntime+nsim);  
 
 //======= Fin MODIFICACIONES =======

  

  
//========================================================================

PRELIMINARY_CALCS_SECTION

//========================================================================
 yrs=column(mdatos,1);
 Desemb=column(mdatos,2);
 CPUE=column(mdatos,4);
 Bcru=column(mdatos,6);
 nm=column(mdatos,8);
 nmc=column(mdatos,9);
 mph=column(mdatos,10);
 dt_C=column(mdatos,12);
 Frms_bloque=column(mdatos,13);
 //cout << "dt_C:" << dt_C << endl;exit(1);

 edades.fill_seqadd(minedad,1);

 cv_index(1)=column(mdatos,3);
 cv_index(2)=column(mdatos,5);
 cv_index(3)=column(mdatos,7);
 cv_index(4)=column(mdatos,11);

 Linf=Par_bio(1);
 k=Par_bio(2);
 M=Par_bio(6);
 h=Par_bio(7);

 Unos_tallas=1;// lo uso en operaciones matriciales con tallas
 Unos_year=1;// lo uso en operaciones matriciales con el año
 //Modificado por Cubillos
 int j;
 for(j=1;j<=ntallas;j++){
  Id(j,j)=1.0;
 }
//========================================================================

RUNTIME_SECTION

//========================================================================
  convergence_criteria 1.e-1,1.e-01,1.e-03,1e-3,1e-5
  maximum_function_evaluations 100,100,200,3000,3500

//========================================================================

PROCEDURE_SECTION

//========================================================================
// se listan las funciones que contienen los calculos
 Eval_Trans_talla_talla();
 Eval_selectividad();
 Eval_mortalidades();
 Eval_abundancia();
 Eval_biomasas();
 Eval_capturas_predichas();
 Eval_indices();
 //Eval_pbr();
 Eval_logverosim();
 Eval_funcion_objetivo();

 //if(last_phase()){Eval_CTP();}
 if(mceval_phase()){
 	Eval_Mop();
 }

FINAL_SECTION
 	 //Eval_CTP();
     Write_R(); //Agregada por Cubillos
	 //Eval_Mop();
 
  time(&finish);
  elapsed_time=difftime(finish,start);
  hour=long(elapsed_time)/3600;
  minute=long(elapsed_time)%3600/60;
  second=(long(elapsed_time)%3600)%60;
  cout<<endl<<endl<<"*********************************************"<<endl;
  cout<<"--Start time:  "<<ctime(&start)<<endl;
  cout<<"--Finish time: "<<ctime(&finish)<<endl;
  cout<<"--Runtime: ";
  cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
  cout<<"*********************************************"<<endl;

FUNCTION Eval_Mop

    dvector epsilon_rec_fut(ntime+1,ntime+nsim);
    //dvector epsilon_CPUEI_fut(ntime+1,ntime+nsim); 
    dvector epsilon_CPUE_fut(ntime+1,ntime+nsim);
    dvector epsilon_Bcru_fut(ntime+1,ntime+nsim);
    dvector epsilon_mph_fut(ntime+1,ntime+nsim);
   	int upk;//actualizar los años contador para cada ciclo
   	int k,l,j,i,a;//contador de los ciclos for
	dvariable numyear=2002;
	dvector yrs(ntime+1,ntime+nsim);
	simname="SAMsau.dat";//datos del estimador
   	dvector CatchNow(1,npbr);//capturas del estimador
   	dvector eval_now(1,5);//son 5 los indicadores del estimador:BT,BD,R,F,B/Bo
   	dvariable grad_tmp;
	dvariable sigr;
       	   	
	epsilon_rec_fut.fill_randn(rng);//error de proceso
    //epsilon_CPUEI_fut.fill_randn(rng);
    epsilon_CPUE_fut.fill_randn(rng);
    epsilon_Bcru_fut.fill_randn(rng);
    epsilon_mph_fut.fill_randn(rng);
    int rv; 
 
      for(int l=1;l<=npbr;l++){
 		//int rv=system("cp MATT2009_real.dat MATT2009.dat");
 		//rv=system("./MATT2009 -nox -nohess");
		rv= system("./paso1.sh");
		N_fut = N(ntime);
 		Sp = S(ntime);
 		BDp = BD(ntime);
 		Zpbr = Z(ntime);
 	    Sel_f = Sel(ntime);
		Sel_c = Selc(ntime);		
 		//cout << "sel"<< sel<< endl;exit(1);
		sigr = sigmaR;
    	SObj = BD(ntime)/(0.55*SO);
    	ssb_fut  = BD(ntime);
		rec_epsilon(ntime)=log_desv_Rt(ntime); //copia la desv del ultimo año
    	
		//cout << "ssb_fut" << ssb_fut << endl;
		for(int i=ntime+1;i<=ntime+nsim;i++){ //Ciclo periodo de proyeccion
		    rec_epsilon(i) = rho_r*rec_epsilon(i-1)+epsilon_rec_fut(i)*sigr*pow(1-rho_r,0.5);				
		    rec_epsilon_fut(i)=rec_epsilon(i);
		 //CPUEI_epsilon_fut(i)=epsilon_CPUEI_fut(i)*0.2;
		 CPUE_epsilon_fut(i)=epsilon_CPUE_fut(i)*0.2;
	     Bcru_epsilon_fut(i)=epsilon_Bcru_fut(i)*0.15;
	     mph_epsilon_fut(i)=epsilon_mph_fut(i)*0.3;
   
		ifstream CTP_tmp("cba.dat");
        CTP_tmp >> CatchNow;
        CTP_tmp.close();
		C_fut(l,i)=CatchNow(l);
	  
        //cout << "C_fut " << C_fut << endl;exit(1);
		//N_fut(2,ntallas)= ++elem_prod(N_fut(1,ntallas-1),Sp(1,ntallas-1));
		
		if(SrType==1){
			R_fut=(BDp/(alfa_sr + beta_sr * BDp))*exp(rec_epsilon_fut(i));
		}
		if(SrType==2){
			R_fut=(alfa_sr*BDp*mfexp(-beta_sr* BDp))*exp(rec_epsilon_fut(i));
		}
		
		N_fut=elem_prod(N_fut,Sp)*T+pre*R_fut;//*N_fut(i,j);  
   
        NVflo_fut= elem_prod(N_fut,Sel_f);
      	//cout << "HASTA AQUI= " << NVflo_fut << endl;exit(1);	   
		expl_biom(i)=0;
		for(j=1;j<=ntallas;j++){
         //expl_biom(i)+=Sel_f_fut(i,j)*N_fut(i,j)*Wm_fut(i,j);
         expl_biom(i)+=Wmed(j)*NVflo_fut(j);
         }
 		if(C_fut(l,i)!=0){
  			//M_tmp = M;
  			dvariable ffpen = 0.0;
  			dvariable SK = posfun((expl_biom(i) - C_fut(l,i))/expl_biom(i),0.05,ffpen);
  			Kobs_tot_catch = expl_biom(i)-SK*expl_biom(i);
  			do_Newton_Raphson_for_mortality(i);
  			F_fut(i)=F_tmp(i);} 
 		else{
 			F_fut(i) = 0;
 		}
		//cout << "F_fut i " << F_fut(i) << endl; exit(1);	  
        Zpbr=F_fut(i)+M;
        Sp=exp(-1.*Zpbr);
 		//Reclacula los indicadores después que la mortalidad fue calculada
 		//Biomasa vulnerable
 		NVflo_fut = elem_prod(elem_prod(N_fut,mfexp(-dt(2)*Zpbr)),Sel_f);
 		BMflo_fut(i) = Wmed*NVflo_fut; //Biomasa vulnerable 				
 		
		
		//Crucero
 		NVcru_fut = elem_prod(elem_prod(N_fut,mfexp(-dt(3)*(Zpbr))),Sel_c);
 		BMcru_fut(i) = Wmed*NVcru_fut;
 		NMDp=elem_prod(elem_prod(N_fut,mfexp(-dt(1)*(Zpbr))),msex); //Numero desovante
		
		BMmph_fut(i)=sum(elem_prod(Wmed,NMDp));
		//cout << "HASTA AQUI= " << BMmph_fut<< endl;exit(1);	
  	 
		
 	    //NMDp=elem_prod(NMDp,outer_prod(Unos_year,msex)); //CONSULTAR
 		BDp=sum(elem_prod(Wmed,NMDp)); //biomasa desovante mph
 		BTp=sum(elem_prod(N_fut,Wmed));
	
		//CPUE
		sim_CPUE(i) = exp(log_qflo(nqbloques))*pow(BMflo_fut(i),exp(log_b))*exp(CPUE_epsilon_fut(i));
		//Biomasa crucero
		sim_Bcru(i) = exp(log_qcru(nqbloquesc))*BMcru_fut(i)*exp(Bcru_epsilon_fut(i)); 
        //MPH		
        sim_mph(i) = exp(log_qmph(nqbloquesmph))*BMmph_fut(i); //*exp(mph_epsilon_fut(i)); //??	
	  
	    //Estructura de tallas - pesqueria y cruceros
		for(int j=1;j<=ntallas;j++){//captura por edad. cambiar caja de º
	    fut_Ctot(i,j) = (F_fut(i,j)/Zpbr(j))*N_fut(j)*(1.-Sel_f(j));
		fut_pcru(i,j)= NVcru_fut(j);
		}   	
//=========== ESCRITURA DEL ARCHIVO *.dat ===============	
		upk=i;
		yrs_fut(i) = numyear+i-1;
		ofstream simdata(simname);
		simdata << "# SIMULACION DE DATOS" << endl;
		simdata << "# Numero de anos actual" << endl;
		simdata << upk << endl;
		simdata << "# Numero de edades " << endl;
		simdata << nedades << endl;
		simdata << "# Numero de clases de talla " << endl;
		simdata << ntallas << endl;
		simdata << "# MATRIZ DE DATOS HISTORICOS " << endl;
		simdata << "#Yr Desemb cv_Desemb CPUEI cv_CPUEI CPUEA cv_CPUEA Bcruc cv_Bcruc Mpdh cv_Mpdh nm_flota nm_cruc" <<endl;
	    simdata << mdatos << endl;
		simdata << "# Agrega datos simulados a la matriz de datos historica" << endl;
		for(k=ntime+1;k<=upk;k++){
			simdata << yrs_fut(k) << "  " << C_fut(l,k) << "  " <<0.1 <<" "<<  sim_CPUE(k) << "  " << 0.26 <<" "<< sim_Bcru(k) << "  " << 0.22 << "  " << 25 << "  " << 15 << "  " << 0 << "  "<< 0 << "  " << 0.33 << " " << 0.31 << endl;
		}
		simdata << "#Estructura de tallas flota historica" << endl;
		simdata << "#Tallas" << endl;
		simdata << Tallas << endl;
		simdata << Ctot << endl;
		simdata << "# flota simulados " << endl;
		for(k=ntime+1;k<=upk;k++){
			simdata << fut_Ctot(k) << endl;
		}
		simdata << "#Crucero historico " << endl;
		simdata <<  Ncru << endl;
		simdata << "#crucero simulados " << endl;
		for(k=ntime+1;k<=upk;k++){
			simdata  <<fut_pcru(k) << endl;
		}		
		simdata.close();
		//rv=system("./MATT2009  -nox -nohess");
		rv=system("./paso2.sh");
		//Lee y retiene los indicadores del estimador en la evaluación actual
		ifstream eval_estatus("estatus.dat");
	
		eval_estatus >> eval_now;
		eval_estatus.close();		
		keep_Btot(l,i)=eval_now(1); //biomasa total
		keep_SSBt(l,i)=eval_now(2); //biomasa desovante
		keep_Rt(l,i) =eval_now(3); //reclutas
		keep_F(l,i) =eval_now(4); //Mortalidad por pesca
		keep_RBB(l,i)=eval_now(5); //B/Bo
		//  
		//Guarda los indicadores del simulador que debío estimar el estimador (valga la redundancia)
				
		des_Btot(l,i)=BTp; //Biomasa total
		des_SSBt(l,i)=BDp; //Biomasa desovante
		des_Rt(l,i) = R_fut; //Reclutamiento 
		des_F(l,i) = max(F_fut(i)); //Mortalidad por pesca 
		des_RBB(l,i) = BDp/SO;
	
		ifstream lee_grad("grad_final.dat");
		lee_grad >> grad_tmp;
		lee_grad.close();
		if(grad_tmp>0.001){check_convergence(l,i)=0;}else{check_convergence(l,i)=1;}
		} //cierra ciclo de proyeccion
	    } //cierra ciclo del numero de HCR		
 
	 //======== REPORTA SIMLACIONES Y ESTIMACIONES
	 //Capturas efectivas		       
	    rep1 << C_fut << endl;
		//Guarda variables de estado del simulador y estimador
		//Simulador
		rep2 << des_Btot << endl;
		rep3 << des_SSBt << endl;
		rep4 << des_Rt << endl;
		rep5 << des_F << endl;
		rep15 << des_RBB <<endl;
		//Estimador
		rep6 << keep_Btot << endl;
		rep7 << keep_SSBt << endl;
		rep8 << keep_Rt << endl;
		rep9 << keep_F << endl;
		rep14 << keep_RBB <<endl; 
	 //Variables en la fase histórica
	    rep10 << fyr << endl; //mortalidad pesca         
	    rep11 << BD << endl; //desovante
	    rep12 << BT << endl; //biomasa total
	    rep13 << Reclutas << endl; //Reclutas
		rep16 << BD/SO <<endl; //Historico
	   //reporta el gradiente CONSULTARRRR
	   rep_grad << check_convergence << endl;
	   

	  
FUNCTION void do_Newton_Raphson_for_mortality(int i)
	
				dvariable Fold = Kobs_tot_catch/expl_biom(i);
				dvariable Fnew;
				for (int ii=1;ii<=5;ii++)
				{
				dvariable ZZ = Fold + M;
				dvariable XX = exp(-1.*ZZ);
				dvariable AA = Fold * (1. - XX);
				dvariable BB = ZZ;
				dvariable CC = 1. + (Fold - 1) * XX;
				dvariable dd = 1.;
				dvariable FX = AA / BB - Kobs_tot_catch/expl_biom(i);
				dvariable FPX = (BB * CC - AA * dd) / (BB * BB);
				Fnew = Fold - FX / FPX;
				Fold = Fnew;
				}
				F_tmp(i)= Fnew*Sel_f; 
				
	



 //#########################
 /*
 FUNCTION Eval_pbr
    
    fpbr   = exp(log_fpbr);
 	Nspr   = pre*1*inv(Id - T*(exp(-1.*M)*Id));
    SPR0   = sum(elem_prod(Nspr*exp(-dt(1)*M),elem_prod(Wmed,msex)));
 	zpbr   = M + Sel(ntime)*fpbr;
    Nsprf=pre*1; //*inv(Id - T*(exp(-1*zpbr)*Id));
 	for(int j=1;j<=10;j++){
		Nsprf=elem_prod(Nsprf,exp(-1.*zpbr))*T+pre*1;
	}
 	//Nsprf  = (elem_prod(Nspr,exp(-zpbr)))*T+pre*1;
     Nsprf  = elem_prod(elem_prod(Nsprf,mfexp(-dt(1)*zpbr)),msex);
     SPRf   = sum(elem_prod(Wmed,Nsprf));
     ratio_spr = SPRf/SPR0;    

 */
//========================================================================

FUNCTION Eval_Trans_talla_talla

//========================================================================
  Linf=exp(log_Linf);
  k=exp(log_k);
  beta=exp(log_beta);

//  if(active(log_k)){k=mfexp(log_k);}
  if(active(log_beta)){beta=mfexp(log_beta);}

 int i, j;
 
// matriz de transicion modelo normal

  delta=(Linf-Tallas)*(1-mfexp(-k));// incremento en tallas
  Lesp=Tallas+delta; // talla esperada luego del crecimiento
  sigmaL=delta*beta;  

  for (i=1;i<=ntallas;i++){
    for (j=1;j<=ntallas;j++){
      if(i==j){
         T(i,j)=1.0;}}
   }


  for (i=1;i<=ntallas;i++){

    for (j=1;j<=ntallas;j++){
     if(sigmaL(i)>0){
     T(i,j)=mfexp(-0.5*square((Lesp(i)-Tallas(j))/sigmaL(i)));}}
   }


  for (j=1;j<=ntallas;j++){
  T(j)/=sum(T(j));
  } 

//========================================================================

FUNCTION Eval_selectividad

//========================================================================
 int i,j;

 // FLOTA...................

 for (j=1;j<=nbloques1;j++){

 S1(j)=exp(-0.5*square(Tallas-exp(log_L50(j)))/square(exp(log_sigma1(j))));


    for (i=1;i<=ntallas;i++){

      if(Tallas(i)>=exp(log_L50(j))){
      S1(j,i)= exp(-0.5*square(Tallas(i)-exp(log_L50(j)))/square(exp(log_sigma2(j))));
      }

 }}

   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques1;j++){
              if (yrs(i)>=ybloques1(j)){
                Sel(i)=S1(j);}
       }
   }

 // CRUCERO...................

 for (j=1;j<=nbloques2;j++){

 S2(j)=exp(-0.5*square(Tallas-exp(log_L50c(j)))/square(exp(log_sigma1c(j))));

    for (i=1;i<=ntallas;i++){

      if(Tallas(i)>=exp(log_L50c(j))){
      S2(j,i)= exp(-0.5*square(Tallas(i)-exp(log_L50c(j)))/square(exp(log_sigma2c(j))));
      }

 }}


   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques2;j++){
              if (yrs(i)>=ybloques2(j)){
                Selc(i)=S2(j);}
       }
   }

//========================================================================

FUNCTION Eval_mortalidades

//========================================================================

 F=elem_prod(Sel,outer_prod(mfexp(log_F),Unos_tallas));
 Z=F+M;
 S=mfexp(-1.0*Z);
 fyr = mfexp(log_F); //Modifica L.A. Cubillos
 
//========================================================================

FUNCTION Eval_abundancia

//========================================================================
 int i, j;
 
 if(opt_Ro<0)
  {
  log_Rmed=log_priorRo;
  }

  Lr=Par_bio(3);
  sr=Par_bio(4);

  if (active(log_Lr)){Lr=mfexp(log_Lr);}
  if (active(log_sr)){sr=mfexp(log_sr);}


// genero la composicion de tallas del reclutamiento
  pre=exp(-0.5*square((Tallas-Lr)/sr));
  pre/=sum(pre);

// genero una estructura inicial en torno a Z del primer año;
  Reclutas=mfexp(log_Rmed+log_desv_Rt);

// genero la poblacion en equilibrio virginal de LP;

  No=pre*exp(log_Rmed);
  for (int j=1;j<=5*nedades;j++){
  No=(No*exp(-1.*M))*T+pre*exp(log_Rmed); }
  
  SSBo    = sum(elem_prod(No*mfexp(-dt(1)*M),elem_prod(Wmed,msex)));
  //alfa_sr = 4*h*exp(log_Rmed+0.5*square(sigmaR))/(5*h-1);//
  //beta_sr = (1-h)*SSBo/(5*h-1);// Reclutamiento
  //########################################################   	
  //Modifica LA Cubillos - MJ Cuevas
  //STOCK-RECLUTA
  RO = mfexp(log_Rmed);//mfexp(log_Rmed + 0.5*square(sigmaR));
  SO = SSBo;
  if(SrType==1){
	  alfa_sr = (SO/RO)*((1-h)/(4*h));
	  beta_sr = (5*h-1)/(4*h*RO);
	  }
	  else
	  {
		  alfa_sr = (RO/SO)*exp(log(5*h)/0.8);
		  beta_sr = log(5*h)/(0.8*SO);
	  }	
  //######################################################	

// -----------------primer año
  Reclutas(1) = mfexp(log_Rmed+log_desv_Rt(1));
  Rpred(1)    = Reclutas(1);
// genero una estructura inicial en torno a Z del primer año;
  Neq=pre*Reclutas(1);
  for (j=1;j<=nedades;j++){
  Neq    = elem_prod(Neq,exp(-1.*Z(1)))*T+pre*exp(log_Rmed+log_desv_No(j));}
  N(1)   = Neq;
  NMD(1) = elem_prod(elem_prod(N(1),mfexp(-dt(1)*Z(1))),msex);
  BD(1)  = sum(elem_prod(Wmed,NMD(1)));

// --------------------dinamica anual
  for (i=2;i<=ntime;i++){
  Reclutas(i) = mfexp(log_Rmed+log_desv_Rt(i));
  Rpred(i)    = Reclutas(i);
  
  if(i>minedad){
  //Rpred(i)    = (alfa_sr*BD(i-minedad)/(beta_sr + BD(i-minedad)));
  if(SrType==1){Rpred(i)=BD(i-minedad)/(alfa_sr + beta_sr*BD(i-minedad));}
  if(SrType==2){Rpred(i)=alfa_sr*BD(i-minedad)*exp(-beta_sr*BD(i-minedad));}
  
  Reclutas(i) = Rpred(i)*mfexp(log_desv_Rt(i)); }

  N(i)        = (elem_prod(N(i-1),S(i-1)))*T+pre*Reclutas(i);
  NMD(i)      = elem_prod(elem_prod(N(i),mfexp(-dt(1)*Z(i))),msex);
  BD(i)       = sum(elem_prod(Wmed,NMD(i)));
  } 

//========================================================================

FUNCTION Eval_biomasas

//========================================================================
 NMD=elem_prod(N,mfexp(-dt(1)*Z));
 NMD=elem_prod(NMD,outer_prod(Unos_year,msex));
 NVflo=elem_prod(elem_prod(N,mfexp(-dt(2)*(Z))),Sel);
 NVcru=elem_prod(elem_prod(N,mfexp(-dt(3)*(Z))),Selc);
 
 for(int i=1;i<=ntime;i++){
  NVcru(i)=elem_prod(elem_prod(N(i),mfexp(-dt_C(i)*(Z(i)))),Selc(i));
 }

// cout << "dt_C:" << dt_C << endl;exit(1);
// TEMP=outer_prod(ntallas,dt_C);
// NVcru=elem_prod(elem_prod(N,mfexp(outer_prod(dt_C,ntallas)*(-Z))),Selc);
// NVmph=elem_prod(elem_prod(N,mfexp(-dt(4)*(Z))),Selc);

// vectores de biomasas derivadas
 BD    = Wmed*trans(NMD);
 BMflo = Wmed*trans(NVflo);
 BMcru = Wmed*trans(NVcru);
 BMmph = Wmed*trans(NVmph);
 BT    = Wmed*trans(N);
 BDo   = sum(elem_prod(No,Wmed));
 RPRlp = BD/SSBo;
 RPR   = BD/(SSBo*0.55);
 Frpr  = elem_div(mfexp(log_F),Frms_bloque);
//========================================================================

FUNCTION Eval_capturas_predichas

//========================================================================

// matrices de capturas predichas por edad y año
 pred_Ctot=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));

// vectores de desembarques predichos por año
 pred_Desemb=Wmed*trans(pred_Ctot);

// matrices de proporcion de capturas por talla y año
 pobs=elem_div(Ctot,outer_prod(rowsum(Ctot+1e-10),Unos_tallas));
 ppred=elem_div(pred_Ctot,outer_prod(rowsum(pred_Ctot+1e-10),Unos_tallas));

 pobsc=elem_div(Ncru,outer_prod(rowsum(Ncru+1e-10),Unos_tallas));
 ppredc=elem_div(NVcru,outer_prod(rowsum(NVcru+1e-10),Unos_tallas));

 Lmed_pred=Tallas*trans(ppred);
 Lmed_obs=Tallas*trans(pobs);

 Lmed_predc=Tallas*trans(ppredc);
 Lmed_obsc=Tallas*trans(pobsc);

//========================================================================

FUNCTION Eval_indices

//========================================================================

   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques;j++){
              if (yrs(i)>=yqbloques(j)){
                 pred_CPUE(i)=exp(log_qflo(j))*pow(BMflo(i),exp(log_b));}
       }
   }


   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloquesc;j++){
              if (yrs(i)>=yqbloquesc(j)){
                 pred_Bcru(i)=exp(log_qcru(j))*BMcru(i);}
       }
   }


  for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloquesmph;j++){
              if (yrs(i)>=yqbloquesmph(j)){
                 pred_mph(i)=exp(log_qmph(j))*BD(i);}
       }
   }


//========================================================================

FUNCTION Eval_logverosim

//========================================================================
// esta funcion evalua el nucleo de las -log-verosimilitudes marginales para
// series con datos 0.
 int i;

 suma1=0; suma2=0; suma3=0; penalty=0;

 for (i=1;i<=ntime;i++)
 {
  if (CPUE(i)>0){
    suma1+=square(log(CPUE(i)/pred_CPUE(i))*1/cv_index(2,i));}

  if (Bcru(i)>0){
    suma2+=square(log(Bcru(i)/pred_Bcru(i))*1/cv_index(3,i));}

   if (mph(i)>0){
    suma3+=square(log(mph(i)/pred_mph(i))*1/cv_index(4,i));}
  }

//========================================================================

FUNCTION Eval_funcion_objetivo

//========================================================================

 suma4=0; suma5=0; penalty=0;

 likeval(1)=0.5*suma1;//CPUE
 likeval(2)=0.5*suma2;//Bcru
 likeval(3)=0.5*suma3;//mph

 likeval(4)=0.5*norm2(elem_div(log(elem_div(Desemb,pred_Desemb)),cv_index(1)));// desemb

 for (int i=1;i<=ntime;i++){
 suma4+=-nm(i)*sum(elem_prod(pobs(i),log(ppred(i)+1e-10)));
 suma5+=-nmc(i)*sum(elem_prod(pobsc(i),log(ppredc(i)+1e-10)));
 }
 
 likeval(5)=suma4;//
 likeval(6)=suma5;//

// lognormal Ninicial y Reclutas
 if(active(log_desv_Rt)){
 likeval(7)=1./(2*square(sigmaR))*norm2(log_desv_Rt);}

 if(active(log_desv_No)){
 likeval(8)=1./(2*square(sigmaR))*norm2(log_desv_No);}

 if(active(log_Lr)){
 likeval(9)=1./(2*square(cv_Par_bio(3)))*square(log_Lr-log_Lr_prior);}

  //if (active(log_F)){
  //pF=1000*norm2(log_F-mean(log_F));}

  penalty(1)=0.5/square(cv_q(1))*norm2(log_qflo-log_prior_qf);
  penalty(2)=0.5/square(cv_q(2))*norm2(log_qcru-log_prior_qc);
  //penalty(3)=0.5/square(cv_q(3))*norm2(log_qmph-log_prior_qmph);
  penalty(4)=0.5/square(0.4)*norm2(log_sigma1-log_s1prior);
  penalty(5)=0.5/square(0.4)*norm2(log_L50-log_L50prior);
  penalty(6)=0.5/square(cv_Par_bio(1))*square(log_Linf-log_Linf_prior);
  penalty(7)=0.5/square(cv_Par_bio(2))*square(log_k-log_k_prior);
  
  //Agrega cubillos
  //if(active(log_fpbr)){
  //	  penpbr=10000*square(SPRf-0.4*SPR0);
  //}

  //f=festim*(sum(likeval)+sum(penalty)+pF);
  f=festim*(sum(likeval)+sum(penalty)+penpbr);
  
  if(last_phase){
  f=festim*(sum(likeval)+sum(penalty));}
 
//========================================================================

FUNCTION  Eval_CTP

//========================================================================

  for (int i=1;i<=npbr;i++){ // ciclo de PBR
  Npact=N(ntime);
  Np=N(ntime);
  Sp=S(ntime);

 // Fpbr=F(ntime)*pbr(i);// este usa multiplicador pbr
  Fpbr=Sel(ntime)*pbr(i);// este usa la selectividad x Fpbr
  Zpbr=Fpbr+M;
  NV      = elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),Npact));
  CTP     = elem_prod(NV,Wmed);
  Yact(i) = sum(CTP);

  for (int j=1;j<=ntime_sim;j++){ // ciclo de años

  if(j<=minedad){
 // Np=(elem_prod(Np,Sp))*T+pre*(alfa_sr*BD(ntime-minedad+1)/(beta_sr+BD(ntime-minedad+1)));} // Estima CTP con R_last
 if(opt_Ro<0)
  {
  log_Rmed=log_priorRo;
  }
  
  if(oprec==1){  
  Np=(elem_prod(Np,Sp))*T+pre*(mfexp(log_Rmed));} // Estima CTP con R_med
  if(oprec==2){  
  Np=(elem_prod(Np,Sp))*T+pre*Reclutas(11);} // Estima CTP con R_alto (2012)
  if(oprec==3){  
  Np=(elem_prod(Np,Sp))*T+pre*Reclutas(17);} // Estima CTP con R_bajo (2018)
  }
 if(j>minedad){
  Np=(elem_prod(Np,Sp))*T+pre*(alfa_sr*Bp(i,j-minedad)/(beta_sr+Bp(i,j-minedad)));} //

  Bp(i,j)=sum(elem_prod(elem_prod(Np,exp(-dt(1)*Zpbr)),elem_prod(msex,Wmed)));
  
  NVp  = elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),Np));
  CTPp = elem_prod(NVp,Wmed);
  Yp(i,j)=sum(CTPp);
  Sp=exp(-1.*Zpbr);
  }}
  
  //Imprime en archivo CTP y Estatus
  ofstream salida("cba.dat");
  salida << 
  Yp(2,1) << 
  endl; 
  salida.close();

  ofstream grad_out("grad_final.dat");
  grad_out<< objective_function_value::pobjfun -> gmax <<endl;
  grad_out.close();

  ofstream out("estatus.dat"); //Biomasa total, Biom Des., Recl., F
  out << BT(ntime)   << "  " <<
  BD(ntime)          << "  " <<
  Reclutas(ntime)       << "  " <<
  fyr(ntime)      << "  " <<
  BD(ntime)/SSBo <<  endl;
  out.close();  
  
  /*
  //============= PROYECCION LA Cubillos - MJ Cuevas
  FUNCTION Proyecta_Pop

    dvector dev_rec_fut(1,nsim);
  	dev_rec_fut.fill_randn(rng);
  	FObj = pbr(2); //Ingreso externo del Fobj
  	for(int i=1;i<=nhcr;i++){
  		N_fut   = N(ntime);
  		S_fut   = S(ntime);
  		//RPR_fut = RPRequ3(nanos);
  		//Nv_fut = Nv(nanos);
  		SObj = BD(ntime)/(0.55*SO);
  		ssb_fut  = BD(ntime);
  		for(int j=1;j<=nsim;j++){ //Ciclo periodo de proyeccion
			
  			//Reclutamiento
  			if(SrType==1){
  				R_fut(i,j)=(ssb_fut/(alfa_sr + beta_sr*ssb_fut))*exp(dev_rec_fut(j)*sigmaR + 0.5*square(sigmaR));
  			}
  			if(SrType==2){
  				R_fut(i,j)=alfa_sr*ssb_fut*exp(-beta_sr*ssb_fut)*exp(dev_rec_fut(j)*sigmaR + 0.5*square(sigmaR));
  			}
			
  			N_fut=elem_prod(N_fut,S_fut)*T + pre*R_fut(i,j);
            //Reglas de control
  			if(i==1){
  				Fref = FObj;				
  			}
  			if(i==2){
  				p1=0.1;
  				p2=0.9;
  				r1=0.9;
  				r2=0.3;
  				if(SObj < 0.9){
  					Fref = r1*FObj;
  					if(SObj < 0.5){
  						Fref = FObj*(SObj*(p2-p1)/0.5+p1);
  						if(SObj < 0.25){
  							Fref = r2*FObj;
  						}
  					}
  				}
  				else{
  					Fref=FObj;
  				}
  			}
  	  	  if(i==3) //0-40
  	  	  {
  	  		  p1=0;
  	  		  p2=0.5;
  	  		  if(SObj<1){
  	  			  //Fref=FObj/(1 + exp(-(log(19)*(SObj-0.5)/(0.7-0.5))));
  	  			  Fref=FObj*((SObj-0.5)*(1-p2)/0.5 + p2);
  	  			  if(SObj<0.5){
  	  				  Fref=FObj*(SObj*(p2-p1)/0.5+p1);
  	  			  }
  	  		  }
  	  		  else{
  	  			  Fref=FObj;
  	  		  }
  	  	  }
  	  	  if(i==4) //10-40
  	  	  {
  	  		  if(SObj<0.9){
  	  			  Fref=FObj*(SObj-0.25)/(0.9-0.25);
  	  			  if(SObj<=0.25){
  					  Fref=0;}
  	  		  }
  	  		  else{
  	  			  Fref=FObj;
  	  		  }
  	  	  }
  		  F_fut(i,j)=Fref;
  		  Fmort = Sel(ntime)*Fref;
  		  Z_fut = Fmort+M;
  		  S_fut = exp(-1*Z_fut);
  		  CBAn = elem_prod(elem_div(Fmort,Z_fut),elem_prod(1-S_fut,N_fut));
  		  BD_fut(i,j) = sum(elem_prod(elem_prod(N_fut,exp(-dt(1)*Z_fut)),elem_prod(msex,Wmed)));
  		  //BDv_fut(i,j)= sum(elem_prod(elem_prod(Nv_fut*mfexp(-dt(3)*M),msex),colsum(Wmed)/nanos));;
  		  Y_fut(i,j)  = sum(elem_prod(CBAn,Wmed));
  		  ssb_fut = BD_fut(i,j);
  		  SObj = BD_fut(i,j)/(0.55*SO);
  		  BT_fut(i,j) = sum(elem_prod(N_fut,Wmed));
  		  Depl(i,j) = BD_fut(i,j)/SO;
  		}
  	}
    rep1 << BD_fut << endl;
  	rep2 << Y_fut << endl;
  	rep3 << BD << endl;
  	rep4 << F_fut << endl;
  	rep5 << exp(log_F) << endl;
  	rep6 << Depl << endl;
  	rep7 << BD/SO << endl;
  	rep8 << R_fut << endl;
  	rep9 << Reclutas << endl;
  */

FUNCTION Write_R

  	    wm =  Wmed;
  	  	sel = Sel(ntime);
  	  	adstring report_name;
  	  	{
  	  		report_name = "For_R_SprattusChiloe.rep";
  	  		ofstream R_report(report_name);
  	  		R_report << "yrs" << endl << yrs << endl;
  	  		R_report << "yld" << endl << Desemb << endl;
  	  		R_report << "indice1" << endl << Bcru << endl;
			R_report << "indice2" << endl << mph << endl;
  	  		R_report << "bt" << endl << BT << endl;
  	  		R_report << "rt" << endl << Reclutas << endl;
  	  		R_report << "ssb" << endl << BD << endl;
  	  		R_report << "ft" << endl << fyr << endl;
  	  		R_report << "CIssb" << endl;
  	  		for(int i=1;i<=ntime;i++){
  	  			double lb = value(BD(i)/exp(2.*sqrt(log(1+square(BD.sd(i))/square(BD(i))))));
  	  			double ub = value(BD(i)*exp(2.*sqrt(log(1+square(BD.sd(i))/square(BD(i))))));
  	  			R_report << yrs(i) << " " << BD(i) << " " << BD.sd(i) << " " << lb << " " << ub << endl;
  	  		}
  	  		R_report << "CIrt" << endl;
  	  		for(int i=1;i<=ntime;i++){
  	  			double lb = value(Reclutas(i)/exp(2.*sqrt(log(1+square(Reclutas.sd(i))/square(Reclutas(i))))));
  	  			double ub = value(Reclutas(i)*exp(2.*sqrt(log(1+square(Reclutas.sd(i))/square(Reclutas(i))))));
  	  			R_report << yrs(i) << " " << Reclutas(i) << " " << Reclutas.sd(i) << " " << lb << " " << ub << endl;			
  	  		}
  	  		R_report << "CIft" << endl;
  	  		for(int i=1;i<=ntime;i++){
  	  			double lb = value(fyr(i)/exp(2.*sqrt(log(1+square(fyr.sd(i))/square(fyr(i))))));
  	  			double ub = value(fyr(i)*exp(2.*sqrt(log(1+square(fyr.sd(i))/square(fyr(i))))));
  	  			R_report << yrs(i) << " " << fyr(i) << " " << fyr.sd(i) << " " << lb << " " << ub << endl;			
  	  		}
  	  		R_report << "m" << endl << M << endl;
  	  		R_report << "lh" << endl;
  	  		for(int j=1;j<=ntallas;j++){
  	  			R_report << Tallas(j) << " " << sel(j) << " " << wm(j) << " " << msex(j) <<  endl;
  	  		}
		
  	  	}



  	//=============== FIN Modificacion ===========================


//========================================================================

REPORT_SECTION

//========================================================================
 report << "Years" << endl;
 report << yrs << endl;
 report << "Bcru_obs" << endl;
 report << Bcru << endl;
 report << "Bcru_pred" << endl;
 report << pred_Bcru << endl;
 report << "CPUE_obs" << endl;
 report << CPUE << endl;
 report << "CPUE_pred" << endl;
 report << pred_CPUE << endl;
 report << "MPH_obs" << endl;
 report << mph << endl;
 report << "MPH_pred" << endl;
 report << pred_mph << endl;
 report << "Desemb_obs" << endl;
 report << Desemb << endl;
 report << "Desemb_pred" << endl;
 report << pred_Desemb << endl;
 report << "Lmf_obs" << endl;
 report << Lmed_obs << endl;
 report << "Lmf_pred" << endl;
 report << Lmed_pred << endl;
 report << "Lmc_obs" << endl;
 report << Lmed_obsc << endl;
 report << "Lmc_pred" << endl;
 report << Lmed_predc << endl;
 report << "Biomasa_desovante" << endl;
 report << BD << endl;
 report << "Biomasa_total" << endl;
 report << BT << endl;
 report << "Biomasa_explotable" << endl;
 report << BMflo << endl;
 report << "Reclutamiento" << endl;
 report << Reclutas<< endl;
 report << "Rpred" << endl;
 report << Rpred<< endl;
 report << "F" << endl;
 report << exp(log_F) << endl;
 report<<"Tallas"<<endl;
 report<<Tallas<<endl;
 rep(SO)
 rep(RO)
 rep(h)
 rep(SSBo)
 rep(alfa_sr)
 rep(beta_sr);
 rep(fpbr)
 rep(SPR0)
 rep(SPRf)
 rep(ratio_spr)
 rep(pre)
 
 report<<"Abundancia_talla"<<endl;
 report<<N<<endl;
 report<<"Selflo_talla"<<endl;
 report<<Sel<<endl;
 report<<"Selcru_talla"<<endl;
 report<<Selc<<endl;
 report << "Propfl_obs" << endl;
 report << pobs<< endl;
 report << "Propfl_pred" << endl;
 report << ppred<< endl;
 report << "Propcru_obs" << endl;
 report << pobsc<< endl;
 report << "Propcru_pred" << endl;
 report << ppredc<< endl;
 report << "BD_virgen_anual" << endl;
 report << BDo << endl;
 report << "BD_virgen_LP" << endl;
 report << SSBo << endl;
 report << "Reduccion_LP " << endl;
 report << RPRlp << endl;
 
 report << "Talla_media_por_grupo" << endl;
 report << Lesp << endl;
 report <<  "desvest_por_grupo" << endl;
 report << sigmaL << endl;
 report << "Fun_rec_talla" << endl;
 report << pre<< endl;
 report << "MatrizTrans" << endl;
 report << T << endl;
 report << "bCPUE  Lr  Sr  beta  h " << endl;
 report << exp(log_b)<<" "<<exp(log_Lr)<<" "<<exp(log_sr)<<" "<<exp(log_beta)<<" "<<h<< endl;
 report << "pred_Ctot" << endl;
 report << pred_Ctot << endl;

//-------------------------------------------------------------------
// ESTIMA nm y CV

  suma1=0; suma2=0;nm1=1;cuenta1=0;

  for (int i=1;i<=ntime;i++){ //

   if (sum(pobs(i))>0){
      suma1=sum(elem_prod(ppred(i),1-ppred(i)));
      suma2=norm2(pobs(i)-ppred(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}

 report << "nm_flota_cru" <<endl;
 report <<pow(nm1,1/cuenta1)<< endl;


 suma1=0; suma2=0;nm1=1;cuenta1=0;

  for (int i=1;i<=ntime;i++){ //

   if (sum(pobs(i))>0){
      suma1=sum(elem_prod(ppredc(i),1-ppredc(i)));
      suma2=norm2(pobsc(i)-ppredc(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}

 report <<pow(nm1,1/cuenta1)<< endl;
 report << "BD_proy" << endl; //biomasa desovante proyectada para cada multiplo de Flast" 
 report << Bp << endl;
 report << "Capt_proy" << endl; // Capturas proyectadas para cada Fpbr
 report << Yp << endl;
 report << "Captura_act" << endl;  //Captura proyectadas año en curso
 report << Yact << endl;
 report<<"Likeval CPUE_BCru_Bmph_Desemb_pfFlota_pfCru_dev_R_devNo_LR"<<endl;
 report << likeval << endl;
 report << "Prioris" << endl;
 report << penalty << endl;
 report << "Rec_proy" << endl; //Reclutamientos proyectadas para cada Fpbr
 report << Rpp << endl;
 report << "N_actual" << endl; // abundancia del ultimo año
 report << Npact << endl;
 report << "N_proyect" << endl; // abundancia del ultimo año
 report << Np << endl;
 report << "Rproy1" << endl; //Rec
 report << mfexp(log_Rmed) << endl;
 report << "Rproy2" << endl; //Rec
 report << Reclutas(11)<< endl;
 report << "Rproy3" << endl; //Rec
 report << Reclutas(17)<< endl;
 
 report << "NV" << endl;
 report << NV << endl;
 report << "CTP" << endl;
 report << CTP << endl;
 report << "NVp" << endl;
 report << NVp << endl;
 report << "CTPp" << endl;
 report << CTPp << endl;

 
