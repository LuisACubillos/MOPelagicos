//##############################################################################
// MODELO ANUAL EN EDADES ANCHOVETA V-X REGIONES
// modelo CTP 2019
//##############################################################################
GLOBALS_SECTION
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc.csv");
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
  
 #undef rep
 #define rep(object) report << #object "\n" << object << endl;	
 #undef depur
 #define depur(object) cout << #object "\n" << object << endl; exit(1);


//##############################################################################
TOP_OF_MAIN_SECTION
//##############################################################################
 time(&start);
 arrmblsize=60000000; // 
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(30000000);
 gradient_structure::set_CMPDIF_BUFFER_SIZE(50000000);
 gradient_structure::set_MAX_NVAR_OFFSET(100000);

//##############################################################################

DATA_SECTION

//##############################################################################
  int iseed
  !!long int lseed=iseed;
  !!CLASS random_number_generator rng(iseed);


// LEER DATOS "MAEanchoveta.dat"
  init_int nanos  
  init_int nedades
  init_int ntallas
  init_matrix matdat(1,nanos,1,9)
  init_matrix Ctot(1,nanos,1,nedades)
  init_matrix Ccru_a(1,nanos,1,nedades)
  init_matrix Ccru_pel(1,nanos,1,nedades)
  init_matrix Ccru_l(1,nanos,1,ntallas)
  init_matrix Wmed(1,nanos,1,nedades)
  init_matrix Win(1,nanos,1,nedades)  
  int reporte_mcmc

//==============================================================================
// LEER controles y opciones
 !! ad_comm::change_datafile_name("controles_op.ctl");
//==============================================================================
// 1. Coeficientes de variación y tamaños de muestra  
  init_matrix error_edad(1,nedades,1,nedades)
  init_vector Wmedp_3(1,5)   //vector de pesos medios promedio de los últimos 5 años de la serie
  init_vector Winip_3(1,5)   // vector de pesos iniciales promedio de los últimos 5 años de la serie  
  init_vector edades(1,nedades)
  init_vector Tallas(1,ntallas)
  init_vector msex(1,nedades) 	  
  init_number sigmaR
  init_number cvpriorq_reclas
  init_number cvpriorq_pelaces
  init_number log_priorRo 
  //init_int fase_Ro
  init_vector nmus(1,4)
  init_vector dt(1,4)
// 2. Fases de selectividad y penalizadores de selectividad no parametrica
  init_int    opt1_fase
  init_int    opt2_fase
  init_int    opt_Scru1
  init_int    opt_Scru2
  init_vector lam(1,4)
// 3. Fases de capturabilidad
  init_int    opt_qrecl
  init_int    opt_qpela
  init_int    opt_qmph
// 4. Parámetros de crecimiento
  init_vector pars_Bio(1,5)

// 7. Fases de estimación Lo y cv edad
  init_int    opt_Lo
  init_int    opt_cv
// 8. Considera la matriz de asignación de error edad
  init_number erredad 
// 9. Fase de estimación de M
  init_int    opt_M
// 10. Fase de estimación condiciones iniciales
  init_int    opt_Ro
  init_int    opt_devR
  init_int    opt_devNo
// 11. Fase de estimación de F
  init_int    opt_F
// 12. Puntos biológicos de referencia
  init_number  Fmedian_ext
  init_int    opt_Fspr       // fase de estimación
  init_int    npbr
  init_vector ratio(1,npbr) //cambio (1,npbr) a (1,2) por que cambie el numero de pbr
 // 13. PROYECCION, CRITERIOS DE EXPLOTACIÓN Y ESCENARIOS DE RECLUTAMIENTO
  init_int nproy
  init_number opProy         // define desde donde proyecto
  init_number oprec          // define escenario de reclutamiento
  init_vector prop(1,2)      // proporción semestral de la captura 
  init_number mF             // multiplicadores de Frms para determinar período de recuperación
  init_number Desembactual
  init_vector prop_est(1,5)
  init_number h
  init_number logSo
  init_int SrType 
  init_int ErrorObs  
	  //!!cout << "SrType " << SrType << endl; exit(1);
//##############################################################################

INITIALIZATION_SECTION

//##############################################################################

// defino un valor inicial de log_reclutamiento promedio (factor de escala)
  log_Ro      12.54
  log_Lo      2
  log_cv_edad -2.52
  log_M       0
  log_qrecl   0
  log_qpela   0

//##############################################################################

PARAMETER_SECTION

//##############################################################################
// selectividad paramétrica
 init_bounded_number A50f(-1,2,opt1_fase)  
 init_bounded_number log_rangof(-4,0,opt1_fase)

 init_bounded_number A50c(-1,2,opt_Scru1)  
 init_bounded_number log_rangoc(-4,0.6,opt_Scru1)

 init_bounded_number A50pela(-1,2,opt1_fase)  
 init_bounded_number log_rangopela(-4,0.6,opt1_fase)

// selectividad no_paramétrica
 init_bounded_vector log_Sflo(1,nedades,-11,0,opt2_fase)  
 init_bounded_vector log_Scru(1,nedades,-11,0,opt_Scru2)
 init_bounded_vector log_Scru_pela(1,nedades,-11,0,opt_Scru2)

// parametros reclutamientos y mortalidades)
 init_bounded_number log_Ro(5,20,opt_Ro)
 init_bounded_vector log_desv_No(1,nedades-1,-10,10,opt_devNo)
 init_bounded_dev_vector log_desv_Rt(1,nanos,-10,10,opt_devR)
 init_bounded_vector log_Ft(1,nanos,-6,1.6,opt_F) // log  mortalidad por pesca por flota

// capturabilidades
 init_number log_qrecl(opt_qrecl)
 init_number log_qpela(opt_qpela)
 init_number log_qmph(opt_qmph)

// Crecimiento y M
 init_bounded_number log_Lo(1,2.4,opt_Lo)
 init_bounded_number log_cv_edad(-4,.69,opt_cv)
 init_bounded_number log_M(-0.3,0.4,opt_M)
 init_bounded_vector log_Fref(1,npbr,0.01,2.,opt_Fspr)
//##############################################################################

// VARIABLES DE ESTADO

//##############################################################################
  // ARREGLOS
  vector anos(1,nanos);
  vector Unos_edad(1,nedades);
  vector Unos_anos(1,nanos);
  vector Unos_tallas(1,ntallas);
  //=====================================
  // 1. SELECTIVIDAD
  //=====================================
  matrix Sel_f(1,nanos,1,nedades)
  matrix Scru(1,nanos,1,nedades)
  matrix Scru_pela(1,nanos,1,nedades)
  vector pSf(1,nedades-2);
  vector pSc(1,nedades-2);
  vector pStotf(1,nanos);
  vector pStotc(1,nanos);
  //=====================================
  // 2. PROBABILIDAD EDAD-TALLA
  //=====================================
  number Linf
  number k
  number Lo
  number cv_edad
  vector mu_edad(1,nedades);
  vector sigma_edad(1,nedades);
  matrix Prob_talla(1,nedades,1,ntallas)
  matrix P1(1,nedades,1,ntallas)
  matrix P2(1,nedades,1,ntallas)
  matrix P3(1,nedades,1,ntallas)
  //=====================================
  // 3. MORTALIDADES
  //=====================================
  number M
  matrix Ftot(1,nanos,1,nedades)
  matrix Z(1,nanos,1,nedades)
  matrix S(1,nanos,1,nedades)
  //=====================================
  // 4. ABUNDANCIAS Y BIOMASAS
  //=====================================
  // 4.1. Abundancia inicial en equilibrio
  vector Neq(1,nedades);
  number SSBo;
  // 4.2. Abundancia
  matrix N(1,nanos,1,nedades)
  matrix NM(1,nanos,1,nedades)
  // 4.3. Matrices y vectores abundancia derivadas
  matrix NVflo(1,nanos,1,nedades)
  matrix NVcru(1,nanos,1,nedades)
  matrix NVpel(1,nanos,1,nedades)
  matrix NVpel_l(1,nanos,1,ntallas)
  matrix NMD(1,nanos,1,nedades)
  sdreport_vector Reclutas(1,nanos);
  // 4.4. Vectores de biomasa derivadas
  sdreport_vector SSB(1,nanos) 
  vector BD(1,nanos);
  sdreport_vector BT(1,nanos) 
  vector BMflo(1,nanos);
  vector BMpel(1,nanos);
  vector Bcru(1,nanos);
  //=====================================
  // 5. CAPTURAS observadas y predichas
  //=====================================
  // 5.1. Matrices de capturas predichas por edad y años
  matrix pred_Ctot(1,nanos,1,nedades)
  // 5.2. Matrices de proporción de capturas por edad y años
  matrix pobs_f(1,nanos,1,nedades)
  matrix ppred_f(1,nanos,1,nedades)
  matrix pobs_crua(1,nanos,1,nedades)
  matrix ppred_crua(1,nanos,1,nedades)
  matrix pobs_pel(1,nanos,1,nedades)
  matrix ppred_pel(1,nanos,1,nedades)
  matrix pobs_crul(1,nanos,1,ntallas)
  matrix ppred_crul(1,nanos,1,ntallas)
  //==================================================
  // 6. INDICES DE ABUNDANCIA observadas y predichas
  //==================================================
  number qrecl
  number qpela
  vector Reclas(1,nanos);
  vector Reclas_pred(1,nanos);
  vector Pelaces(1,nanos);
  vector Pelaces_pred(1,nanos);
  vector MPH(1,nanos);
  vector MPH_pred(1,nanos);
  vector Desemb(1,nanos);
  vector Desemb_pred(1,nanos);
  //==================================================
  // 7. REDUCCION DE LA BIOMASA DESOVANTE
  //==================================================
  matrix Nv(1,nanos,1,nedades)
  matrix NDv(1,nanos,1,nedades)
  matrix NMDv(1,nanos,1,nedades)
  vector BDo(1,nanos);
  vector RPRdin(1,nanos)  // 
  vector RPRequ(1,nanos)  // 
  vector RPRequ2(1,nanos) // 
  sdreport_vector RPRequ3(1,nanos) // 
  sdreport_vector Frpr(1,nanos)
  //==================================================
  // 7.1. PBRS
  //==================================================
  vector Fspr(1,nedades)
  vector Zspr(1,nedades)
  vector Nspro(1,nedades)
  vector Nspr(1,nedades)
  vector Nmed(1,nedades)
  number Bspro
  number Bspr
  vector ratio_spr(1,npbr)
  
  number Fmedian   
  vector Fmed(1,nedades)
  vector Zmed(1,nedades)
  number Bsprmed
  number ratio_Fmed
  number Bmed
  number Bo
  sdreport_vector Brms(1,npbr)
  number Blim
  //==================================================
  // 8. LOGVEROSIMILITUD
  //==================================================
  matrix cvar(1,4,1,nanos)
  vector likeval(1,15);
  objective_function_value f
  
  //==================================================
  // 9. CBA año biológico sin proyectar
  //==================================================
  number Fref_r0
  vector Frms_r0(1,nedades);
  vector Zrms_r0(1,nedades);
  vector CTP_r0(1,nedades);
  number YTP_r0
  vector NMD_r0(1,nedades);
  number BD_r0
  number RPR_r0
  //==================================================
  // 10. PROYECCION DEL STOCK
  //==================================================
  
  vector Np(1,nedades);
  vector Sp(1,nedades);
  vector Nvp(1,nedades);
  vector RPRp(1,nproy)
  vector Npp(1,nedades); 
  vector Wmedp(1,nedades)
  vector Winp(1,nedades)        
  
  number Fref_p0
  vector Frms_p0(1,nedades);
  vector Zrms_p0(1,nedades);
  vector CTP_p0(1,nedades);
  vector YTP_p0(1,nproy);
  sdreport_vector BD_p0(1,nproy);
  sdreport_vector RPR_p0(1,nproy);
 
  vector NVrecl_p0(1,nedades);
  vector NVpel_p0(1,nedades);
  vector Brecl_p0(1,nproy)
  vector Bpel_p0(1,nproy)
  number BDp;
  number BTp;
  vector Zpbr(1,nedades);
	  
  number alpha;
  number SPRO;
  number RO;
  number SO;  

  sdreport_number CBA_c0
    
  number suma1
  number suma2
  number suma3
  number suma4
  number nm1
  number nm2
  number nm3
  number nm4
  number alfa
  number beta

  number cuenta1
  number cuenta2
  number cuenta3
  number cuenta4

//###########################################################################  
// VARIABLES PARA EL MODELO OPERATIVO
// Proyeccion futura
	  vector yrs_fut(nanos+1,nanos+nproy);
	  //number alpha
	  //number beta
	  number Kobs_tot_catch;
	  matrix catch_future(1,npbr,nanos+1,nanos+nproy);
	  matrix future_biomass(1,npbr,nanos+1,nanos+nproy);//matrix o vector??
	  vector expl_biom(nanos+1,nanos+nproy);	  
	  matrix future_ssbiom(1,npbr,nanos+1,nanos+nproy);
	  vector N_fut(1,nedades);
	  matrix F_fut(nanos+1,nanos+nproy,1,nedades);
	  //matrix Zpbr(nanos+1,nanos+nproy,1,nedades);
	  //vector Zpbr(2,nedades);
	  //vector Sp(1,nedades);
	  matrix check_convergence(1,npbr,nanos+1,nanos+nproy);
	  matrix Fyr_future(1,npbr,nanos+1,nanos+nproy);
	  vector NVflo_fut(1,nedades);
	  vector NVcru_fut(1,nedades);
	  vector NVpel_fut(1,nedades);
	  matrix F_tmp(nanos+1,nanos+nproy,1,nedades);
      vector BMflo_fut(nanos+1,nanos+nproy);
	  vector Bcru_fut(nanos+1,nanos+nproy);
	  vector BMpel_fut(nanos+1,nanos+nproy);
	  vector BMmph_fut(nanos+1,nanos+nproy);
	  vector NMD_fut(1,nedades);
	  //number BMflo_fut;
	 //Datos simulados
	  matrix catage_fut(nanos+1,nanos+nproy,1,nedades);
	  matrix reclas_age_fut(nanos+1,nanos+nproy,1,nedades);
	  matrix pelaces_age_fut(nanos+1,nanos+nproy,1,nedades);
	  matrix pelaces_len_fut(nanos+1,nanos+nproy,1,ntallas);
	  vector wm(1,nedades);
	  vector Wm_fut(1,nedades);
	  vector Win_fut(1,nedades);
	  vector error_edad_fut(1,nedades);	
   
	  //vector Win_fut(1,nanos);	
	  //matrix Scru_fut(nanos+1,nanos+nproy,1,nedades);
	  vector Scru_fut(1,nedades);
	  vector Scru_pela_fut(1,nedades);
	  //matrix Scru_pela_fut(nanos+1,nanos+nproy,1,nedades);
	  vector Sel_f_fut(1,nedades);
	  vector sel(1,nedades);
	  //vector Np(1,nanos);
	  vector rec_epsilon_fut(nanos+1,nanos+nproy);
	  vector Reclas_epsilon_fut(nanos+1,nanos+nproy);
	  vector Pelaces_epsilon_fut(nanos+1,nanos+nproy);
	  vector mpdh_epsilon_fut(nanos+1,nanos+nproy);
	  vector sim_Reclas(nanos+1,nanos+nproy);
	  vector sim_Pelaces(nanos+1,nanos+nproy);
	  vector sim_MPDH(nanos+1,nanos+nproy);
	  matrix rpr_fut(1,npbr,nanos+1,nanos+nproy);
	  matrix des_Btot(1,npbr,nanos+1,nanos+nproy);
	  matrix des_SSBt(1,npbr,nanos+1,nanos+nproy);
	  matrix des_Rt(1,npbr,nanos+1,nanos+nproy);
	  matrix des_F(1,npbr,nanos+1,nanos+nproy);
	  matrix des_Yt(1,npbr,nanos+1,nanos+nproy);
	  matrix des_rpr(1,npbr,nanos+1,nanos+nproy);
	  matrix keep_Btot(1,npbr,nanos+1,nanos+nproy);
	  matrix keep_SSBt(1,npbr,nanos+1,nanos+nproy);
	  matrix keep_Rt(1,npbr,nanos+1,nanos+nproy);
	  matrix keep_F(1,npbr,nanos+1,nanos+nproy);
	  matrix keep_RBB(1,npbr,nanos+1,nanos+nproy);
	  matrix des_RBB(1,npbr,nanos+1,nanos+nproy);


PRELIMINARY_CALCS_SECTION

//##############################################################################
// leo la matriz de indices
  anos=column(matdat,1);// asigno la 1 columna de indices a "anos"
  Reclas=column(matdat,2);
  cvar(1)=column(matdat,3);
  Pelaces=column(matdat,4);
  cvar(2)=column(matdat,5);
  MPH=column(matdat,6);
  cvar(3)=column(matdat,7);
  Desemb=column(matdat,8);
  cvar(4)=column(matdat,9);

  Unos_edad=1;;// lo uso en  operaciones matriciales con la edad
  Unos_anos=1;// lo uso en operaciones matriciales con el año
  Unos_tallas=1;// lo uso en operaciones matriciales con el año
  reporte_mcmc=0;

//##############################################################################

RUNTIME_SECTION

//##############################################################################

  maximum_function_evaluations 200,1000,5000
  convergence_criteria  1e-3,1e-5,1e-6

//##############################################################################

PROCEDURE_SECTION

//##############################################################################
// se listan las funciones que contienen los calculos

  if(opt1_fase>0)
  {
      Eval_selectividad_logis();
  }

  if(opt2_fase>0)
  {
      Eval_selectividad_nopar();
  }

  Eval_prob_talla_edad();
  Eval_mortalidades();
  Eval_abundancia();
  Eval_biomasas();
  Eval_SRR(); //agrega funcion para parametros relacion stock-recluta
  Eval_capturas_predichas();
  Eval_indices();
  Eval_PBR();
  Eval_Estatus();
  Eval_logverosim();
  Eval_funcion_objetivo();
  //Eval_CTP();
  //Eval_mcmc();
  if(mceval_phase()){ //se ejecuta solo cuando se corre mcmc
  	   Eval_Mop();
    }

 //Con fines de hacer el MOP en la última fase	
 //FINAL_SECTION
	     //Write_R();
	//Eval_Mop();

FUNCTION Eval_Mop

	    //Error de Proceso
	   	dvector epsilon_rec_vect(nanos+1,nanos+nproy);
        //Error de observacion
	    dvector epsilon_Reclas(nanos+1,nanos+nproy);
	   	dvector epsilon_Pelaces(nanos+1,nanos+nproy);
	   	dvector epsilon_MPDH(nanos+1,nanos+nproy);
	   	int upk;//actualizar los años contador para cada ciclo
	   	int k,l,j,i,a;//contador de los ciclos for
		dvariable numyear=1997;
		dvector yrs(nanos+1,nanos+nproy);
		simname="MAE0920b.dat";//datos del estimador
	   	dvector CatchNow(1,npbr);//capturas del estimador
	   	dvector eval_now(1,5);//Indicadores de estatus del estimador:BT,BD,R,F,B/Bo
	   	dvariable grad_tmp; //Gradiente
		dvariable sigr=0.6; //Sigma R para el operativo (puede ser menor que del estimador)
		epsilon_rec_vect.fill_randn(rng);//error de proceso
		epsilon_Reclas.fill_randn(rng);//errores de observacion
		epsilon_Pelaces.fill_randn(rng);
		epsilon_MPDH.fill_randn(rng);
		int rv;
		int indice_cuota; //Indice para obtener ya sea la cuota inicial o revisada para el estimador, SOLO revision incluyendo
		
	    //cout << "ctp"<< CatchNow << endl;exit(1);
	
	 	for(int l=1;l<=npbr;l++){
	 		//int rv=system("cp MAE0920_real.dat MAE0920b.dat");
	 		//rv=system("./MAE0920b -nox -nohess");
			indice_cuota = 1;
			rv = system("./paso1.sh");
	 		N_fut=N(nanos);
	 		Sp=S(nanos);
	 		BDp=BD(nanos);
	 		Zpbr=Z(nanos);
	 		Win_fut=Win(nanos);
			Wm_fut=Wmed(nanos);//asume el peso medio del año mas reciente
	 	    Scru_fut=Scru(nanos);
	 	    Scru_pela_fut=Scru_pela(nanos);
	 	    sel=Sel_f(nanos);

	 		for(i=nanos+1;i<=nanos+nproy;i++)
	 		{
	 			rec_epsilon_fut(i)=epsilon_rec_vect(i)*sigr; //Modelo operativo con sigma_r definido en linea 524 
	 			Reclas_epsilon_fut(i)=epsilon_Reclas(i)*0.3;
	 			Pelaces_epsilon_fut(i)=epsilon_Pelaces(i)*0.3;
	 			mpdh_epsilon_fut(i)=epsilon_MPDH(i)*0.3; //????
	 	        //lee la cuota inicial del estimador
	 	        ifstream CTP_tmp("ctp.dat");
	 	        CTP_tmp >> CatchNow;
	 	        CTP_tmp.close();
	 			catch_future(l,i) = CatchNow(l);
				//Sobrevivencia
	 			N_fut(2,nedades)= ++elem_prod(N_fut(1,nedades-1),Sp(1,nedades-1));
	 			if(SrType==1){
	 				N_fut(1)=(BDp/(alpha + beta*BDp))*mfexp(rec_epsilon_fut(i) + 0.5*square(sigr));
	 			}
	 			if(SrType==2){
	 				N_fut(1)=(alpha*BDp*mfexp(-beta*BDp))*mfexp(rec_epsilon_fut(i) + 0.5*square(sigr));
	 			}
		        //Seccion para estimar la F que genera la CBA inicial leida desde el estimador		   
		        NVflo_fut= elem_prod(N_fut,sel);
	          	//cout << "HASTA AQUI= " << NVflo_fut << endl;exit(1);	   
				expl_biom(i)=0;
				for(j=1;j<=nedades;j++){
					//expl_biom(i)+=Sel_f_fut(i,j)*N_fut(i,j)*Wm_fut(i,j);
					expl_biom(i)+=Wm_fut(j)*NVflo_fut(j);
		         }
				 
				 //cout << "BDp: " << BDp << endl;
				 //cout << "Rt: " << N_fut(1) << endl; exit(1);

		        if(catch_future(l,i)!=0)
				{
					dvariable ffpen=0.0;
					dvariable SK=posfun((expl_biom(i)-catch_future(l,i))/expl_biom(i),0.05,ffpen);
					Kobs_tot_catch=expl_biom(i)-SK*expl_biom(i);
					do_Newton_Raphson_for_mortality(i);	
					F_fut(i)=F_tmp(i);
				    }else{
					F_fut(i)=0;
				    }		
	            Zpbr=F_fut(i)+M;
		        Sp=exp(-1.*Zpbr);
	         	
				//Recalcula los indicadores despues que la mortalidad fue calculada 
		        NVflo_fut = elem_prod(elem_prod(N_fut,mfexp(-dt(4)*Zpbr)),sel); //explotable flota		
		        BMflo_fut(i) = sum(elem_prod(Win_fut,NVflo_fut));
		        //Crucero 
		        NVcru_fut = elem_prod(elem_prod(N_fut,mfexp(-dt(1)*Zpbr)),Scru_fut); 
				//Pelaces
				NVpel_fut = elem_prod(elem_prod(N_fut,mfexp(-dt(2)*Zpbr)),Scru_pela_fut); 
			    if(erredad==1)
			    {
			    NVcru_fut = NVcru_fut*error_edad;
				NVpel_fut =	NVpel_fut*error_edad;
			    }
				Bcru_fut(i)=sum(elem_prod(Wm_fut,NVcru_fut));
				BMpel_fut(i) = sum(elem_prod(Win_fut,NVpel_fut));		
		        //Desovante
		        NMD_fut = elem_prod(elem_prod(N_fut,mfexp(-dt(3)*(Zpbr))),msex);	//numero desovante
		        BMmph_fut(i)=sum(elem_prod(Win_fut,NMD_fut));
					
		        BDp=sum(elem_prod(NMD_fut,Win_fut));	
		        BTp=sum(elem_prod(N_fut,Win_fut));	
					
		        //======== SIMULA LOS DATOS PARA ACTUALIZAR LA EVALUACION CON EL ESTIMADOR =====
		        if(ErrorObs==1){
					sim_Reclas(i) = exp(log_qrecl)*Bcru_fut(i);
					sim_Pelaces(i) = exp(log_qpela)*BMpel_fut(i);
					sim_MPDH(i) = exp(log_qmph)*BMmph_fut(i);
					}
				if(ErrorObs==2){
					sim_Reclas(i) = exp(log_qrecl)*Bcru_fut(i)*mfexp(Reclas_epsilon_fut(i));
					sim_Pelaces(i) = exp(log_qpela)*BMpel_fut(i)*mfexp(Pelaces_epsilon_fut(i));
					sim_MPDH(i) = exp(log_qmph)*BMmph_fut(i)*mfexp(mpdh_epsilon_fut(i));
				}				
		        //Ahora simula la composicion por edad y por talla
		        for(j=1;j<=nedades;j++){
		            catage_fut(i,j)   = (F_fut(i,j)/Zpbr(j))* N_fut(j)*(1.-Sp(j));	
		            reclas_age_fut(i,j) = NVcru_fut(j);			
		            pelaces_age_fut(i,j)=NVpel_fut(j);
				}
				pelaces_len_fut(i) = pelaces_age_fut(i)*Prob_talla;//SACAR DEL  
		        //cout << "HASTA AQUI5= " << NVflo_fut << endl;exit(1);	   
					
		upk=i;
		yrs_fut(i)=numyear+i-1;
		ofstream simdata(simname);
		simdata <<"#SIMULACION DE DATOS"<<endl;
		simdata << "#Numero de anos actual" << endl;
		simdata << upk << endl;
		simdata << "#Numero de edades" << endl;
		simdata << nedades << endl;
		simdata << "#Numero de tallas" << endl;
		simdata << ntallas << endl;
		simdata << "#Datos observados"<<endl;
		simdata << "#Yr Reclas cv_reclas Pelaces cv_pelaces Mpdh cv_mpdh Desemb cv_desemb  "<<endl;
		simdata << matdat << endl;//DUDA
		simdata << "#Agrega datps simulados a la matrix de datos historicos"<<endl;
		for(k=nanos+1;k<=upk;k++){
		simdata << yrs_fut(k) << " " << sim_Reclas(k) << " " << 0.3 << " " << sim_Pelaces(k) << " " << 0.3 << "  " << sim_MPDH(k) << " " << 100 << " " << catch_future(l,k) << " " << 0.01 <<  endl; 	
		}
		
		simdata << "#Estructura de tallas flota historica"<< endl;
		simdata << Ctot << endl;
		simdata << "# flota simulados" << endl;
		for(k=nanos+1;k<=upk;k++){
		simdata << catage_fut(k) << endl;	
	    }
		simdata << "#Reclas age historico"<<endl;
		simdata << Ccru_a << endl;
		simdata << "# Reclas age simulados" << endl;
		for(k=nanos+1;k<=upk;k++){
		simdata << reclas_age_fut(k) << endl;
		}	
		simdata << "#Pelaces age historico"<<endl;
		simdata << Ccru_pel << endl;
		simdata << "# Pelaces age simulados" << endl;
		for(k=nanos+1;k<=upk;k++){
		simdata << pelaces_age_fut(k) << endl;
		}	
		simdata << "#Pelaces length historico"<<endl;
		simdata << Ccru_l << endl;
		simdata << "#Pelaces length simulado"<<endl; 
		for(k=nanos+1;k<=upk;k++){
		simdata << pelaces_len_fut(k) << endl;
		}
		simdata << "#Peso medio observados" << endl;
		simdata << Wmed << endl;
		simdata << "#Peso medio simulacion"<< endl;
	
		for(k=nanos+1;k<=upk;k++){
		simdata << Wm_fut << endl;
	    }
		simdata << "#Peso inicial observados" << endl;
		simdata << Win << endl;
		simdata << "#Peso inicial simulacion"<< endl;
		for(k=nanos+1;k<=upk;k++){
		simdata << Win_fut << endl;	
	    }	
		
			 
		simdata.close();		
							
		rv=system("./paso2.sh");
		//Lee y retiene los indicadores del estimador en la evaluación actual
		ifstream eval_estatus("estatus.dat");
	
		eval_estatus >> eval_now;
		eval_estatus.close();		
		keep_Btot(l,i)=eval_now(1); //biomasa total
		keep_SSBt(l,i)=eval_now(2); //biomasa desovante
		keep_Rt(l,i) =eval_now(3); //reclutas
		keep_F(l,i) =eval_now(4); //Mortalidad por pesca
		keep_RBB(l,i)=eval_now(5);
		//  
		//Guarda los indicadores del simulador que debío estimar el estimador (valga la redundancia)
				
		des_Btot(l,i)=BTp; //Biomasa total
		des_SSBt(l,i)=BDp; //Biomasa desovante
		des_Rt(l,i) = N_fut(1); //Reclutamiento 
		des_F(l,i) = max(F_fut(i)); //Mortalidad por pesca 
		des_RBB(l,i) = BDp/Bo;
	
		ifstream lee_grad("grad_final.dat");
		lee_grad >> grad_tmp;
		lee_grad.close();
		if(grad_tmp>0.001){check_convergence(l,i)=0;}else{check_convergence(l,i)=1;}
		} //cierra ciclo de proyeccion
		//cout << "SSB =" << des_SSBt << endl;
		//cout << "Rt =" << des_Rt << endl; exit(1);
	    } //cierra ciclo del numero de HCR		
 
	 //======== REPORTA SIMLACIONES Y ESTIMACIONES
	 //Capturas efectivas		       
	    rep1 << catch_future << endl;
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
	    rep10 << mfexp(log_Ft) << endl; //mortalidad pesca         
	    rep11 << BD << endl; //desovante
	    rep12 << BT << endl; //biomasa total
	    rep13 << Reclutas << endl; //Reclutas
		rep16 << BD/Bo <<endl;
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
				F_tmp(i)= Fnew*sel; 
				
	
	
  
//===============================================================================

FUNCTION Eval_selectividad_logis

//===============================================================================
   Sel_f     = outer_prod(Unos_anos,(elem_div(Unos_edad,(1+exp(-1.0*log(19)*(edades-A50f)/exp(log_rangof))))));
   Scru_pela = outer_prod(Unos_anos,(elem_div(Unos_edad,(1+exp(-1.0*log(19)*(edades-A50pela)/exp(log_rangopela))))));
    
   if (opt_Scru1>0)// evaluo si el indice es >0
    {
       Scru  = outer_prod(Unos_anos,(elem_div(Unos_edad,(1+exp(-1.0*log(19)*(edades-A50c)/exp(log_rangoc))))));
    }    
   else
    {
       Scru  = 1.;
    }

//===============================================================================

FUNCTION Eval_selectividad_nopar

//===============================================================================
  int i;

  for (i=1;i<=nanos;i++)
   {
          Sel_f(i)(1,nedades)=mfexp(log_Sflo);

      if (opt_Scru2>0)
      {
          Scru(i)(1,nedades)=mfexp(log_Scru);
          Scru_pela(i)(1,nedades)=mfexp(log_Scru_pela);
      }
      else
      {
          Scru=1.;
      }
// controla que S(older)<=S(older-1)
      if (Sel_f(i)(nedades)>Sel_f(i)(nedades-1))
       {
           Sel_f(i)(nedades)=Sel_f(i)(nedades-1);
       }
      if (Scru(i)(nedades)>Scru(i)(nedades-1))
       {
           Scru(i)(nedades)=Scru(i)(nedades-1);
       }
   }

//===============================================================================

FUNCTION Eval_prob_talla_edad

//===============================================================================
// se supone proporcionalidad entre la integral de una pdf y la pdf

  Linf    = pars_Bio(1);
  k       = pars_Bio(2);
  Lo      = pars_Bio(3);
  cv_edad = pars_Bio(4);

  if(active(log_Lo))
  {
   Lo=exp(log_Lo);
  }
  			if(active(log_cv_edad))
  				{
     				cv_edad=exp(log_cv_edad);
  				}

// genero una clave edad-talla para otros calculos. Se modela desde L(1)
  int i, j;
  mu_edad(1)=Lo;

  for (i=2;i<=nedades;i++)
   {
      mu_edad(i)=Linf*(1-exp(-k))+exp(-k)*mu_edad(i-1);
   }
      sigma_edad=cv_edad*mu_edad;

  for (i=1;i<=nedades;i++)
   {
       P1(i)=(Tallas-mu_edad(i))/sigma_edad(i); //  standarizated deviation of the length (respect to expected length-at-age)

     for (j=1;j<=ntallas;j++) // cumulative pdf (normal distribution)
     {
       P2(i,j)=cumd_norm(P1(i,j));
     }
   } 
   
  for (i=1;i<=nedades;i++)// estimation of probabilities
  {
     for (j=2;j<=ntallas;j++)
     {
       P3(i,j)=P2(i,j)-P2(i,j-1);
     }
  } 
  /*
  P1=elem_div(outer_prod(Unos_edad,Unos_tallas),sqrt(2*3.1416)*outer_prod(sigma_edad,Unos_tallas));
  P2=mfexp(elem_div(-square(outer_prod(Unos_edad,Tallas)-outer_prod(mu_edad,Unos_tallas)),2*square(outer_prod(sigma_edad,Unos_tallas))));
  P3=elem_prod(P1,P2);
  */
  Prob_talla = elem_div(P3+1e-16,outer_prod(rowsum(P3+1e-16),Unos_tallas));// normalizo para que la suma sobre las edades sea 1.0
// cout<<Prob_talla<<endl;exit(1);
//===============================================================================

FUNCTION Eval_mortalidades

//===============================================================================

  if(opt_M>0)
  {
  M = exp(log_M);
  }
  else
  {
  M = pars_Bio(5);
  }

  Ftot = elem_prod(Sel_f,outer_prod(mfexp(log_Ft),Unos_edad));
  Z    = Ftot+M;
  S    = mfexp(-1.0*Z);

//===============================================================================

FUNCTION Eval_abundancia

//===============================================================================
 int i, j;
// reclutas anuales a la edad 2
  if(opt_Ro<0)
  {
   log_Ro=log_priorRo;
  }
 
  for (i=1;i<=nanos;i++)
  {
    N(i,1) = mfexp(log_Ro+log_desv_Rt(i)+0.5*square(sigmaR)); 
  }
    Neq(1) = mfexp(log_Ro+0.5*square(sigmaR));


// Abundancia inicial en equilibrio
  for (i=2;i<=nedades;i++)
  { 
    Neq(i) = Neq(i-1)*exp(-1*M);
  }
    Neq(nedades)=Neq(nedades)/(1-exp(-1*M));

    SSBo=sum(elem_prod(Neq*exp(-dt(4)*M),elem_prod(msex,colsum(Win)/nanos)));
    
// Abundancia inicial
  for (i=2;i<=nedades;i++)
  {
  N(1)(i)=Neq(i)*exp(log_desv_No(i-1)+0.5*square(sigmaR));
  }
// se estima la sobrevivencia por edad(a+1) y año(t+1)
  for (i=2;i<=nanos;i++)
  {
      N(i)(2,nedades)=++elem_prod(N(i-1)(1,nedades-1),S(i-1)(1,nedades-1));
      N(i,nedades)+=N(i-1,nedades)*S(i-1,nedades); 
  }
//===============================================================================
  

FUNCTION Eval_biomasas

//===============================================================================
// matrices y vectores de abundancias derivadas
  NVcru    = elem_prod(elem_prod(N,mfexp(-dt(1)*Z)),Scru);// Crucero Reclas
  NVpel    = elem_prod(elem_prod(N,mfexp(-dt(2)*Z)),Scru_pela);// Pelaces
// corrección por error de asignación de la edad
  if(erredad==1)
  {
  NVcru    = NVcru*error_edad;
  NVpel    = NVpel*error_edad;
  }
  NMD      = elem_prod(elem_prod(N,mfexp(-dt(3)*Z)),outer_prod(Unos_anos,msex));// desovante y MPH
  NVflo    = elem_prod(elem_prod(N,mfexp(-dt(4)*Z)),Sel_f);// explotable
  Reclutas = column(N,1);
// vectores de biomasas derivadas
  BD       = rowsum(elem_prod(NMD,Win));      // Desovante
  BT       = rowsum(elem_prod(N,Win));        // Total inicios de año biol
  BMflo    = rowsum(elem_prod(NVflo,Win));    // Biomasa explotable
  BMpel    = rowsum(elem_prod(NVpel,Win));    // pelaces
  Bcru     = rowsum(elem_prod(NVcru,Wmed));   // Reclas, mitad año biol
//===============================================================================

FUNCTION Eval_SRR
    
  RO = mean(Reclutas);
  SO = Bo;
  if(SrType==1){
	  alpha = (SO/RO)*((1-h)/(4*h));
	  beta = (5*h-1)/(4.*h*RO);
  }
  else
  {
	  alpha = (RO/SO)*exp(log(5*h)/0.8);
	  beta = log(5*h)/(0.8*SO);
   }	

  

FUNCTION Eval_capturas_predichas

//===============================================================================
// matrices de capturas predichas por edad y año
  pred_Ctot    = (elem_prod(elem_div(Ftot,Z),elem_prod(1.-S,N)));
// corrección por error de asignación de la edad
  if(erredad==1)
  {
  pred_Ctot    = pred_Ctot*error_edad;
  }
// vectores de desembarques predichos por año
  Desemb_pred  = rowsum(elem_prod(pred_Ctot,Wmed));
// matrices de proporcion de capturas por edad y año
  pobs_f       = elem_div(Ctot,outer_prod(rowsum(Ctot+1e-10),Unos_edad));
  ppred_f      = elem_div(pred_Ctot,outer_prod(rowsum(pred_Ctot),Unos_edad));
// matrices de capturas predichas por talla y año
// RECLAS EN EDADES
  pobs_crua    = elem_div(Ccru_a,outer_prod(rowsum(Ccru_a+1e-10),Unos_edad));
  ppred_crua   = elem_div(NVcru,outer_prod(rowsum(NVcru),Unos_edad));
// PELACES EN EDADES
  pobs_pel     = elem_div(Ccru_pel,outer_prod(rowsum(Ccru_pel+1e-10),Unos_edad));
  ppred_pel    = elem_div(NVpel,outer_prod(rowsum(NVpel),Unos_edad));
// PELACES EN TALLAS
  pobs_crul    = elem_div(Ccru_l,outer_prod(rowsum(Ccru_l+1e-10),Unos_tallas));
  ppred_crul   = elem_div(NVpel*Prob_talla,outer_prod(rowsum(NVpel),Unos_tallas));
  
//===============================================================================

FUNCTION Eval_indices

//===============================================================================
 Reclas_pred   = exp(log_qrecl)*Bcru;
 Pelaces_pred  = exp(log_qpela)*BMpel;
 MPH_pred      = exp(log_qmph)*BD;

 qrecl         = exp(log_qrecl);
 qpela         = exp(log_qpela);
 
 
 //===============================================================================
 
FUNCTION Eval_PBR

  //=============================================================================== 
  if(opt_Ro<0)
  {
    log_Ro=log_priorRo;
  }
  
  // Frms proxy (60%SPR y otros) y xx%SPR de Fmediana histórica
  for(int i=1;i<=npbr;i++){
    
    Fspr  = Sel_f(nanos)*log_Fref(i);
    Zspr  = Fspr+M;
    
    //Fmedian = Fmedian_ext;
    Fmedian = exp(mean(log_Ft));
    Fmed    = Sel_f(nanos)*Fmedian;
    Zmed    = Fmed+M;
    
    Nspro(1)=mfexp(log_Ro+0.5*square(sigmaR)); 
    Nspr(1)=mfexp(log_Ro+0.5*square(sigmaR)); 
    Nmed(1)=mfexp(log_Ro+0.5*square(sigmaR)); 
    
    for (int j=2;j<=nedades;j++)
    { 
      Nspro(j) = Nspro(j-1)*exp(-1*M);
      Nspr(j)  = Nspr(j-1)*exp(-Zspr(j-1));
      Nmed(j)  = Nmed(j-1)*exp(-Zmed(j-1));
    }
    Nspro(nedades) = Nspro(nedades)/(1-exp(-1*M));
    Nspr(nedades)  = Nspr(nedades)/(1-exp(-Zspr(nedades)));
    Nmed(nedades)  = Nmed(nedades)/(1-exp(-Zmed(nedades)));
    
    Bspro   = sum(elem_prod(Nspro*exp(-dt(3)*M),elem_prod(msex,colsum(Win)/nanos)));
    Bspr    = sum(elem_prod(elem_prod(elem_prod(Nspr,mfexp(-dt(3)*Zspr)),msex),colsum(Win)/nanos));
    Bsprmed = sum(elem_prod(elem_prod(elem_prod(Nmed,mfexp(-dt(3)*Zmed)),msex),colsum(Win)/nanos));
    
    ratio_spr(i) = Bspr/Bspro;	
    ratio_Fmed   = Bsprmed/Bspro;// xx%SPR de Fmediana
    
    // Bo y Brms proxy  según metodología Taller PBRs 2014
    //Bmed    = mean(BD(1,12));
    Bmed=mean(BD(1,nanos));
    Bo      = Bmed/(ratio_Fmed-0.05);
    Brms(i) = Bo*(ratio_spr(i)-0.05);
  }    	
  
 	
 	 
//===============================================================================

FUNCTION Eval_Estatus

//===============================================================================
  SSB=BD(1,nanos);   // variables de interés para mcmc 

// Rutina para calcular RPR dinámico
  Nv    = N;// solo para empezar los calculos
  
 for (int i=2;i<=nanos;i++)
  {
      Nv(i)(2,nedades)=++Nv(i-1)(1,nedades-1)*exp(-1.0*M);
      Nv(i)(nedades)+=Nv(i-1)(nedades)*exp(-1.0*M);
  }

  NDv  = elem_prod(Nv*exp(-dt(3)*M),outer_prod(Unos_anos,msex));
  BDo  = rowsum(elem_prod(NDv,Win));
  
  // INDICADORES DE REDUCCIÓN DEL STOCK
  RPRdin =  elem_div(BD,BDo);                       // RPR BDspr_t, dinámico
  RPRequ =  BD/Bspro;                               // RPR con BDspro
  RPRequ2 = BD/Bo;                                 // RPR con Bo proxy
  RPRequ3 = BD/Brms(1);                            // Razón para diagrama de fase
  Frpr    = exp(log_Ft)/log_Fref(1);
//===============================================================================

FUNCTION Eval_logverosim

//===============================================================================
// esta funcion evalua el nucleo de las -log-verosimilitudes marginales para
// series con datos 0.
  int i;
  suma1=0;
  suma2=0;
  suma3=0;
  suma4=0;

  for (i=1;i<=nanos;i++)
  {
      if (Reclas(i)>0)
      {
         suma1   += square((log(Reclas(i))-log(Reclas_pred(i)))/cvar(1,i));
      }
      if (Pelaces(i)>0)
      {
         suma2   += square((log(Pelaces(i))-log(Pelaces_pred(i)))/cvar(2,i));
      }
      if (MPH(i)>0)
      {
         suma3   += square((log(MPH(i))-log(MPH_pred(i)))/cvar(3,i));
      }
      if (Desemb(i)>0)
      {
         suma4   += square((log(Desemb(i))-log(Desemb_pred(i)))/cvar(4,i));}
  }
//===============================================================================

FUNCTION Eval_funcion_objetivo

//===============================================================================
// se calcula la F.O. como la suma de las -logver
// lognormal
  likeval(1)   = 0.5*suma1;//Reclas
  likeval(2)   = 0.5*suma2;//pelaces
  likeval(3)   = 0.5*suma4;//Desemb
  likeval(4)   = 0.5*suma3;//MPH
// multinomial
  likeval(5)   = -nmus(1)*sum(elem_prod(pobs_f,log(ppred_f)));
  likeval(6)   = -nmus(2)*sum(rowsum(elem_prod(pobs_crua,log(ppred_crua))));
  likeval(7)   = -nmus(3)*sum(rowsum(elem_prod(pobs_pel,log(ppred_pel))));
// tallas del pelaces
  likeval(8)   = -nmus(4)*sum(rowsum(elem_prod(pobs_crul,log(ppred_crul))));
//  Reclutas
  likeval(9)   = 1./(2*square(sigmaR))*norm2(log_desv_Rt);
// q cruceros
  likeval(10)  = 1./(2*square(cvpriorq_reclas))*square(log_qrecl);
  likeval(11)  = 1./(2*square(cvpriorq_pelaces))*square(log_qpela);

// Penaliza F 1991-1992 al promedio
  //likeval(12) = 1000*(square(log_Ft(2)-mean(log_Ft))+square(log_Ft(3)-mean(log_Ft)));  // S12

// Penaliza Fspr
  if(active(log_Fref)){
  likeval(12) = 1000*norm2(ratio_spr-ratio);}
 
// penalizacion de selectividades no-parametricas

  if (opt2_fase>0)
  {
    for (int i=1;i<=nedades-2;i++)
    {
       pSf(i) = lam(1)*sum(square(log(Sel_f(i))-2*log(Sel_f(i+1))+log(Sel_f(i+2))));
       pSc(i) = lam(2)*sum(square(log(Scru(i))-2*log(Scru(i+1))+log(Scru(i+2))));
    }
// penalizacion de las selectividades de la flota maximo en 1
    for (int i=1;i<=nanos;i++)
    {
       pStotf(i)=lam(3)*sum(square(log(Sel_f(i))));
    }
// penalizacion de las selectividades del crucero maximo en 1
    for (int i=1;i<=nanos;i++)
    {
    pStotc(i)=lam(4)*sum(square(log(Scru(i))));
    }
  }

// total
   if (opt1_fase>0)
   {
   f = sum(likeval);
   }

   if (opt2_fase>0)
   {
   f = sum(likeval)+sum(pSf)+sum(pSc)+sum(pStotf)+sum(pStotc);
   }
   if(mceval_phase())
   {
   Eval_mcmc();
   }
//==========================================================================================
   /*
 FUNCTION  Eval_CTP

//================================================================================================

//************************************************************************************************
// calcula la CTP para el ultimo año dado los PBR entregados CBA actual
//************************************************************************************************

  // Estimación de CBA AÑO BIOLÓGICO
  //for (int i=1;i<=npbr;i++) // ciclo de PBR
  //{
   // Fpbr   = Sel_f(nanos)*log_Fref(1)*1;
    //Zpbr   = Fpbr+M;
    //CTP    = elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),N(nanos)));
    //YTP(i) = sum(elem_prod(CTP,Wmed(nanos)));                              
  //}
  
  //************************************************************************************************
// calcula la CTP para el ultimo año dado los PBR entregados
// 1era y 2da revisión de CBA
// La estrategia de explotación o regla de control es igual a Frms constante
//************************************************************************************************

  // Estimación de CBA AÑO BIOLÓGICO sin proyectar//revisar último año!!!
    Fref_r0 = exp(log_Ft(nanos)); //log_Fref(1);//aquí debería ir F del último año
	Frms_r0 = Sel_f(nanos)*Fref_r0;
    Zrms_r0 = Frms_r0+M;
    CTP_r0  = elem_prod(elem_div(Frms_r0,Zrms_r0),elem_prod(1.-exp(-1.*Zrms_r0),N(nanos)));
    YTP_r0  = sum(elem_prod(CTP_r0,Wmed(nanos)));      
	NMD_r0  = elem_prod(elem_prod(N(nanos),mfexp(-dt(3)*Zrms_r0)),msex);
    BD_r0   = sum(elem_prod(NMD_r0,Win(nanos)));	
    RPR_r0  = BD_r0/Brms(1);
    

//================================================================
//----------------------------------------------------------------
//                   PROYECCIÓN DEL STOCK
//----------------------------------------------------------------
//================================================================
  if(opt_Ro<0){ log_Ro=log_priorRo; }//Esto corre cuando se hace perfil de verosimilitud
// Variables correspondientes al último año de evaluación
    Np     = N(nanos);
    Sp     = S(nanos);
    Nvp    = Nv(nanos);
    RPRp(1)= RPRequ3(nanos);
// Proyección del stock a partir del último año de evaluación
   for (int j=1;j<=nproy;j++){ // ciclo de 5 años
    Np(2,nedades)=++elem_prod(Np(1,nedades-1),Sp(1,nedades-1));
    Np(nedades)+=Np(nedades)*Sp(nedades);
//--------------------------------------------------------
 // Escenarios de reclutamiento promedio
    if(oprec==1){Np(1)=mean(Reclutas(1,nanos-11));} // reclutamiento promedio 1997-2009 (primeros 13 años)
    if(oprec==2){Np(1)=mean(Reclutas);} // reclutamiento promedio histórico (1997-año más reciente)
    if(oprec==3){Np(1)=mean(Reclutas(nanos-10,nanos));}  // reclutamiento promedio últimos años a partir del 2010 (2010-año más reciente)
    Npp = elem_prod(prop_est,Np);
    Wmedp=Wmedp_3;
	Winp=Winip_3;
//**************************************************************
// Proyección (p) con Regla de decisión N°0-Fconstante=Frms (0)
//**************************************************************
   Fref_p0 = mF*log_Fref(1);
   Frms_p0 = Sel_f(nanos)*Fref_p0;
   Zrms_p0 = Frms_p0+M;
  
   CTP_p0    = elem_prod(elem_div(Frms_p0,Zrms_p0),elem_prod(1.-exp(-1.*Zrms_p0),Npp)); 
   YTP_p0(j) = sum(elem_prod(CTP_p0,Wmedp)); 
   BD_p0(j)  = sum(elem_prod(elem_prod(elem_prod(Npp,mfexp(-dt(3)*Zrms_p0)),msex),Winp)); 
   RPR_p0(j) = BD_p0(j)/Brms(1);
		   
   //Nap(j)   = Npp;    
   Sp       = exp(-1.*Zrms_p0); 
//Proyección de cruceros
   NVrecl_p0   = elem_prod(elem_prod(Npp,mfexp(-dt(1)*Zrms_p0)),Scru(nanos));//considerar sólo mortalidad natural- Crucero Reclas!!!
   NVpel_p0    = elem_prod(elem_prod(Npp,mfexp(-dt(2)*Zrms_p0)),Scru_pela(nanos));
   Brecl_p0(j) = qrecl*sum(elem_prod(NVrecl_p0,Wmedp));
   Bpel_p0(j)  = qpela*sum(elem_prod(NVpel_p0,Winp)); 
	}		   
  //----------------------------------------------------------------
 // CÁLCULO DE CBA EN AÑO CALENDARIO
 // opciones según hito de estimación
 //----------------------------------------------------------------
 // regla Fconstante=Frms (0 = Fconst, 1 = regla mixta, r = mismo año, p=proyectado)
    if(opProy==1){CBA_c0=prop(1)*YTP_r0+prop(2)*YTP_p0(1);}     // Opción 1: 1era y 2da revisión (para el mismo año)
    if(opProy==2){CBA_c0=prop(1)*YTP_p0(1)+prop(2)*YTP_p0(2); } // Opción 2: CBA inicial (proyección de un año calendario)
       
  */
//##############################################################################

REPORT_SECTION

//##############################################################################
//----------------------------------------------------------------------------------
// VEROSIMILITUDES
//----------------------------------------------------------------------------------
// REPORTA SRR
	rep(SrType)
	rep(h)
	rep(SO)
	rep(RO)
	rep(alpha)
	rep(beta)
	  
//----------------------------------------------------------------------------------
// INDICES DE ABUNDANCIA
//----------------------------------------------------------------------------------
  rep(anos);
  rep(Reclas);
  rep(Reclas_pred);
  rep(Pelaces);
  rep(Pelaces_pred);
  rep(MPH);
  rep(MPH_pred);
  rep(Desemb);
  rep(Desemb_pred);
//-------------------------------------
// INDICADORES POBLACIONALES
//-------------------------------------
  rep(Reclutas);
  rep(log_desv_Rt);
  rep(BD);
  rep(BT);
  //Ftot;
  rep(exp(log_Ft));
//-------------------------------------
// INDICES DE REDUCCION
//-------------------------------------
  //BD_Bspro_t;
  rep(RPRdin);
  //BD_Bspro;
  rep(RPRequ);
  //BD_Bo;
  rep(RPRequ2);
  //BD_Brms;
  rep(RPRequ3);
//-------------------------------------
// SELECTIVIDADES
//-------------------------------------
  rep(Sel_f);
  rep(Scru);
  rep(Scru_pela);
//-------------------------------------
// PROPORCIÓN DE LAS CAPTURAS
//-------------------------------------
  rep(pobs_f);
  rep(ppred_f);
  rep(pobs_crua);
  rep(ppred_crua);
  rep(pobs_pel);
  rep(ppred_pel);
  rep(pobs_crul);
  rep(ppred_crul);
  
  //----------------------------------------
  // PUNTOS BIOLÓGICOS DE REFERENCIA TALLER
  //----------------------------------------
  //report << "pSPR Fmed_Fpbrs"<<endl;
  rep(ratio_Fmed);
  
  //report << "Fs Fmed_Fpbrs"<<endl;
  rep(Fmedian);
  rep(log_Fref);
  
  //report << "SSBpbr Bo_Bmed_Bpbrs"<<endl;
  rep(Bo); 
  rep(Bmed); 
  rep(Brms);
  
  //report << "Ro"<<endlNspro(1);
  rep(Nspro);
  
  //report << "SPR SPRFo_SPRFmed_SPRFpbrs"<<endl;
  rep(Bspro);
  rep(Bsprmed);
  rep(Bspr);
  

//----------------------------------------------------------------------------------
  //report << "log_Ro" << endl;
  rep(log_Ro);
  //report << "likeval ReclasPelacesDesembMPH_pf_preclas_ppelaces_ptallas_desvR_qrecl_qpela" << endl;
  rep(likeval);
//----------------------------------------------------------------------------------
// CAPTURABILIDADES
//----------------------------------------------------------------------------------
  //report << "q_reclas_q_pelaces" << endl;
  rep(exp(log_qrecl));
  rep(exp(log_qpela));
  //report << "M"<<endl;
  rep(M);
  //report << "N" << endl;
  rep(N);
  //report << "pred_Ctot" << endl;
  rep(pred_Ctot);
  //report << "F" << endl;
  rep(Ftot);

  
//##############################################################################  
//------------------------------------------------------------------------------

// ESTIMA nm y CV

//------------------------------------------------------------------------------
//############################################################################## 
  /*
 suma1=0; suma2=0;nm1=1;cuenta1=0;cuenta2=0;

  for (int i=1;i<=nanos;i++){ //

   if (sum(pobs_f(i))>0){
      suma1=sum(elem_prod(ppred_f(i),1-ppred_f(i)));
      suma2=norm2(pobs_f(i)-ppred_f(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}

  suma1=0; suma2=0;nm2=1;cuenta2=0;

  for (int i=1;i<=nanos;i++){ //

   if (sum(pobs_crua(i))>0){
      suma1=sum(elem_prod(ppred_crua(i),1-ppred_crua(i)));
      suma2=norm2(pobs_crua(i)-ppred_crua(i));
      nm2=nm2*suma1/suma2;
      cuenta2+=1;
   }}


  suma1=0; suma2=0;nm3=1;cuenta3=0;

  for (int i=1;i<=nanos;i++){ //

   if (sum(pobs_pel(i))>0){
      suma1=sum(elem_prod(ppred_pel(i),1-ppred_pel(i)));
      suma2=norm2(pobs_pel(i)-ppred_pel(i));
      nm3=nm3*suma1/suma2;
      cuenta3+=1;
   }}


  suma1=0; suma2=0;nm4=1; cuenta1=0;cuenta4=0;

  for (int i=1;i<=nanos;i++){ //

   if (sum(pobs_crul(i))>0){
      suma1=sum(elem_prod(ppred_crul(i),1-ppred_crul(i)));
      suma2=norm2(pobs_crul(i)-ppred_crul(i));
      nm4=nm4*suma1/suma2;
      cuenta4+=1;
   }}
  
  report << "nmprior"<<endl;
  report <<nmus<<endl;
  report << "nm_flota  nm_reclas  nm_pelaces    nm_pelaL" << endl;
  report<<pow(nm1,1/cuenta1)<<" "<<pow(nm2,1/cuenta2)<<" "<<pow(nm3,1/cuenta3)<<" "<<pow(nm3,1/cuenta4)<<endl;

  
  suma1=0;  suma2=0;  suma3=0;   suma4=0;  cuenta1=0;     cuenta2=0;   cuenta3=0;

  for (int i=1;i<=nanos;i++)
  {
   if (Reclas(i)>0){
    suma1+=square(log(Reclas(i))-log(Reclas_pred(i)));
    cuenta1+=1;}
   if (Pelaces(i)>0){
    suma2+=square(log(Pelaces(i))-log(Pelaces_pred(i)));
    cuenta2+=1;}
   if (MPH(i)>0){
    suma4+=square(log(MPH(i))-log(MPH_pred(i)));
   cuenta4+=1;}
  }


 report << "cv_recla  cv_pelaces  cv_mph" << endl;
 report<<sqrt(suma1/cuenta1)<<" "<<sqrt(suma2/cuenta2)<<" "<<sqrt(suma4/cuenta4)<<endl;
 

  if(erredad==1){
  report << " ------------------------------------------------" << endl;
  report << "Matriz de error "<< endl;
  report << error_edad << endl;}

  report << " ------------------------------------------------" << endl;
  report << "Talla a la edad & desviación "<< endl;
  report << mu_edad << endl;
  report << sigma_edad << endl;

   */
FUNCTION Eval_mcmc

  if(reporte_mcmc == 0)
  mcmc_report<<"f, RPR, RBV, F_last, CTP_60, CTP_obj, Recl_last"<<endl;
  //mcmc_report<<f<<","<<RPRdin(nanos)<<","<<RPRequ(nanos)<<","<<max(Ftot(nanos))<<"," <<YTP(1)<<","<<YTP(2)<<","<<Reclutas(nanos)<<endl;
  reporte_mcmc++;

  



 
