# MOPelagicos

Evaluación de Estrategias de Manejo para las pesquerías de peces pelágicos en Chile: a) anchoveta (Atacama-Coquimbo), b) anchoveta (Valparaíso-Los Lagos), y c) sardina común (Valparaíso-Los Lagos), y d) sardina austral (Los Lagos).

La codificación original (TPL y DAT) ha sido documentada por investigadores del [Instituto de Fomento Pesquero](https://www.ifop.cl/busqueda-de-informes/) y en las actas e informes del 2020 al 2021 del [Comité Científico de Pesquerías de Pequeños Pelágicos](https://www.subpesca.cl/portal/616/w3-propertyvalue-51142.html#collapse03) de la [Subsecretaría de Pesca y Acuicultura](https://www.subpesca.cl/portal/616/w3-channel.html).

Los modelos originales fueron transformados en modelos operativos por L. Cubillos y M.J. Cuevas a objeto de evaluar el desempeño de los estimadores respecto de la aplicación de la mortalidad por pesca objetivo (PBR) dada la dinámica del reclutamiento que constituye la principal incertdidumbre estructural en estas pesquerías.

Esta actividad es parte de los resultados del Objetivo 2 del proyecto FIPA 2019-17 "Asesoría para la revisión de PBRs y consideraciones ecosistémicas asociados a pesquerías pelágicas.", liderado por Sergio Neira, Centro COPAS Sur-Austral, Departamento de Oceanografía, Universidad de Concepción, Concepción, Chile.

Encargados de la evaluación monoespecífica:

* María José Cuevas: Codificación, evaluación y análisis de resultados.
* Luis A. Cubillos: Gestor del proyecto, conceptualización, evaluación de estrategias de manejo, codificación, evaluación de desempeño y análisis de resultados.


## Aspectos técnicos

Los códigos están escritos en [ADMB](http://www.admb-project.org/) templates ([Fournier et al. 2012](https://doi.org/10.1080/10556788.2011.597854)), y alojados en carpetas para cada pesquería que contienen:

Est1: Modelo de evaluación de stock o estimador.

OP1, OP2, OP3, y OP4: Modelos operativos con diferentes configuraciones para la dinámica del reclutamiento; a saber:

OP1: La dinámica del reclutamiento considera una relación stock-recluta del tipo Beverton-Holt con variabilidad estocástica ~N(E(R),SigmaR)

OP2: La dinámica del reclutamiento considera una relación stock-recluta del tipo Ricker con variabilidad estocástica ~N(E(R),SigmaR)

OP3: La dinámica del reclutamiento considera una relación stock-recluta del tipo Beverton-Holt con autocorrelación en los desvios del reclutamiento e(t)=r*e(t-1)+sqrt(1-r^2 )*n(t), donde e(t) es el desvío de reclutamiento en el año t, r es la correlación de primer orden (fija en 0.8), y n(t) es un numero aleatorio desde una distribución normal con media cero y desviación estándar SigmaR.

OP4: Cambios de régimen en el reclutamiento, según dos estados determinados con la técnica Hidden Markov Model.

* Anchoveta (Atacama-Coquimbo)

*Anzcn2020:* Contiene el modelo de evaluación de stock o estimador (Est1), donde se identifican tres archivos: MATT2009.tpl (plantilla de código ADMB), MATT2009_real.data (respaldo de los datos utilizados en la evaluación), y control_est.ctl que contiene el control de la estimación. Note que MATT2009.dat contiene datos simulados por el modelo operativo.

Además, están los modelos operativos que consideran las siguientes configuraciones: OP1, OP2, OP3

* Anchoveta (Valparaíso-Los Lagos)

*Anzcs2020:* Contiene el modelo de evaluación de stock o estimador (Est1), donde se identifican tres archivos: MAE0920b.tpl (plantilla de código ADMB), MAE0920_real.data (respaldo de los datos utilizados en la evaluación), y control_est.ctl que contiene el control de la estimación. *Nota:* A modo de ejemplo, cuando el numero de simulaciones termina queda en carpeta el último archivo *.dat; en este caso MAE0920b.dat contiene los últimos datos simulados por el modelo operativo.

Además, están los modelos operativos que consideran las siguientes configuraciones: OP1, OP2, OP3, y OP4.

* Sardina común (Valparaíso-Los Lagos)

*Sczcs2020:* Contiene el modelo de evaluación de stock o estimador (Est1), donde se identifican tres archivos: MAE0920.tpl (plantilla de código ADMB), MAE0920_real.data (respaldo de los datos utilizados en la evaluación), y control_est.ctl que contiene el control de la estimación. *Nota:* A modo de ejemplo, cuando el numero de simulaciones termina queda en carpeta el último archivo *.dat; en este caso MAE0920.dat contiene los últimos datos simulados por el modelo operativo.

Además, están los modelos operativos que consideran las siguientes configuraciones: OP1, OP2, OP3, y OP4.

* Sardina austral (Los Lagos)

*Saus2020:* Contiene el modelo de evaluación de stock o estimador (Est1), donde se identifican tres archivos: SAMSau.tpl (plantilla de código ADMB), SAMsau_real.data (respaldo de los datos utilizados en la evaluación), y SAMsau_control.txt que contiene el control de la estimación.

Además, están los modelos operativos que consideran las siguientes configuraciones: OP1, OP2, OP3.

### Requerimientos

* Un compilador C++
* ADMB-12.0 o superior

### Clonación

	cd ~
	git clone https://github.com/luisacubillos/MOPelagicos
	cd MOPelagicos

### Instucciones

1) Compile los archivos *.tpl

2) Revise las rutas contenidas en los archivos "paso1.sh" y "paso2.sh" que contienen códigos "Shell Script (Bach)" que son ejecutados por el modelo operativo para comunicarse con el estimador. Estos archivos son para Mac, pero pueden sis se utiliza sistema operativo Windows deben ser transformados a archivos *.bat y utilizar el lenguaje de la terminal de windows.

3) Las simulaciones se basan en muestrear un número dado de valores alternativos e igualmente probables desde la estimación posterior del modelo operativo mediante MCMC; por ejemplo: "./OP1anzcn -mcmc 10000 -mcsave 50", generará 200 valores alternativos que determinarán la trayectoria futura de la dinámica poblacional sujeta a la simulación de datos (OP1) y al impacto de una captura biológicamente aceptable que calcula cada año de proyección el estimador (Est1).

4) El proceso se activa desde el modelo operativo utilizando mceval; por ejemplo: "./OP1anzcn - mceval". Esto activará la comunicación entre estimador y operativo.


5) Los resultados del procedimiento genera los siguientes archivos individuales:

00rep_convergencia.txt: reporte de convergencia del estimador (0: falta de convergencia; 1: convergencia).

01Capturas_proyectadas.txt: capturas biológicamente aceptables por el estimador.

02BiomasaTotal_est.txt: biomasa total del estimador
02BiomasaTotal_op.txt: biomasa total del modelo operativo

03Desovante_est.txt: biomasa desovante del estimador
03Desovante_op.txt: biomasa desovante del modelo operativo

04Reclutamiento_est.txt: reclutamiento del estimador. 
04Reclutamiento_op.txt: reclutamiento del modelo operativo.

05FMort_est.txt: mortalidad por pesca del estimador.
05FMort_op.txt: mortalidad por pesca impactada por la CBA en el modelo operativo.

08Desovante_hist.txt: Biomasa desovante en el periodo de evaluación del modelo operativo.
09BiomTotal.txt: Biomasa total en el periodo de evaluación del modelo operativo. 
10Reclutas_hist.txt : Reclutamiento en el periodo de evaluación del modelo operativo.
11RPR_hist.txt: razón desovante/biomasa inexplotada del modelo operativo en el periodo de evaluación.
12RPR_est.txt: razón desovante/biomasa inexplotada del estimador
12RPR_op.txt: razón desovante/biomasa inexplotada del modelo operativo

6) Para realizar una estiamción simultánea de diferentes OP, se sugiere copiar Est1 en carpetas como Est2, Est 3. Aunque esto no es necesario, dependerá de las capacidades del computador que se utilice. Se destaca que un proceso completo puede tardar 1 a 2 días de trabajo en el computador con 32 G de RAM

