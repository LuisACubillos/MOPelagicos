# MOPelagicos

Evaluación de Estrategias de Manejo para las pesquerías de peces pelágicos en Chile: a) anchoveta (Atacama-Coquimbo), b) anchoveta (Valparaíso-Los Lagos), y c) sardina común (Valparaíso-Los Lagos), y d) sardina austral (Los Lagos).

La codificación original (TPL y DAT) ha sido documentada por investigadores del [Instituto de Fomento Pesquero](https://www.ifop.cl/busqueda-de-informes/) y en las actas e informes del 2020 al 2021 del [Comité Científico de Pesquerías de Pequeños Pelágicos](https://www.subpesca.cl/portal/616/w3-propertyvalue-51142.html#collapse03) de la [Subsecretaría de Pesca y Acuicultura](https://www.subpesca.cl/portal/616/w3-channel.html).

Los modelos originales fueron transformados en modelos operativos por L. Cubillos y M.J. Cuevas a objeto de evaluar el desempeño de los estimadores respecto de la aplicación de la mortalidad por pesca objetivo (PBR) dada la dinámica del reclutamiento que constituye la principal incertdidumbre estructural en estas pesquerías.

Esta actividad es parte de los resultados del Objetivo 2 del proyecto FIPA 2019-17 "Asesoría para la revisión de PBRs y consideraciones ecosistémicas asociados a pesquerías pelágicas.", liderado por Sergio Neira, Centro COPAS Sur-Austral, Departamento de Oceanografía, Universidad de Concepción, Concepción, Chile.

Encargados de la evaluación monoespecífica:

* María José Cuevas: Codificación, evaluación y análisis de resultados.
* Luis A. Cubillos: Gestor del proyecto, conceptualización, evaluación de estrategias de manejo, codificación, evaluación de desempeño y análisis de resultados.


## Nota

Los códigos están escritos en [ADMB](http://www.admb-project.org/) templates ([Fournier et al. 2012](https://doi.org/10.1080/10556788.2011.597854)), y alojados en las siguientes carpetas:

*Anzcn2020:* Contiene el modelo de evaluación de stock o estimador (Est1), donde se identifican tres archivos: MATT2009.tpl (plantilla de código ADMB), MATT2009_real.data (respaldo de los datos utilizados en la evaluación), y control_est.ctl que contiene el control de la estimación.

Además, están los modelos operativos que consideran las siguientes configuraciones:

OP1: La dinámica del reclutamiento considera una relación stock-recluta del tipo Beverton-Holt con variabilidad estocástica ~N(E(R),SigmaR)

OP2: La dinámica del reclutamiento considera una relación stock-recluta del tipo Ricker con variabilidad estocástica ~N(E(R),SigmaR)

OP3: La dinámica del reclutamiento considera una relación stock-recluta del tipo Beverton-Holt con autocorrelación en los desvios del reclutamiento.





### Requerimientos

* Un compilador C++
* ADMB-12.0 o superior

### Clonación

	cd ~
	git clone https://github.com/luisacubillos/MOPelagicos
	cd MOPelagicos

### Instucciones

1) Compile los archivos *.tpl

2) Revise las rutas contenidas en los archivos "paso1.sh" y "paso2.sh" que contienen códigos "Shell Script (Bach)" que son ejecutados por el modelo operativo para comunicarse con el estimador. Estos archivos son para Mac, pero pueden sis se utiliza sistema operativo Windows deben ser transformados a archivos *.bat.

3) Las simulaciones se basan en muestrear un número dado de valores alternativos e igualmente probables desde la estiamción posterior del modelo operativo mediante MCMC; por ejemplo: "./OP1anzcn -mcmc 10000 -mcsave 50", generará 200 valores alternativos que determinarán la trayectoria futura de la dinámica poblacional sujeta a la simulación de datos (OP1) y al impacto de una captura biológicamente aceptable que calcula cada año de proyección el estimador (Est1).

4) El proceso se activa desde el modelo operativo utilizando mceval; por ejemplo: "./OP1anzcn - mceval". Esto activará la comunicación entre estimador y operativo.

5) Para realizar una estiamción simultánea de diferentes OP, se sugiere copiar Est1 en carpetas como Est2, Est 3. Aunque esto no es necesario, dependerá de las capacidades del computador que se utilice. Se destaca que un proceso completo puede tardar 1 a 2 días de trabajo en el computador con 32 G de RAM

