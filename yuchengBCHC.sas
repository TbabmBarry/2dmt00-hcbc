libname SAS "/folders/myshortcuts/sharefolders";

data BCHC0;
set SAS.bchc_2dmt00;
run;

data BCHC;
set bchc0;
where Race_Ethnicity<>'All';
run;

proc contents data=bchc0;
run;

proc freq data=bchc;
table place*Race_Ethnicity; 
run;

data BCHC1;
set bchc0;
where Race_Ethnicity='All';
run;

/*check the content of the data*/
proc contents data=bchc;
run;

/*ANOVA model*/
proc mixed data=bchc0 method=type3;
class Race_Ethnicity place year;
model mortality=Race_Ethnicity/solution cl outp=pred;
random  place year/solution;
run;

proc mixed data=bchc method=type3;
class Race_Ethnicity place;
model mortality=Race_Ethnicity/solution cl outp=pred;
random place/solution;
run;
proc mixed data=bchc method=type3;
class  place;
model mortality=/solution cl outp=pred;
random place/solution;
run;

proc mixed data=bchc1 method=type3;
class Race_Ethnicity place;
model mortality=Race_Ethnicity/solution cl outp=pred;
random place/solution;
run;

/*test the normality of the residual*/
proc univariate data=pred normal;
var resid;
histogram resid/normal;
probplot resid/normal(MU=est SIGMA=est);
run;

/*check where the normality is violated*/
DATA approx;
N=187;
G1=0.99606809;
G2=3.32063222;
b1=(N-2)*G1/(sqrt(N*(N-1)));
b2=G2*((N-2)*(N-3))/((N+1)*(N-1))+3*(N-1)/(N+1);
*JB=N*(b1**2/6+(G2-3)**2/24);
Cn=(3*(N**2+27*N-70)*(N+1)*(N+3))/((N-2)*(N+5)*(N+7)*(N+9));
Wn2=-1+SQRT(2*(Cn-1));
Alphan=SQRT(2/(Wn2-1));
Dn=1/sqrt(log(sqrt(Wn2)));
Bn=sqrt((N+1)*(N+3)/(6*(N-2)))*b1;
Ts=Dn*log(Bn/Alphan+sqrt(1+(Bn/Alphan)**2));
Mun=3*(N-1)/(N+1);
Sigman=sqrt((24*N*(N-2)*(N-3))/((N+3)*(N+5)*(N+1)**2));
Gamma1n=((6*(N**2-5*N+2))/((N+7)*(N+9)))*sqrt(6*(N+3)*(N+5)/(N*(N-2)*(N-3)));
An=6+(8/(Gamma1n))*(2/Gamma1n+sqrt(1+4/(Gamma1n**2)));
Un=(b2-Mun)/Sigman;
Tk=sqrt(9*An/2)*((9*An-2)/(9*An)-((1-2/An)/(1+Un*sqrt(2/(An-4))))**(1/3));
K2=Tk**2+Ts**2;
Ps=2*min(cdf('Normal',Ts,0,1),1-cdf('Normal',Ts,0,1));
Pk=2*min(cdf('Normal',Tk,0,1),1-cdf('Normal',Tk,0,1));
PK2=1-cdf('chisq',K2,2);
*PJB=1-cdf('chisq',JB,2);
run;

proc print data=approx;
var Ts Tk K2 Ps Pk PK2 ;
run;

/*BOX-COX Transformation*/
DATA BOXCOX;
 	SET BCHC;
 	MORTALITYMINUS2= (-1/2)*((mortality)**-2-1);
 	MORTALITYMINUS1= (-1)*((mortality)**-1-1);
 	MORTALITYMINUS12= (-2)*((mortality)**-(0.5)-1);
 	MORTALITY0= log(mortality);
 	MORTALITY13= (3)*((mortality)**(1/3)-1);
 	MORTALITY12= (2)*((mortality)**(1/2)-1);
 	MORTALITY2= (0.5)*((mortality)**(2)-1);
RUN;

proc univariate data=BOXCOX noprint;
  	histogram MORTALITYMINUS2 /normal;
 	histogram MORTALITYMINUS1 /normal;
 	histogram MORTALITYMINUS12 /normal;
 	histogram MORTALITY0 /normal;
 	histogram MORTALITY13 /normal;
 	histogram MORTALITY12 /normal;
	histogram MORTALITY2 /normal;
run;

data BCHC_log;
set bchc;
logm= log(mortality);
/* drop _TYPE_ _Freq_; */
run;

proc univariate data=BCHC_log cipctldf;
	var logm;
run;

DATA TUKEY;
SET BCHC_log;
p75=3.45947;
p25=2.67415;
IQR=p75-p25;
LOWERT = p25 - 1.5*IQR;
UPPERT = p75 + 1.5*IQR;
RUN;

DATA TUKEY;
SET TUKEY;
T=(logm>UPPERT OR logm<LOWERT);
RUN;

PROC SORT data=TUKEY;
by descending T;
run;

PROC PRINT DATA=TUKEY;
RUN;





data BCHC_logAll;
set bchc1;
logm= log(mortality);
run;

proc univariate data=BCHC_logAll cipctldf;
	var logm;
run;

DATA TUKEY;
SET BCHC_logAll;
p75=3.13766;
p25=2.83005;
IQR=p75-p25;
LOWERT = p25 - 1.5*IQR;
UPPERT = p75 + 1.5*IQR;
RUN;

DATA TUKEY;
SET TUKEY;
T=(logm>UPPERT OR logm<LOWERT);
RUN;

PROC SORT data=TUKEY;
by descending T;
run;

PROC PRINT DATA=TUKEY;
RUN;

data BCHC_logAll_NoOuterlier;
	set tukey;
	where t=0;*delete outliers;
run;

proc univariate data=bchc_logall_noouterlier normal;
var logm;
histogram logm/normal;
probplot logm/normal(MU=est SIGMA=est);
run;


/*New ANOVA model*/
proc mixed data=bchc_log method=type3;
class Race_Ethnicity place;
model logm=Race_Ethnicity/solution cl outp=pred;
random place/solution;
lsmeans race_ethnicity/DIFF ADJUST=TUKEY CL;
run;

proc mixed data=BCHC_logAll_NoOuterlier method=type3;
class place;
model logm=/solution cl outp=pred;
random place/solution;
run;

/*test the normality*/
proc univariate data=pred normal;
var resid;
histogram resid/normal;
probplot resid/normal(MU=est SIGMA=est);
run;

/*test the homogeneity of variance of residuals*/
proc glm data=pred;
class race_ethnicity;
model resid= race_ethnicity;
means  race_ethnicity/hovtest=bf;
run;

proc glm data=pred;
class place;
model resid= place;
means  place/hovtest=bf;
run;

proc freq data=bchc_log;
table place*Race_Ethnicity; 
run;
/*MSE=0.056789*/
/* MSB*/
proc iml;
	sigA=0.06838;
	sigE=0.05679;
	sigB=9*0.056789/4;
	icc=sigA/(sigA+sigB+sigE);
	print(ICC);
run;


/* sigmaGALL= 0.04937 resid=0.006196;*/
/* sigmaRace 0.06838 resid=0.05679*/
/* sigma0.06067 0.1951 */
proc iml;
	sigAll=0.04937;
	sigRace=0.06838;
	sigE1=0.006196;
	sigE2=0.05679;
	icc1=sigAll/(sigE1+sigAll);
	icc2=sigRace/(sigE2+sigRace);	
	A=icc1||icc2||1-(icc2/icc1);
	create Indexes from A[colname={'icc1' 'icc2' 'index'}];
	append from A;
	close Indexes;
quit;
proc print data=Indexes;
run;
