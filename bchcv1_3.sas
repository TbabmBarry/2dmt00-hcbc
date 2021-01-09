/**
 *Raw Data Analysis
 */

/* Load BCHC dataset*/
data bchc;
	set "/folders/myshortcuts/sharefolders/bchc_2dmt00.sas7bdat";
run;

/* View basic infos of the dataset */
proc contents data=bchc;
run;

/* Check whether the data is balanced, if there is repeated measure*/
proc freq data=bchc;
	table place*Race_Ethnicity;
run;

/*Split dataset into 2 datasets*/
/* New dataset after excluding Race_Ethnicity=All */
data bchc_race;
	set bchc;
	id= _n_;
	where Race_Ethnicity <> 'All';
    rename Race_Ethnicity = race;
run;

/* New dataset that only include Race_Ethnicity=All */
data bchc_all;
	set bchc;
	id= _n_;
	where Race_Ethnicity = 'All';
    rename Race_Ethnicity = race;
run;

/*Test Normality of data*/
proc univariate data=bchc_race normal;
var mortality;
histogram mortality/normal;
probplot mortality/normal(MU=est SIGMA=est);
run;
/*check where the normality is violated*/
DATA approx;
N=187;
G1=1.44341126;
G2=2.75482979;
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

/*Test Normality of data*/
proc univariate data=bchc_all normal;
var mortality;
histogram mortality/normal;
probplot mortality/normal(MU=est SIGMA=est);
run;
/*check where the normality is violated*/
DATA approx;
N=48;
G1=3.92338088;
G2=21.8925919;
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
 	SET bchc_race;
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

/*log transformation*/
data bchc_race_log;
set bchc_race;
logm= log(mortality);
run;

/*check normality*/
proc univariate data=bchc_race_log normal;
var logm;
histogram logm/normal;
probplot logm/normal(MU=est SIGMA=est);
run;




/*BOX-COX Transformation*/
DATA BOXCOX;
 	SET bchc_all;
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


data bchc_all_log;
set bchc_all;
logm= log(mortality);
run;

proc univariate data=bchc_all_log normal;
var logm;
histogram logm/normal;
probplot logm/normal(MU=est SIGMA=est);
run;

/*check where the normality is violated*/
DATA approx;
N=48;
G1=0.89745605;
G2=4.60108033;
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

/*Remove Outliers*/
proc univariate data=bchc_all_log cipctldf;
	var logm;
run;

DATA TUKEY;
SET bchc_all_log;
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

data bchc_all_log;
	set tukey;
	where t=0;*delete outliers;
run;

****************************************************************************;



























/**
 * Model Selection
 */

/*Model selction for dataset bchc_race_log*/
/*one-way ANOVA model */
proc mixed data=bchc_race_log method=type3 cl;
    class race;
    model logm = race / solution outp=pred;
run;

proc mixed data=bchc_race_log method=type3 cl;
    class place;
    model logm = / solution outp=pred;
    random place;
run;

proc mixed data=bchc_race_log method=type3 cl;
    class year;
    model logm = / solution outp=pred;
    random year;
run;

/*two-way ANOVA model*/
proc mixed data=bchc_race_log method=type3 cl;
    class place race;
    model logm = race / solution outp=pred_race;
    random place;
run;

proc mixed data=bchc_race_log method=type3 cl;
    class race year;
    model logm = race / solution outp=pred_race;
    random year;
run;

proc mixed data=bchc_race_log method=type3 cl;
    class place year;
    model logm = / solution outp=pred_race;
    random place year;
run;

/*Three-way ANOVA model*/
proc mixed data=bchc_race_log method=type3 cl;
    class place race year;
    model logm = race / solution outp=pred_race;
    random place year;
run;

/*Sum Squares of year is small*/

/*Final model after selection*/
proc mixed data=bchc_race_log method=type3 cl;
    class place race;
    model logm = race / solution outp=pred_race;
    random place;
run;



/*Model selction for dataset bchc_all_log*/
proc mixed data=bchc_all_log method=type3 cl;
    class place;
    model logm = / solution outp=pred_all;
    random place;
run;

******************************************************************;


































/**
 * Verify Model's Assumption
 */
/*Test the assumption of model where dataset is bchc_all_log*/

/* check if our residual follows a normal distribution */
proc univariate data=pred_all normal;
    var resid;
    histogram resid;
    probplot resid /normal(mu=est sigma=est);
run;

/* Residual homogeneity of variance test again */
proc glm data=pred_race;
    class race;
    model resid = race;
    means race /hovtest=Bartlett;
run;
/*Randomness Test*/

%MACRO RUNSCUC(data=,var=,alpha=);
PROC IML;
use &data;
read all var {&var};
close &data;

X=&var;
n=nROW(X);
MED=median(X);

XC=X;
DO i=1 to n by 1;
	IF (XC[i] >= MED) then XC[i]=1;
	ELSE XC[i]=0;
END;

n1C=sum(XC);
n2C=n-n1C;

RC=1;
DO i=2 to n by 1;
	if(XC[i] ^= XC[i-1]) then RC=RC+1;
END;

MUC=1+(2*n1C*n2C)/(n1C+n2C);
VARC=2*n1C*n2C*(2*n1C*n2C-n1C-n2C)/((n1C+n2C-1)*(n1C+n2C)**2);

SC=(RC-MUC)/SQRT(VARC);
TC=QUANTILE('NORMAL',&alpha/2);
TCU=QUANTILE('NORMAL',1-&alpha/2);
PC=(1-CDF('NORMAL',abs(SC)))*2;

XUC=REPEAT(0,n-1,1);
TIES=0;
DO i=1 to (n-1) by 1;
	IF (X[i+1] > X[i]) then XUC[i]=1;
	IF (X[i+1] < X[i]) then XUC[i]=0;
	IF (X[i+1] = X[i]) then XUC[i]=XUC[i-1];
	IF (X[i+1] = X[i]) then TIES=TIES+1;
END;

RUC=1;
DO i=2 to (n-1) by 1;
	if(XUC[i] ^= XUC[i-1]) then RUC=RUC+1;
END;

MUUC=(2*(n-TIES)-1)/3;
VARUC=(16*(n-TIES)-29)/90;

SUC=(RUC-MUUC)/SQRT(VARUC);
TUC=QUANTILE('NORMAL',&alpha/2);
TUCU=QUANTILE('NORMAL',1-&alpha/2);
PUC=(1-CDF('NORMAL',abs(SUC)))*2;

PRINT("Median based (conditional) runs test");
PRINT(RC||MUC||sqrt(VARC)||PC||SC||TC||TCU||n);
PRINT("(unconditional) runst test for serial randomness");
PRINT(TIES);
PRINT(RUC||MUUC||sqrt(VARUC)||PUC||SUC||TUC||TUCU||(n-TIES));
quit;
%MEND;

data rm;
    set rm;
    residm = resid;
    predm = pred;
    drop resid pred;
run;
data predmerge;
    merge rm rc;
    by id place;
run;
data predmerge;
    set predmerge;
    eblup = pred-predm;
run;
proc sort data=predmerge;
    by eblup;
run;

proc print data=predmerge;run;

proc sgplot data=predmerge;
    scatter X=eblup Y=residm;
run;

%RUNSCUC(data=predmerge,var=residm,alpha=0.05);

proc sgplot data=predmerge;
    scatter x = eblup y = resid;
run;

%RUNSCUC(data=predmerge,var=resid,alpha=0.05);




/*Test the assumption of model where dataset is bchc_race_log*/
/* check if our residual follows a normal distribution */
proc univariate data=pred_race normal;
    var resid;
    histogram resid;
    probplot resid /normal(mu=est sigma=est);
run;
/* We don't have to check residual homogeneity of variance here*/
/* We don’t have to test for homogeneity of variance among the groups defined by random effects */

/*Randomness Test*/
proc mixed data=bchc_all_log method=type3 cl;
    class place;
    model logm = / solution outpm=rm outp=rc;
    random place;
run;

%MACRO RUNSCUC(data=,var=,alpha=);
PROC IML;
use &data;
read all var {&var};
close &data;

X=&var;
n=nROW(X);
MED=median(X);

XC=X;
DO i=1 to n by 1;
	IF (XC[i] >= MED) then XC[i]=1;
	ELSE XC[i]=0;
END;

n1C=sum(XC);
n2C=n-n1C;

RC=1;
DO i=2 to n by 1;
	if(XC[i] ^= XC[i-1]) then RC=RC+1;
END;

MUC=1+(2*n1C*n2C)/(n1C+n2C);
VARC=2*n1C*n2C*(2*n1C*n2C-n1C-n2C)/((n1C+n2C-1)*(n1C+n2C)**2);

SC=(RC-MUC)/SQRT(VARC);
TC=QUANTILE('NORMAL',&alpha/2);
TCU=QUANTILE('NORMAL',1-&alpha/2);
PC=(1-CDF('NORMAL',abs(SC)))*2;

XUC=REPEAT(0,n-1,1);
TIES=0;
DO i=1 to (n-1) by 1;
	IF (X[i+1] > X[i]) then XUC[i]=1;
	IF (X[i+1] < X[i]) then XUC[i]=0;
	IF (X[i+1] = X[i]) then XUC[i]=XUC[i-1];
	IF (X[i+1] = X[i]) then TIES=TIES+1;
END;

RUC=1;
DO i=2 to (n-1) by 1;
	if(XUC[i] ^= XUC[i-1]) then RUC=RUC+1;
END;

MUUC=(2*(n-TIES)-1)/3;
VARUC=(16*(n-TIES)-29)/90;

SUC=(RUC-MUUC)/SQRT(VARUC);
TUC=QUANTILE('NORMAL',&alpha/2);
TUCU=QUANTILE('NORMAL',1-&alpha/2);
PUC=(1-CDF('NORMAL',abs(SUC)))*2;

PRINT("Median based (conditional) runs test");
PRINT(RC||MUC||sqrt(VARC)||PC||SC||TC||TCU||n);
PRINT("(unconditional) runst test for serial randomness");
PRINT(TIES);
PRINT(RUC||MUUC||sqrt(VARUC)||PUC||SUC||TUC||TUCU||(n-TIES));
quit;
%MEND;

data rm;
    set rm;
    residm = resid;
    predm = pred;
    drop resid pred;
run;
data predmerge;
    merge rm rc;
    by id place;
run;
data predmerge;
    set predmerge;
    eblup = pred-predm;
run;
proc sort data=predmerge;
    by eblup;
run;

proc print data=predmerge;run;

proc sgplot data=predmerge;
    scatter X=eblup Y=residm;
run;

%RUNSCUC(data=predmerge,var=residm,alpha=0.05);

proc sgplot data=predmerge;
    scatter x = eblup y = resid;
run;

%RUNSCUC(data=predmerge,var=resid,alpha=0.05);

*******************************************************************;






























/**
 * Problem Solving
 */

/*Question1*/

/*Should we use Satterthwaite approach?*/
proc mixed data=bchc_race_log method=type3 cl;
    class place race;
    model logm = race / solution outp=pred_race;
    random place;
    lsmeans race/DIFF ADJUST=TUKEY CL;
run;
/*The back transformation is not done yet*/

/*Question2*/
proc mixed data=bchc_all_log method=type3 cl;
    class place;
    model logm = / solution outpm=rm outp=rc;
    random place;
run;

/*One way ANOVA sigmaG= 0.04937 resid=0.006196;*/
/*Two way ANOVA sigmaG= 0.06838 resid=0.05679*/
/*icc1 is the one way ANOVA model's intraclass correlation coefficient*/
/*icc2 is the two way ANOVA model's intraclass correlation coefficient*/
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

/*Question3*/
/*By choosing 3 cities whose absolute value of EBLUP is significantly larger than others*/
