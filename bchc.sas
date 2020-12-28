
/* Load BCHC dataset */
data bchc;
	set "/folders/myfolders/bchc_2dmt00.sas7bdat";
run;

/* View basic infos of the dataset */
proc contents data=bchc;
run;

/* Check if there are repeated observations */
proc freq data=bchc;
	table place*Race_Ethnicity;
run;

/* Create formats to convert Year, Place, and Race_Ethnicity into categorical numeric value */
proc format;
	invalue numPlace
		'Chicago, Il' = 1
		'Las Vegas (Clark County), NV' = 2
		'Long Beach, CA' = 3
		'New York City, NY' = 4
		'Oakland (Alameda County), CA' = 5
		'Phoenix, AZ' = 6
		'San Diego County, CA' = 7
		'San Francisco, CA' = 8
		'San Jose, CA' = 9
		'Seattle, WA' = 10
		;
run;

proc format;
	invalue numYear
		2013 = 1
		2014 = 2
		2015 = 3
		2016 = 4
		2017 = 5
		2018 = 6
		2019 = 7
		;
run;
		

proc format;
	invalue numRace
		'All' = 0
		'Asian/PI' = 1
		'Black' = 2
		'Hispanic' = 3
		'White' = 4
	;
run;

/* New dataset excluding Race_Ethnicity=All */
data bchc_race;
	set bchc;
    id_new = _n_;
	*race=input(Race_Ethnicity,numRace.);
	*pl=input(Place,numPlace.);
	*yr=input(cats(Year),numYear.);
	where Race_Ethnicity <> 'All';
    rename Race_Ethnicity = race;
run;

/* Check potential interactions */
/* race and place */
ods graphics on;
proc glm data=bchc_race;
   class race place;
   model mortality=race place race*place;
run;
ods graphics off;

/* race and year */
ods graphics on;
proc glm data=bchc_race;
   class race year;
   model mortality=race year race*year;
run;
ods graphics off;

/* place and year */
ods graphics on;
proc glm data=bchc_race;
   class year place;
   model mortality=year place year*place;
run;
ods graphics off;

/* Initial ANOVA model */
proc mixed data=bchc_race method=type3 cl;
    class place race year;
    model mortality = race / solution outp=pred;
    random place year race*year place*race year*place year*place*race;
run;

/* 
   We found that two-way interactions race*Year and Place*Year are not significant. Also,
   we do not have enough degrees of freedom. Therefore, the higher order interaction Place*race*Year
   should be removed.
*/

/* New ANOVA model */
proc mixed data=bchc_race method=type3 cl;
    class place race;
    model mortality = race / solution outp=pred;
    random place place*race;
run;

/* Check the normality of residuals */
proc univariate data=pred normal;
    var resid;
    histogram resid;
    probplot resid /normal(mu=est sigma=est);
run;

/* 
    Since the current residuals do not obey a normal distribution,
    we consider a transformation of the dependent variable mortality.
*/

/* Find suitable transformation function using Box-Cox */
data bchc_race_boxcox;
	set bchc_race;
	mminus2 = (-1/2) * (mortality ** -2 -1);
	mminus1= (-1) * (mortality ** -1 -1);
	mminus12= (-2) * (mortality ** -(0.5)-1);
	m0= log(mortality);
	m12= (2) * (mortality ** (1/2) -1);
	m2= (0.5) * (mortality ** (2) -1);
run;

proc univariate data=bchc_race_boxcox normal;
	histogram mminus2 /normal;
	histogram mminus1 /normal;
	histogram mminus12 /normal;
	histogram m0 /normal;
	histogram m12 /normal;
	histogram m2 /normal;
	histogram mortality /normal;
run;

/* 
    It seems log transformation is the most appropriate
*/

/* Apply Log-transformation. */
/* New dataset excluding Race_Ethnicity=All */
data bchc_race;
	set bchc_race;
	logm = log(mortality);
run;

/* New ANOVA model after log transformation */
proc mixed data=bchc_race method=type3 cl;
    class place race;
    model logm = race / solution outp=pred;
    random place place*race;
run;

/* Re-check if our residual follows a normal distribution */
proc univariate data=pred normal;
    var resid;
    histogram resid;
    probplot resid /normal(mu=est sigma=est);
run;

/* 
    The current residual meets the premise of normality
*/

/* Residual homogeneity of variance test */
proc glm data=pred;
    class race;
    model resid = race;
    means race /hovtest=Bartlett;
run;

/* 
    The null hypothesis that the variance of residual cross all groups equal is rejected.
    So we check if there is any outlier causing such a result.
*/

/* Tukey's test */
proc univariate data=bchc_race cipctldf;
    var logm;
run;

data bchc_tukey;
    set bchc_race;
    p75 = 3.45947;
	p25 = 2.67415;
	iqr = p75 - p25;
	lowert = p25 - 1.5 * iqr;
	uppert = p75 + 1.5 * iqr;
run;
data bchc_tukey;
    set bchc_tukey;
    t = (logm > uppert or logm < lowert);
run;
proc sort data=bchc_tukey;
    by descending t;
run;

proc sort data=bchc_tukey;
    by logm;
run;
proc print data=bchc_tukey;run;

/* 
    However, we could not find any outlier using Tukey's test. Now I will adjust for
    unqual variances (heteroscedasticity). I decide drop the interaction effect.
*/

/* New ANOVA model without interaction effects */
proc mixed data=bchc_race method=type3 cl;
    class race place year;
    model logm = race/ solution outp=pred;
    random place year;
run;

/* 
    It is noteworthy that the variable Year is not significant and the 
    level of between group variability is not sufficient to warrant 
    incorporating random effect in the model. So I decide to drop
    Year as a random effect.
*/

/* New New ANOVA model */
proc mixed data=bchc_race method=type3 cl;
    class race place;
    model logm = race/ solution outp=pred;
    random place;
run;

/* Re-check if our residual follows a normal distribution */
proc univariate data=pred normal;
    var resid;
    histogram resid;
    probplot resid /normal(mu=est sigma=est);
run;

/* Residual homogeneity of variance test again */
proc glm data=pred;
    class race;
    model resid = race;
    means race /hovtest=Bartlett;
run;

/* Test if the residuals and random effect coefficients are independent and identically distributed */
proc mixed data=bchc_race method=type3 cl;
    class race place;
    model logm = race/ solution outpm=rm outp=rc;
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
    by id_new place;
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

/* 
    The above results tell us that the residual satified
    the homogeneity of variance and the residuals and random
    effect coefficients are independent and identically distributed. 
    So far, all the prerequisites required by the ANOVA model
    have been satisfied. Based on the results of the existing
    ANOVA model, we obtain F-value: 144.66 and p-value < 0.0001. 
    We should reject the null hypothesis that there is no difference
    among various race groups.
    
    Our next job is to quantify the difference in expected diabetes
    disease mortality rate between the ethnicity groups studied as
    well as the uncertainty of the estimates.
*/

/* Apply Tukey's studentized range test */
proc mixed data=bchc_race method=type3;
    class race place;
    model logm = race /solution;
    random place;
    lsmeans race /diff adjust=tukey cl;
run;

/* 
    Back-transform the estimates of the difference in means in the log scale
*/

/* 
    quantify how much of the remaining variability can be explained by the
    difference between cities.
*/

/* Compute the intra-correlation coefficient (ICC) */
proc iml;
	sigG=0.06838;
	sigE=0.05679;
	icc=sigG/(sigG+sigE);
	print(ICC);
quit;

/* 
    This result tells us that there is quite moderate agreement between these cities.
*/

/* 
    How much of the variability between cities has been explained by
    a (unknown) difference in ethnicity distribution among cities.
*/

/* Check potential interactions */
/* race and place */
ods graphics on;
proc glm data=bchc_race;
   class place race;
   model logm=place race place*race;
   random place place*race;
run;
ods graphics off;












