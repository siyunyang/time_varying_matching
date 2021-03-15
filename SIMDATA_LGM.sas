
*Simulation resembling Taylor et al. 2013; 
*Set parameters;
*To reduce noise and get results similar to Taylor et al. use N=10000;
*To run quickly use N<1000;
%let N = 1000; 
%let alpha0 = 1.635;  
%let alpha11 = 2.443;  
%let alpha12 = 0.217;
%let alpha13 = 0.249;
%let alpha21 = 0.242;
%let alpha22 = 0.224;
%let alpha23 = 0.547;
%let sigma2 = 0.061; 
%let beta0 = -7.258;
%let beta11 = -0.036;
%let beta12 = -0.022; 
%let beta3 = 0.740;
%let lambda0 = 0.007503;
%let theta01 = 0.812;
%let theta02 = 0.918;
%let theta1 = 0.050;
%let theta2 = 2.018;
%let gamma = -1.5;
%let endtime = 120;

*generate random effects for subjects;
data SIGMA (type=cov);
   input _TYPE_ $ 1-4 _NAME_ $ 6-8 a0 a1 a2;
   datalines;
cov  a0 1.084 1.065 0.148
cov  a1 1.065 2.658 0.456
cov  a2 0.148 0.456 0.322
mean     0     0     0
run;
*Perfect number is 500 and seed=654;
proc simnormal data=SIGMA out=temp numreal=&N. seed=654;
	var a0 a1 a2;
run;

*create a subject identifier;
data RE; set temp;
retain subject;
if _n_ = 1 then subject=0;
subject = subject+1;
run;

data simdata; set RE;
	call streaminit(654); 
	age_0 = 70 + 6 * rand("normal");   
	if rand("uniform") < 0.309 then X = 1;
		else if rand("uniform") < 0.617 / 0.691 then X = 2;
		else X = 3;
	S_t=0;
	S_start=.;
	do time = 1 to &endtime. by 1;
		T = time/10;

		*Model for PSA values over time;
		log_P_t = (&alpha0. + a0) + 
				  (&alpha11. + &alpha12.*(X=2) + &alpha13.*(X=3) + a1) * ((1+T)**(-1.5)-1) +
				  (&alpha21. + &alpha22.*(X=2) + &alpha23.*(X=3) + a2) * T +
				   sqrt(&sigma2.)*rand("normal"); 

		*Model for treatment with SADT, once treated always treated;
		if S_t =0 then do;
	 	  	lpS =  &beta0. + &beta11.*(X=2) + &beta12.*(X=3) + &Beta3.*log_P_t;
	      	pS = 1/(1 + exp(-lpS));
	      	S_t = rand('bern',pS);
			if S_t = 1 then S_start = time;
		end;

		*Model for recurrence;
		slope_P_t = -1.5*(1+T)**(-2.5)*(&alpha11. + &alpha12.*(X=2) + &alpha13.*(X=3) + a1) + 
					(&alpha21. + &alpha22.*(X=2) + &alpha23.*(X=3) + a2);
		*Discrete time adaptation of the survival model in paper;
	 	lpR =  -7.19 + &theta01.*(X=2) + &theta02.*(X=3) + &theta1.*log_P_t + 
					&theta2.*slope_P_t + &gamma.*S_t;
	  	pR = 1/(1 + exp(-lpR));
		R_t = rand('bern',pR);

		if time=1 then do; log_P_0 = log_P_t; slope_P_0 = slope_P_t; end;
        time0=time-1;
		output;
		if R_t = 1 then leave; 
	end; 
run;
proc print data = simdata (obs=500);
run;
 
*Output status at last time point;
data summary;
	set simdata;
	by subject T;
	if last.subject;
run;
*Compare to Taylor et al. Section 4.5; 
proc means data=summary mean sum maxdec=3;
	var S_t R_t;
run; 

***********************************************************;
*Note: The simulation of Taylor et al. augmented matching methods
by adding a random effects model for longitudinal PSA values, 
slope and patient-specific intercept. Estimates from these models
replaced the longitudinal data points. In Section 3.3 they discuss the fact
that this is an addition, and not a necessary feature. Instead, it would
also be possible to use the raw values of PSA. Here we use 
the known value slope_P_t for simplicity though in practice it would need
to be estimated. Results are similar to estimated version in Taylor et al.;
***********************************************************;
 
********************************************************************************************;
* Longitudinal matching as in Section 5.2.1 of Thomas et al. 2020
* In application to simulated data like Taylor et al. 2013
*********************************************************************************************;

*repeat data over multiple experiments j;
	*SADT holds the treatment status as baseline for experiment j;
	*All subsequent records are output for all patients in experiment j;
	*log_P_j and slope_P_j hold the baseline covariate information at experiment j;
data long; set simdata;
do j = 1 to &endtime. by 1;
	if (S_start >= j or S_start = .) then do;
		if time = j then do;
			SADT=S_t;
			X_j = X;
			log_P_j = log_P_t;  
			slope_P_j = slope_P_t;
			output;
		end;
		if time > j then do;  
			output;
		end;
	end;
end;
keep j subject time time0 S_t S_start SADT log_P_j X_j R_t log_P_t X slope_P_j slope_P_t;
run;
proc sort data = long; by j subject time; run;
proc print data = long (obs=500); where time>60 and time<70 or time=1; run;


*Retain baseline values at time j;
data long2; set long;
by j subject time;
retain SADT0 log_P_j0 slope_P_j0 X_j0;
if not missing(SADT) then SADT0=SADT;
else SADT=SADT0;
if not missing(log_P_j) then log_P_j0=log_P_j;
else log_P_j=log_P_j0;
if not missing(slope_P_j) then slope_P_j0=slope_P_j;
else slope_P_j=slope_P_j0;
if not missing(X_j) then X_j0 = X_j;
else X_j = X_j0;
drop SADT0 log_P_j0 slope_P_j0 X_j0;
run;
proc print data = long2 (obs=500); where time>60 and time<70 or time=1; run;


*******************************************************;
*Intent to treat analysis similar to Section 5.2.1   ;
proc phreg data = long2 covs(aggregate);  
strata j;
effect Pspl=spline(log_P_j / naturalcubic details knotmethod=percentiles(3));   
effect Sspl=spline(slope_P_j / naturalcubic details knotmethod=percentiles(3));  
class X_j;
model (time0,time)*R_t(0) =  SADT Pspl Sspl X_j ;
id subject;
run; 
 
*******************************************************;
*Per-protocol analysis similar to Section 5.2.1 (but without adding weights);
*censoring follow-up if a control patient starts treatment (note that treated can not crossover);
data pp_data; set long2;
if SADT=0 and S_t = 1 then delete;
run; 

proc phreg data = pp_data covs(aggregate);  
strata j;
effect Pspl=spline(log_P_j / naturalcubic details knotmethod=percentiles(3));   
effect Sspl=spline(slope_P_j / naturalcubic details knotmethod=percentiles(3));  
class j X_j;
model (time0,time)*R_t(0) =  SADT Pspl Sspl X_j ;
id subject;
run; 

*******************************************************;
*Per-protocol analysis with IPCW weights similar to Section 5.2.1 ;
* START EDITING HERE;
*Prep data for modeling censoring (due to switching) with time-varying covariates;
*i.e. the probability of remaining untreated given that the patient was untreated at the 
start of experiment;
data censor; set long2;
where SADT=0; *only this group has non-administrative censoring;
C=0;
if S_start =time then C=1;
if S_start ne . and S_start<time then delete;
run;
proc print data = censor (obs=500); where time>60 and time<70 or time=1; run;

*Baseline model for stabilization;
proc logistic data=censor ; 
   effect timespl=spline(time / naturalcubic details knotmethod=percentiles(4));   
   effect Pspl=spline(log_P_j / naturalcubic details knotmethod=percentiles(3));   
   effect Sspl=spline(slope_P_j / naturalcubic details knotmethod=percentiles(3));  
   class j X_j / param=ref ref=first;
   model C = j timespl Pspl Sspl X_j ;
   output out=model0 (keep= j subject time pC_0) p=pC_0; *probability of remaining uncensored;
run;
*Model for time-dependent censoring;
proc logistic data=censor ; 
   effect timespl=spline(time / naturalcubic details knotmethod=percentiles(4));   
   effect Pspl=spline(log_P_t / naturalcubic details knotmethod=percentiles(3));   
   effect Sspl=spline(slope_P_t / naturalcubic details knotmethod=percentiles(3));  
   class j X / param=ref ref=first;
   model C = j timespl Pspl Sspl X ;
   output out=model1 (keep= j subject time pC_w) p=pC_w; *probability of remaining uncensored;
run;
*Creation of time-varying weights;
proc sort data = model0; by j subject time; run;
proc sort data = model1; by j subject time; run;
proc sort data = pp_data; by j subject time; run;
data main_w;
merge model0 model1 pp_data (in=a);  
by j subject time;
if a;
if pC_0=. then pC_0=1; *for the treated group that had no censoring;
if pC_w=. then pC_w=1; *for the treated group that had no censoring;
/* variables ending with _0 refer to the numerator of the weights
variables ending with _w refer to the denominator of the weights */
if first.subject then do;
	stabw=1;
end;
retain stabw;
/* Inverse probability of censoring weights */
ratio= pC_0/pC_w;
stabw = stabw*ratio;
run;
proc print data = main_w (obs=500); where time>60 and time<70 or time=1; run;


*Review the weights for extreme values;
proc means data = main_w min q1 median q3 max maxdec=1;
class time;
var stabw;
run; 

*Per-protocol analysis with IPCW weights;
proc phreg data = main_w covs(aggregate);  
strata j;
effect Pspl=spline(log_P_j / naturalcubic details knotmethod=percentiles(3));   
effect Sspl=spline(slope_P_j / naturalcubic details knotmethod=percentiles(3));  
class j X_j;
model (time0,time)*R_t(0) =  SADT Pspl Sspl X_j ;
id subject;
weight stabw;
run; 


********************************************************************************************;
* Longitudinal matching as in Section 5.2.2 of Thomas et al. 2020
* In application to simulated data like Taylor et.al.
*********************************************************************************************;
*NOTE: We create strata based on eligibility (being at risk) and exact agreement on
time and tumor type X;
*NOTE: PSA values are accounted for in the model; 

*Creating strata over time and our strata variable of interest X;
data long2; set long2;
length strata $32.;
PSAcat = round(log_P_j/0.5,1);
strata = j || "_" || PSAcat;
run;
 
*Exclude strata without any comparators;
	*Treated without any controls in the same strata;
	*Untreated without any treated in the same strata;
*Only strata with both treated and untreated will remain;
proc sort data=long2 out=temp ; 
where j=time;
by strata SADT; 
run; 
data treated; set temp;  
by strata SADT;
if last.strata and SADT =1 then output;  
rename SADT=SADT1;
rename subject=subject1;
keep subject strata SADT j;
run;
data control; set temp;  
by strata SADT;
if first.strata and SADT=0 then output;  
rename SADT=SADT0;
rename subject=subject0;
keep subject strata SADT j;
run;
proc sort data = control; by strata; run;
proc sort data = treated; by strata; run;
data combine; merge control (in=a) treated (in=b); 
by strata ;
if a and b ;
run;  
proc sort data = long2; by strata; run;
proc sort data = combine; by strata; run;
*(1) including all possible eligible controls;
data ss_data; merge long2 combine (in=a); 
by strata;
if a;
run; 

*Different options for how many controls to allow in a strata;
*This was not proposed in the original method but is considered here as an alternative;

* (2) Further limiting to only 1 control and 1 treated;
data ss_data2; set ss_data;
where subject = subject0 or subject=subject1;
run;
*NOTE: Inclusion of many controls makes the result much worse 
*this is not a general feature of the methods but can be explained in this specific context;
*Essentially more controls means more model-based adjustment and dependence on model-specification;
 
* (3) Including up to 2 controls per case;
proc sort data = combine; by strata ; run;
proc sort data = temp; by strata; run;
data count; merge combine (in=a) temp;
by strata;
if a;
run;
proc sort data= count; by strata descending SADT subject ;
run;
data temp; set count;
by strata descending SADT subject;
if first.strata then count=0;
count+1;
run; 
data limit; set temp;
where count<4;
keep strata subject;
run;
proc sort data = long2; by strata subject; run;
proc sort data = limit; by strata subject; run;
data ss_data3; merge long2 limit (in=a); 
by strata subject;
if a;
run;  


*******************************************************;
*Intent to treat analysis similar to Section 5.2.2  ;
*Using all available controls;
proc phreg data = ss_data covs(aggregate);  
strata strata;  
class strata X_j;
model (time0,time)*R_t(0) =  SADT log_P_j slope_P_j  X_j;
id subject;
run; 
  
*This version is more similar to Taylor et al. 2013 because that paper adapted
*sequential stratification in a number of unique ways so that it became more like 
*propensity score matching and each case was matched with only a few controls;
*These adaptations were unique to Taylor et al. 2013 and not a fundamental part of the proposed method;
proc phreg data = ss_data3 covs(aggregate);  
strata strata;  
effect tspl=spline(time0 / naturalcubic details knotmethod=percentiles(3)); 
class strata X_j;
model (time0,time)*R_t(0) =  SADT log_P_j log_P_j*tspl slope_P_j slope_P_j*tspl X_j;
id subject;
run; 

*******************************************************;
*Per-protocol analysis similar to Section 5.2.2 (but without adding weights);
*censoring follow-up if a control patient starts treatment (note that treated can not crossover);
data pp_data; set ss_data3;
if SADT=0 and S_t=1 then delete;
run;
proc phreg data = pp_data covs(aggregate);  
strata strata;  
effect tspl=spline(time0 / naturalcubic details knotmethod=percentiles(3)); 
class strata X_j;
model (time0,time)*R_t(0) =  SADT log_P_j log_P_j*tspl slope_P_j slope_P_j*tspl X_j;
id subject;
run; 

*******************************************************;
*Per-protocol analysis with IPCW weights similar to Section 5.2.2 ;

*Prep data for modeling censoring with time-varying covariates;
*This only occurs on the untreated side, otherwise censoring is only at 12 years; 
data censor; set ss_data;
where SADT=0;
C=0;
if S_start =time then C=1;
if S_start ne . and S_start <time then delete;
run; 
*Baseline model for stabilization;
proc logistic data=censor ; 
   effect timespl=spline(time / naturalcubic details knotmethod=percentiles(3));   
   effect Pspl=spline(log_P_j / naturalcubic details knotmethod=percentiles(3));   
   effect Sspl=spline(slope_P_j / naturalcubic details knotmethod=percentiles(3));  
   class j X_j / param=ref ref=first;
   model C = j timespl SADT Pspl Sspl X_j ;
   output out=model0 (keep= j subject time pC_0) p=pC_0; *pC_0 is the probability of being uncensored;
run;
*Model for time-dependent censoring;
proc logistic data=censor ; 
   effect timespl=spline(time / naturalcubic details knotmethod=percentiles(3));   
   effect Pspl=spline(log_P_t / naturalcubic details knotmethod=percentiles(3));   
   effect Sspl=spline(slope_P_t / naturalcubic details knotmethod=percentiles(3));  
   class j X / param=ref ref=first;
   model C = j timespl SADT Pspl Sspl X ;
   output out=model1 (keep= j subject time pC_w) p=pC_w; *pC_w is the probability of being uncensored;
run;
*Creation of time-varying weights;
proc sort data = model0; by j subject time; run;
proc sort data = model1; by j subject time; run;
proc sort data = pp_data; by j subject time; run;
data main_w;
merge model0 model1 pp_data (in=a);  
by j subject time;
if a;
if pC_w=. then pC_w = 1;
if pC_0=. then pC_0 = 1;
/* variables ending with _0 refer to the numerator of the weights
variables ending with _w refer to the denominator of the weights */
if first.subject then do;
	stabw=1;
end;
retain stabw;
/* Inverse probability of censoring weights */
ratio= pC_0/pC_w;
stabw = stabw*ratio;
run;
*Review the weights for extreme values;
proc means data = main_w min q1 median q3 max maxdec=1;
class SADT time;
var stabw;
run; 

*Per-protocol analysis with IPCW weights;
proc phreg data = main_w covs(aggregate);  
strata strata;  
effect tspl=spline(time0 / naturalcubic details knotmethod=percentiles(3)); 
class strata X_j;
model (time0,time)*R_t(0) =  SADT log_P_j log_P_j*tspl slope_P_j slope_P_j*tspl X_j;
weight stabw;
id subject;
run; 
 

********************************************************************************************;
* Longitudinal matching as in Section 5.2.3 of Thomas et al. 2020
* In application to simulated data like Taylor et.al.
*********************************************************************************************;

*time-dependent propensity model to start treatment;
proc phreg data = simdata; 
where S_start >= time or S_start=.; *excluding times post-treatment initiation;
effect Pspl=spline(log_P_t / naturalcubic details knotmethod=percentiles(3));   
effect Sspl=spline(slope_P_t / naturalcubic details knotmethod=percentiles(3));  
class X;
model (time0,time)*S_t(0) =  Pspl Sspl X;
output out=propensity xbeta=xbeta;
run;
*potential matches over j (output at time j only);
data long3; set propensity;
do j = 1 to &endtime. by 1;
	if (S_start >= j or S_start = .) then do;
		if time = j then do;
			SADT=S_t;
			X2_j = (X=2);
			X3_j = (X=3);
			X_j = X;
			log_P_j = log_P_t;  
			slope_P_j = slope_P_t;
			output;
		end; 
	end;
end;
keep j subject time time0 S_start SADT log_P_j X2_j X3_j X_j R_t log_P_t X slope_P_j slope_P_t prior xbeta;
run;

*SAS macro for matching requires a ``propensity score'' between 0 and 1;
proc sql noprint;
       select max(xbeta) into: max_xbeta from long3; 
quit;
data long3; set long3; 
xbeta2 = (&max_xbeta.+1-xbeta)/100;
run; 
*Matching on Xbeta from the hazard model;
	*Match exactly on j;
	*Match within a caliper on the linear predictor (propensity);
proc psmatch data=long3 ; 
   class j SADT X2_j X3_j;
   psdata treatvar=SADT(Treated='1') ps=xbeta2;
   match exact=j method=GREEDY(K=1 ORDER=descending ) stat=lps caliper=.25;
   assess ps var=(log_P_j slope_P_j X2_j X3_j) / weight=none;
   output out(obs=match)=matched_sample  ps=ps matchid=_MatchID;
run;
*Merging the outcome data back onto the matched sample;
proc sort data = matched_sample; by j subject time; run;
proc sort data = long2; by j subject time; run;
data match_data; merge matched_sample (in=a) long2;
by j subject;
if a;
run;

*******************************************************;
*Intent to treat analysis (doubly robust version) similar to Section 5.2.3 ;
proc phreg data = match_data covs(aggregate);  
strata _MatchID;   
effect tspl=spline(time0 / naturalcubic details knotmethod=percentiles(3)); 
class _MatchID X_j;
model (time0,time)*R_t(0) =  SADT log_P_j log_P_j*tspl slope_P_j slope_P_j*tspl X_j;
id subject;
run; 
 
 
*******************************************************;
*Per-protocol analysis similar to Section 5.2.3 (but without adding weights);
*censoring follow-up if a control patient starts treatment (note that treated can not crossover);
data pp_data; set match_data;
if SADT=0 and S_t=1 then delete;
run; 
proc phreg data = pp_data covs(aggregate);  
strata _MatchID;
effect tspl=spline(time0 / naturalcubic details knotmethod=percentiles(3)); 
class _MatchID X_j;
model (time0,time)*R_t(0) =  SADT log_P_j log_P_j*tspl slope_P_j slope_P_j*tspl X_j;
id subject;
run; 
