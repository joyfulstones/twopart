libname liulei "C:\liulei\paper\paper3";
proc mixed data=liulei.all covtest;
class subject;
model logdrink=male age trt| week_offset drink_base /solution;
random intercept week_offset / type=un sub=subject;
where drink>0;
run;

* Part I model, with random intercept only;

proc nlmixed data=liulei.all;
parms alpha0=2.58 alpha1=-.215 alpha2=.0709 alpha3=-.38 alpha4=.014 alpha5=-.0799 alpha6=-.136 var1=2;
bounds var1>=0;
teta=alpha0 + a + alpha1 * male + alpha2 * age + alpha3 * drink_base + alpha4 * trt + 
	alpha5 * week_offset + alpha6 * trt* week_offset ;
expteta=exp(teta);
p=expteta / (1+expteta);
model drinkyes ~ binomial(1, p);
random a ~ normal(0, var1) subject=subject;
run;
* Part I model, with random intercept and slope;

proc nlmixed data=liulei.all qpoints=5;
parms alpha0=2.58 alpha1=-.215 alpha2=.0709 alpha3=-.38 alpha4=.014 alpha5=-.0799 alpha6=-.136 vara=6.8 varb=.5 covab=0;
bounds vara varb>=0;
teta=alpha0 + a + b* trt* week_offset + alpha1 * male + alpha2 * age + alpha3 * drink_base + alpha4 * trt + 
	alpha5 * week_offset + alpha6 * trt* week_offset ;
expteta=exp(teta);
p=expteta / (1+expteta);
model drinkyes ~ binomial(1, p);
random a b~ normal([0, 0], [vara, covab, varb]) subject=subject;
run;
* Part II model, with random intercept and slope;
proc nlmixed data=liulei.all;
parms beta0=1.5 beta1=.2 beta2=-.0 beta3=.07 beta4=-.03 beta5=-.04 beta6=-.03 vara=.065 varb=.01 covab=0 lambda=.5
		gamma0=-1 gamma1=0 gamma2=0 gamma3=0 gamma4=.2 gamma5=0 	gamma6=0;
bounds vara varb>=0;
mu=beta0 + a + b* week_offset + beta1 * male + beta2 * age + beta3 * drink_base + beta4 * trt + 
	beta5 * week_offset + beta6 * trt* week_offset ;

y=(drink ** lambda - 1) /lambda;

seta=gamma0 + gamma1 * male + gamma2 * age + gamma3 * drink_base + gamma4 * trt + 
	gamma5 * week_offset + gamma6 * trt* week_offset ;
sigma=exp(seta/2);

model y ~ normal(mu, sigma**2);
random a b~ normal([0, 0], [vara, covab, varb]) subject=subject;
where drink>0;
run;
* Box-Cox model with heteroscedasticity, random intercept and slope;
proc nlmixed data=liulei.all qpoints=5;
parms  alpha0=-.5 alpha1=-.27 alpha2=.047 alpha3=.17 alpha4=-.22 alpha5=-.1 alpha6=-.15
		beta0=2.7 beta1=0.32 beta2=-0.00049 beta3=0.005219 beta4=-0.19 beta5=-0.06 beta6=-0.03 
		gamma0=.2 gamma1=.11 gamma2=-0.02 gamma3=.04 gamma4=.09 gamma5=-.02 	gamma6=-0.01
		var1=6.7938 var2=0.51 var3=0.06293 var4=0.007 cov12=0.21 cov13=-0. cov14=-0.00037 
		cov23=-0.00311 cov24=-0.00022 cov34=0.002716 
	lambda=.37;

bounds var1 var2 var3 var4>=0;
teta=alpha0 + a + c* week_offset + alpha1 * male + alpha2 * age + alpha3 * drink_base + alpha4 * trt + 
	alpha5 * week_offset + alpha6 * trt* week_offset ;
expteta=exp(teta);
p=expteta / (1+expteta);

mu=beta0 + b + d* week_offset + beta1 * male + beta2 * age + beta3 * drink_base + beta4 * trt + 
	beta5 * week_offset + beta6 * trt* week_offset ;

if drink=0 then loglik=log(1-p);
if drink>0 then do;

	seta=gamma0 + gamma1 * male + gamma2 * age + gamma3 * drink_base + gamma4 * trt + 
	gamma5 * week_offset + gamma6 * trt* week_offset ;
	sigma=exp(seta/2);

	y=(drink ** lambda - 1) /lambda;

	loglik=log(p) - .5*((y-mu)/sigma)**2-.5*log(2*3.1415926*sigma**2) + (lambda -1) * log(drink);
end;
model drink ~ general(loglik);
random a b c d~ normal([0, 0, 0, 0], [var1, cov12, var2, cov13, cov23, var3, cov14, cov24, cov34, var4]) subject=subject;
estimate 'rho12' cov12/sqrt(var1*var2);
estimate 'rho13' cov13/sqrt(var1*var3);
estimate 'rho14' cov14/sqrt(var1*var4);
estimate 'rho23' cov23/sqrt(var2*var3);
estimate 'rho24' cov24/sqrt(var2*var4);
estimate 'rho34' cov34/sqrt(var3*var4);
ods output ParameterEstimates=est1_all FitStatistics=fit1_all Additionalestimates=est2_all; 
run;

* 2 part random effects model, with random intercept;
proc nlmixed data=liulei.all qpoints=5;
parms alpha0=1.5 alpha1=-.27 alpha2=.03 alpha3=.05 alpha4=-.54 alpha5=-.14 alpha6=-.08 
beta0=1.32 beta1=.1 beta2=-.0 beta3=.09 beta4=-.03 beta5=-.02 beta6=-.02 
sigma_e=.35 var1=6.8 var2=.08 cov12=0;
bounds var1 var2>=0;
teta=alpha0 + a + alpha1 * male + alpha2 * age + alpha3 * drink_base + alpha4 * trt + 
	alpha5 * week_offset + alpha6 * trt* week_offset ;
expteta=exp(teta);
p=expteta / (1+expteta);
mu=beta0 + b + beta1 * male + beta2 * age + beta3 * drink_base + beta4 * trt + 
	beta5 * week_offset + beta6 * trt* week_offset ;
if drink=0 then loglik=log(1-p);
if drink>0 then loglik=log(p)-.5*((log(drink)-mu)/sigma_e)**2-.5*log(2*3.14159*sigma_e**2);
model drink ~ general(loglik);
random a b ~ normal([0, 0], [var1, cov12, var2]) subject=subject;
ods output ParameterEstimates=est1_all FitStatistics=fit1_all; 
estimate 'sigma_a' sqrt(var1);
estimate 'sigma_b' sqrt(var2);
estimate 'rho' cov12/sqrt(var1*var2);
run;
* 2 part random effects model, with random intercept and slope, qpoints=5;
proc nlmixed data=liulei.all qpoints=5;
parms alpha0=.6 alpha1=-.18 alpha2=.048 alpha3=.1 alpha4=-.28 alpha5=-.14 alpha6=-.155 
beta0=1.33 beta1=.156 beta2=-.0 beta3=.084 beta4=-.069 beta5=-.02 beta6=-.19 
sigma_e=.59 var1=6.7 var2=.19 var3=.13 var4=.002 cov12=-.33 cov13=.013 cov14=0 cov23=0 cov24=0 cov34=0 ;
bounds var1 var2 var3 var4>=0;
teta=alpha0 + a + b* week_offset + alpha1 * male + alpha2 * age + alpha3 * drink_base + alpha4 * trt + 
	alpha5 * week_offset + alpha6 * trt* week_offset ;
expteta=exp(teta);
p=expteta / (1+expteta);
mu=beta0 + c + d* week_offset + beta1 * male + beta2 * age + beta3 * drink_base + beta4 * trt + 
	beta5 * week_offset + beta6 * trt* week_offset ;
if drink=0 then loglik=log(1-p);
if drink>0 then loglik=log(p)-.5*((log(drink)-mu)/sigma_e)**2-.5*log(2*3.14159*sigma_e**2);
model drink ~ general(loglik);
random a b c d~ normal([0, 0, 0, 0], [var1, cov12, var2, cov13, cov23, var3, cov14, cov24, cov34, var4]) subject=subject;
ods output ParameterEstimates=est1_all FitStatistics=fit1_all; 
estimate 'sigma_a' sqrt(var1);
estimate 'sigma_b' sqrt(var2);
estimate 'rho' cov12/sqrt(var1*var2);
run;
* 2 part random effects model, with random intercept and slope, qpoints=10;
proc nlmixed data=liulei.all qpoints=10;
parms alpha0=.6 alpha1=-.18 alpha2=.043 alpha3=.12 alpha4=-.28 alpha5=-.11 alpha6=-.18 
beta0=1.33 beta1=.12 beta2=-.0 beta3=.092 beta4=-.058 beta5=-.026 beta6=-.028 
sigma_e=.57 var1=6.7 var2=.067 var3=.13 var4=.002 cov12=-.11 cov13=-.033 cov14=0 cov23=0 cov24=0.0057 cov34=-.005 ;
bounds var1 var2 var3 var4>=0;
teta=alpha0 + a + b* week_offset + alpha1 * male + alpha2 * age + alpha3 * drink_base + alpha4 * trt + 
	alpha5 * week_offset + alpha6 * trt* week_offset ;
expteta=exp(teta);
p=expteta / (1+expteta);
mu=beta0 + c + d* week_offset + beta1 * male + beta2 * age + beta3 * drink_base + beta4 * trt + 
	beta5 * week_offset + beta6 * trt* week_offset ;
if drink=0 then loglik=log(1-p);
if drink>0 then loglik=log(p)-.5*((log(drink)-mu)/sigma_e)**2-.5*log(2*3.14159*sigma_e**2);
model drink ~ general(loglik);
random a b c d~ normal([0, 0, 0, 0], [var1, cov12, var2, cov13, cov23, var3, cov14, cov24, cov34, var4]) subject=subject;
ods output ParameterEstimates=est1_all FitStatistics=fit1_all; 
estimate 'sigma_a' sqrt(var1);
estimate 'sigma_b' sqrt(var2);
estimate 'rho' cov12/sqrt(var1*var2);
run;

