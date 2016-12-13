libname liulei "C:\liulei\paper\paper3";
proc mixed data=liulei.all covtest;
class subject;
model logdrink=male age trt| week_offset drink_base /solution;
random intercept week_offset / type=un sub=subject;
where drink>0;
run;

* Part I model, with random intercept only;

proc nlmixed data=liulei.all qpoints=5;
parms alpha0=2.58 alpha1=-.215 alpha2=.0709 alpha3=-.38 alpha4=.014 alpha5=-.0799 alpha6=-.136 var1=2;
bounds var1>=0;
teta=alpha0 + a + alpha1 * male + alpha2 * age + alpha3 * drink_base + alpha4 * trt + 
	alpha5 * week_offset + alpha6 * trt* week_offset ;
expteta=exp(teta);
p=expteta / (1+expteta);
model drinkyes ~ binomial(1, p);
random a ~ normal(0, var1) subject=subject;
predict a out=a1;
run;
* Part I model, with random intercept and slope;

proc nlmixed data=liulei.all qpoints=5;
parms alpha0=2.58 alpha1=-.215 alpha2=.0709 alpha3=-.38 alpha4=.014 alpha5=-.0799 alpha6=-.136 vara=6.8 varb=.5 covab=0;
bounds vara varb>=0;
teta=alpha0 + a + b * week_offset + alpha1 * male + alpha2 * age + alpha3 * drink_base + alpha4 * trt + 
	alpha5 * week_offset + alpha6 * trt* week_offset ;
expteta=exp(teta);
p=expteta / (1+expteta);
model drinkyes ~ binomial(1, p);
random a b~ normal([0, 0], [vara, covab, varb]) subject=subject;
predict a out=a;
predict b out=b;
ods output FitStatistics=fit;
ods output ParameterEstimates=par;
run;
* Part I model, with random intercept and slope;
proc nlmixed data=liulei.all qpoints=5;
parms alpha0=-0.5202 alpha1=-0.2693 alpha2=0.04828 alpha3=0.1640 alpha4=-0.2167 
alpha5=-0.1009 alpha6=-0.1744 vara=5.7012 varb=0.06633 covab=-0.08345;
bounds vara varb>=0;
teta=alpha0 + a + b * week_offset + alpha1 * male + alpha2 * age + alpha3 * drink_base + alpha4 * trt + 
	alpha5 * week_offset + alpha6 * trt* week_offset ;
expteta=exp(teta);
p=expteta / (1+expteta);
model drinkyes ~ binomial(1, p);
random a b~ normal([0, 0], [vara, covab, varb]) subject=subject;
predict a out=a;
predict b out=b;
ods output FitStatistics=fit;
ods output ParameterEstimates=par;
run;
* Part II model, with random intercept and slope;
proc nlmixed data=liulei.all;
parms beta0=1.32 beta1=.1 beta2=-.0 beta3=.09 beta4=-.03 beta5=-.02 beta6=-.02 vara=.13 varb=.002 covab=0;
bounds vara varb>=0;
mu=beta0 + a + b* week_offset + beta1 * male + beta2 * age + beta3 * drink_base + beta4 * trt + 
	beta5 * week_offset + beta6 * trt* week_offset ;
model logdrink ~ normal(mu, sigma_e**2);
random a b~ normal([0, 0], [vara, covab, varb]) subject=subject;
where drink>0;
run;
* generalized Gamma model with heteroscedasticity, random intercept only;
proc nlmixed data=liulei.all qpoints=5;
parms  alpha0=1.5 alpha1=-.27 alpha2=.03 alpha3=.05 alpha4=-.54 alpha5=-.14 alpha6=-.08
		beta0=1.32 beta1=.1 beta2=-.0 beta3=.09 beta4=-.03 beta5=-.02 beta6=-.02 
		gamma0=.5 gamma1=0 gamma2=0 gamma3=0 gamma4=0 gamma5=0 	gamma6=0
		var1=6.8 var2=.08 cov12=.024 k=.629;
bounds var1 var2 >=0;
teta=alpha0 + a + alpha1 * male + alpha2 * age + alpha3 * drink_base + alpha4 * trt + 
	alpha5 * week_offset + alpha6 * trt* week_offset ;
expteta=exp(teta);
p=expteta / (1+expteta);

mu=beta0 + b + beta1 * male + beta2 * age + beta3 * drink_base + beta4 * trt + 
	beta5 * week_offset + beta6 * trt* week_offset ;

eta=abs(k) ** (-2);
if drink=0 then loglik=log(1-p);
if drink>0 then do;

	seta=gamma0 + gamma1 * male + gamma2 * age + gamma3 * drink_base + gamma4 * trt + 
	gamma5 * week_offset + gamma6 * trt* week_offset ;
	sigma=exp(seta/2);
	value1=eta *log (eta) - log(drink) - log(sigma) -.5 * log(eta) - lgamma(eta);

	u=sign(k)*(log(drink)-mu)/sigma; 
	loglik=log(p) + value1 + u *sqrt(eta) - eta * exp(abs(k)* u);
end;
model drink ~ general(loglik);
random a b ~ normal([0, 0], [var1, cov12, var2]) subject=subject;
estimate 'gamma' sigma-k;
estimate 'Weibull' k-1;
estimate 'lognormal' k;
estimate 'Correlation' cov12/sqrt(var1*var2);
estimate 'shape' k ** (-2);
predict a out=a;
predict b out=b;
ods output ParameterEstimates=est1_all FitStatistics=fit1_all Additionalestimates=est2_all; 
run;
* generalized Gamma model with heteroscedasticity, random intercept and slope;
proc nlmixed data=liulei.all qpoints=5;
parms  alpha0=-.5 alpha1=-.27 alpha2=.047 alpha3=.17 alpha4=-.22 alpha5=-.1 alpha6=-.15
		beta0=1.8 beta1=.21 beta2=-.0 beta3=.066 beta4=.054 beta5=-.03 beta6=-.037 
		gamma0=-1.12 gamma1=-.13 gamma2=-.03 gamma3=-.03 gamma4=0.21 gamma5=-.01 gamma6=.03
		var1=5.7 var2=.10 var3=.064 var4=.002 cov12=-.055 cov13=-.095 cov14=0 cov23=0 cov24=0 cov34=0.006 
k=.71;

bounds var1 var2 var3 var4>=0;
teta=alpha0 + a + c* week_offset + alpha1 * male + alpha2 * age + alpha3 * drink_base + alpha4 * trt + 
	alpha5 * week_offset + alpha6 * trt* week_offset ;
expteta=exp(teta);
p=expteta / (1+expteta);

mu=beta0 + b + d* week_offset + beta1 * male + beta2 * age + beta3 * drink_base + beta4 * trt + 
	beta5 * week_offset + beta6 * trt* week_offset ;

eta=abs(k) ** (-2);
if drink=0 then loglik=log(1-p);
if drink>0 then do;

	seta=gamma0 + gamma1 * male + gamma2 * age + gamma3 * drink_base + gamma4 * trt + 
	gamma5 * week_offset + gamma6 * trt* week_offset ;
	sigma=exp(seta/2);
	value1=eta *log (eta) - log(drink) -log(sigma) -.5 * log(eta) - lgamma(eta);

	u=sign(k)*(log(drink)-mu)/sigma; 
	loglik=log(p) + value1 + u *sqrt(eta) - eta * exp(abs(k)* u);
end;
model drink ~ general(loglik);
random a b c d~ normal([0, 0, 0, 0], [var1, cov12, var2, cov13, cov23, var3, cov14, cov24, cov34, var4]) subject=subject;
estimate 'gamma' sigma-k;
estimate 'Weibull' k-1;
estimate 'lognormal' k;
estimate 'rho12' cov12/sqrt(var1*var2);
estimate 'rho13' cov13/sqrt(var1*var3);
estimate 'rho14' cov14/sqrt(var1*var4);
estimate 'rho23' cov23/sqrt(var2*var3);
estimate 'rho24' cov24/sqrt(var2*var4);
estimate 'rho34' cov34/sqrt(var3*var4);
estimate 'shape' k ** (-2);
predict a out=a;
predict b out=b;
ods output ParameterEstimates=est1_all FitStatistics=fit1_all Additionalestimates=est2_all; 
run;
