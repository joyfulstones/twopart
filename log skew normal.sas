libname liulei "C:\liulei\paper\paper3";
* log skew normal model with heteroscedasticity, random intercept only;
proc nlmixed data=liulei.all qpoints=5;
parms  alpha0=1.5 alpha1=-.27 alpha2=.03 alpha3=.05 alpha4=-.54 alpha5=-.14 alpha6=-.08
		beta0=1.32 beta1=.1 beta2=-.0 beta3=.09 beta4=-.03 beta5=-.02 beta6=-.02 
		gamma0=.5 gamma1=0 gamma2=0 gamma3=0 gamma4=0 gamma5=0 	gamma6=0
		var1=6.8 var2=.08 cov12=.024 lambda=0;
bounds var1 var2 >=0;
teta=alpha0 + a + alpha1 * male + alpha2 * age + alpha3 * drink_base + alpha4 * trt + 
	alpha5 * week_offset + alpha6 * trt* week_offset ;
expteta=exp(teta);
p=expteta / (1+expteta);

mu=beta0 + b + beta1 * male + beta2 * age + beta3 * drink_base + beta4 * trt + 
	beta5 * week_offset + beta6 * trt* week_offset ;

if drink=0 then loglik=log(1-p);
if drink>0 then do;

	seta=gamma0 + gamma1 * male + gamma2 * age + gamma3 * drink_base + gamma4 * trt + 
	gamma5 * week_offset + gamma6 * trt* week_offset ;
	sigma=exp(seta/2);

	value1=sigma**2 + lambda **2;
	value2=(log(drink)-mu) / sqrt(value1);
	value3=value2 * lambda/sigma;

	logden=log(2/drink) -.5 * log(value1) - .5 *log(2*3.1415926) -.5 *value2**2 + log(CDF('NORMAL',value3)) ;

	loglik=log(p) + logden;
end;
model drink ~ general(loglik);
random a b ~ normal([0, 0], [var1, cov12, var2]) subject=subject;
*estimate 'gamma' sigma-k;
estimate 'lognormal' lambda;
estimate 'Correlation' cov12/sqrt(var1*var2);
predict a out=a;
predict b out=b;
ods output ParameterEstimates=est1_all FitStatistics=fit1_all Additionalestimates=est2_all; 
run;

* log skew normal model with heteroscedasticity, random intercept and slope;
proc nlmixed data=liulei.all qpoints=5;
parms  alpha0=-.5 alpha1=-.27 alpha2=.047 alpha3=.17 alpha4=-.22 alpha5=-.1 alpha6=-.15
		beta0=1.7 beta1=.17 beta2=-.0 beta3=.066 beta4=-.044 beta5=-.035 beta6=-.03 
		gamma0=-1.06 gamma1=-.12 gamma2=-.03 gamma3=-.03 gamma4=0.2 gamma5=-.01 	gamma6=0.03
		var1=5.7 var2=.12 var3=.065 var4=.002 cov12=-.055 cov13=-.095 cov14=-.007 cov23=-.002 cov24=-.002 cov34=.006 
lambda=0;

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

	value1=sigma**2 + lambda **2;
	value2=(log(drink)-mu) / sqrt(value1);
	value3=value2 * lambda/sigma;

	logden=log(2/drink) -.5 * log(value1) - .5 *log(2*3.1415926) -.5 *value2**2 + log(CDF('NORMAL',value3)) ;

	loglik=log(p) + logden;
end;
model drink ~ general(loglik);
random a b c d~ normal([0, 0, 0, 0], [var1, cov12, var2, cov13, cov23, var3, cov14, cov24, cov34, var4]) subject=subject;
estimate 'lognormal' lambda;
estimate 'rho12' cov12/sqrt(var1*var2);
estimate 'rho13' cov13/sqrt(var1*var3);
estimate 'rho14' cov14/sqrt(var1*var4);
estimate 'rho23' cov23/sqrt(var2*var3);
estimate 'rho24' cov24/sqrt(var2*var4);
estimate 'rho34' cov34/sqrt(var3*var4);
predict a out=a;
predict b out=b;
ods output ParameterEstimates=est1_all FitStatistics=fit1_all Additionalestimates=est2_all; 
run;

