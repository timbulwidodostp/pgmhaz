version 11.2
mata:

void pgmhaz_gf2mataviewo(transmorphic M,
					real scalar todo, 
					real rowvector b, 
					real colvector fv, 
					real matrix S,
					real matrix H)
{

// y1=1 if an individual reports failure, 0 if right censored

// PARAMETERS TO BE ESTIMATED
// plva: log variance scalar
// ptco: time constance covariates column vector with
//       each row equal to x_i * \hat B_tc
// ptv[k]: the kth time varying covariate parameter estimate
//         which is a scalar   

y1   = moptimize_util_depvar(M, 1)
ptco = moptimize_util_xb(M, b, 1)

// The following three lines are in the pgmhazgf2.ado file
// that runs this mata file.
// mata: st_view(v1=.,.,"`id'", "`touse'") 
// mata: st_view(v2=.,.,"`3'", "`touse'")
// mata: idinfo=panelsetup(v1,1)

external v1
external v2
moptimize_init_by(M,v1)
idinfo=panelsetup((*moptimize_util_by(M)),1)

// moptimize for each time varying covariate, put these in a matrix
ptv=J(cols(v2),1,NULL)
aptv=J(1,cols(v2),.)
for (i=1; i<=cols(v2); i++) {
j=i+1
ptv[i] = &moptimize_util_xb(M, b, j) 
aptv[1,i]=(*ptv[i])
}
plva = moptimize_util_xb(M, b, cols(v2)+2)

thetac=J(panelstats(idinfo)[1],1,.)
thetal=J(panelstats(idinfo)[1],1,.)
thetacij=exp(v2*aptv')
for (k=1; k<=panelstats(idinfo)[1]; k++) {
	// aptv is a vector of the time varying covariate parameter estimates
    tv=panelsubmatrix(thetacij,k,idinfo)
	thetac[k,1]=sum(tv, 0)
	thetal[k,1]=thetac[k,1]-tv[rows(tv),1]
}
Lc=(1:+exp(plva:+ptco):*thetac)
Ll=(1:+exp(plva:+ptco):*thetal)

fv  = y1:*ln(Ll:^(-exp(-plva))-Lc:^(-exp(-plva)))-(1:-y1):*exp(-plva):*ln(Lc)
// quadcolsum(fv)

if (todo==0) return
// gradient
mis = 0 
// covariates that are constant over time
gbu  = y1:*exp(ptco):*(	      Lc:^(-exp(-plva)-1) :*thetac ///
							 -Ll:^(-exp(-plva)-1) :*thetal	)  ///
						:/(	  Ll:^(-exp(-plva)) ///
							 -Lc:^(-exp(-plva))	) ///
-(1:-y1):*exp(ptco):*thetac:/Lc
Gb = gbu:*moptimize_util_indepvars(M, 1)

// gamma variance
gvu = y1:* ///
 		( ///
		 (exp(-plva):*ln(Ll)-exp(ptco):*thetal:/Ll):*Ll:^(-exp(-plva)) ///
		-(exp(-plva):*ln(Lc)-exp(ptco):*thetac:/Lc):*Lc:^(-exp(-plva)) ///
		) ///
		:/ (	   Ll:^(-exp(-plva)) ///
				  -Lc:^(-exp(-plva))	) ///
	+(1:-y1):*(exp(-plva):*ln(Lc) ///
	          -exp(ptco):*thetac:/Lc)
Gv = gvu:*moptimize_util_indepvars(M, 2)

// each time varying covariate
thetacgij=v2:*thetacij
thetacg=J(panelstats(idinfo)[1],cols(v2),.)
thetalg=J(panelstats(idinfo)[1],cols(v2),.)
Gtva=J(panelstats(idinfo)[1],cols(v2),.)
	for (i=1; i<=panelstats(idinfo)[1]; i++) {
		tv=panelsubmatrix(thetacgij,i,idinfo)
			thetacg[i,.]=quadcolsum(tv)
			thetalg[i,.]=thetacg[i,.]-tv[rows(tv),.]
	}	   

gtbu = y1:*exp(ptco):*(	      Lc:^(-exp(-plva)-1) :*thetacg ///
							 -Ll:^(-exp(-plva)-1) :*thetalg	)  ///
						:/(	  Ll:^(-exp(-plva)) ///
							 -Lc:^(-exp(-plva))	) ///
-(1:-y1):*exp(ptco):*thetacg:/Lc

for (tvi=1; tvi<=cols(v2); tvi++) {
j=tvi+2
Gtva[.,tvi] = gtbu[.,tvi]:*moptimize_util_indepvars(M,j)
}

S = Gb, Gtva, Gv
// quadcolsum(S) // This displays the gradient vector

if (todo==1) return
// Hessian
// time constant
h11=y1:*(gbu /// 
	     +(Lc:^(-exp(-plva)-2) :*(-exp(2:*ptco)-exp(plva:+2:*ptco)):*thetac:^2 ///
	      -Ll:^(-exp(-plva)-2) :*(-exp(2:*ptco)-exp(plva:+2:*ptco)):*thetal:^2  )  ///
	 :/(   Ll:^(-exp(-plva)) ///
		  -Lc:^(-exp(-plva))	   ) /// 
		- gbu:^2) ///  
 +(1:-y1):*(gbu - exp(plva:+ptco):*thetac:*gbu:/Lc  )
H[(1..cols(Gb)),(1..cols(Gb))] = moptimize_util_matsum(M, 1, 1, h11, mis)

// gamma variance
h22 =	y1:* ( ///
				-gvu ///
				+(  exp(plva:+2:*ptco):*thetal:^2 :*Ll:^(-exp(-plva)) :/ Ll:^2 ///
				   -exp(plva:+2:*ptco):*thetac:^2 :*Lc:^(-exp(-plva)) :/ Lc:^2 ) ///
			       :/(   Ll:^(-exp(-plva)) ///
				        -Lc:^(-exp(-plva))	   ) ///
			    + (	  (		exp(-plva):*ln(Ll)-exp(ptco):*thetal:/Ll ):^2 :*Ll:^(-exp(-plva)) ///
				     -(		exp(-plva):*ln(Lc)-exp(ptco):*thetac:/Lc ):^2 :*Lc:^(-exp(-plva))		) ///
			       :/(   Ll:^(-exp(-plva)) ///
				        -Lc:^(-exp(-plva))	   ) ///
				-gvu:^2
			  ) ///
		+(1:-y1):*(-gvu:+  thetac:^2 :*exp(plva:+2:*ptco):/Lc:^2)
H[(cols(Gb)+cols(v2)+1),(cols(Gb)+cols(v2)+1)] = moptimize_util_matsum(M, cols(v2)+1, cols(v2)+1, h22, mis)

// time varying
thetachij=v2:^2:*thetacij
thetach=J(panelstats(idinfo)[1],cols(v2),.)
thetalh=J(panelstats(idinfo)[1],cols(v2),.)
	for (i=1; i<=panelstats(idinfo)[1]; i++) {
		tv=panelsubmatrix(thetachij,i,idinfo)
		thetach[i,.]=quadcolsum(tv)
		thetalh[i,.]=thetach[i,.]-tv[rows(tv),.]
		}
htb = y1:*(((exp(2:*ptco)+exp(plva:+2:*ptco)):*Ll:^(-exp(-plva)-2) :*thetalg:^2-exp(ptco):*Ll:^(-exp(-plva)-1) :*thetalh ///
		   -(exp(2:*ptco)+exp(plva:+2:*ptco)):*Lc:^(-exp(-plva)-2) :*thetacg:^2+exp(ptco):*Lc:^(-exp(-plva)-1) :*thetach  ) ///
			:/ (Ll:^(-exp(-plva)) ///
			   -Lc:^(-exp(-plva)) ) ///
		  -gtbu:^2) ///
  +(1:-y1):*(-exp(ptco):*thetach:/Lc+gtbu:^2 :*exp(plva))
  for (tvi=1; tvi<=cols(v2); tvi++) {
j=tvi+1
k=tvi+cols(Gb)
H[k,k] = moptimize_util_matsum(M, j, j, htb[.,tvi], mis)
}

// gamma variance, time constant
h21 = y1:*(		///
		   ( ///
			 exp(plva:+2:*ptco):*thetal:^2 :*Ll:^(-exp(-plva)-2) - exp(-plva:+ptco):*ln(Ll):*Ll:^(-exp(-plva)-1) :*thetal + exp(2:*ptco):*thetal:^2 :*Ll:^(-exp(-plva)-2) ///
			-exp(plva:+2:*ptco):*thetac:^2 :*Lc:^(-exp(-plva)-2) + exp(-plva:+ptco):*ln(Lc):*Lc:^(-exp(-plva)-1) :*thetac - exp(2:*ptco):*thetac:^2 :*Lc:^(-exp(-plva)-2) ///
			) ///
			:/ (Ll:^(-exp(-plva)) ///
			   -Lc:^(-exp(-plva)) ) ///
		   -gvu:*gbu ///
		) ///
	+(1:-y1):*gbu:^2 :*exp(plva)
H[cols(Gb)+cols(v2)+1,(1..cols(Gb))] = moptimize_util_matsum(M, cols(v2)+2, 1, h21, mis) 

// time varying, time constant
h31 = y1:*(gtbu+exp(ptco):^2 :*((1+exp(plva)):*Ll:^(-exp(-plva)-2) :*thetal:*thetalg ///
							   -(1+exp(plva)):*Lc:^(-exp(-plva)-2) :*thetac:*thetacg ) ///
							:/ (Ll:^(-exp(-plva)) ///
							   -Lc:^(-exp(-plva)) ) ///
		  -gtbu:*gbu) ///
	+(1:-y1):*(gtbu+gbu:*gtbu:*exp(plva))
for (tvi=1; tvi<=cols(v2); tvi++) {
H[(cols(Gb)+tvi),(1..cols(Gb))] = moptimize_util_matsum(M, tvi+1, 1, h31[.,tvi], mis)
}

// gamma, time varying
h32  = y1:*(( (exp(plva)+1):*exp(2:*ptco):*thetal:*thetalg:*Ll:^(-exp(-plva)-2) - exp(-plva:+ptco):*ln(Ll):*Ll:^(-exp(-plva)-1) :*thetalg ///
			:-(exp(plva)+1):*exp(2:*ptco):*thetac:*thetacg:*Lc:^(-exp(-plva)-2) + exp(-plva:+ptco):*ln(Lc):*Lc:^(-exp(-plva)-1) :*thetacg  ) ///
			:/ (Ll:^(-exp(-plva)) ///
			   -Lc:^(-exp(-plva)) ) ///
			:-gvu:*gtbu) ///
	+(1:-y1):*(gtbu:*gbu:*exp(plva))
for (tvi=1; tvi<=cols(v2); tvi++) {
H[(cols(Gb)+cols(v2)+1),(cols(Gb)+tvi)] = moptimize_util_matsum(M, cols(v2)+2, tvi+1, h32[.,tvi], mis) 
}

// combinations of time varying
for (i=1; i<=cols(v2)-1; i++) {
if (i==1) {
j=2
thetach12ij=v2[.,i]:*v2[.,(2::cols(v2))]:*thetacij
}
else {
if (i<cols(v2)) {
j=j+1
thetach12ij=thetach12ij,v2[.,i]:*v2[.,(j::cols(v2))]:*thetacij
}
} 
}

thetach12=J(panelstats(idinfo)[1],cols(thetach12ij),.)
thetalh12=J(panelstats(idinfo)[1],cols(thetach12ij),.)
for (i=1; i<=panelstats(idinfo)[1]; i++) {
		thetach12[i,.]=quadcolsum(panelsubmatrix(thetach12ij,i,idinfo))
		thetalh12[i,.]=thetach12[i,.]-panelsubmatrix(thetach12ij,i,idinfo)[rows(panelsubmatrix(thetach12ij,i,idinfo)),.]
}

ct=0
for (tvj=1; tvj<=cols(v2)-1; tvj++) {
	for (tvi=2; tvi<=cols(v2); tvi++) {
	if (tvi>tvj) {
			ct=ct+1
			h32=y1:*(exp(ptco):*((1+exp(plva)):*exp(ptco):*Ll:^(-exp(-plva)-2) :*thetalg[.,tvj]:*thetalg[.,tvi] - Ll:^(-exp(-plva)-1) :*thetalh12[.,ct]  ///
								-(1+exp(plva)):*exp(ptco):*Lc:^(-exp(-plva)-2) :*thetacg[.,tvj]:*thetacg[.,tvi] + Lc:^(-exp(-plva)-1) :*thetach12[.,ct]  ) ///
								:/(Ll:^(-exp(-plva)) ///
							      -Lc:^(-exp(-plva)) ) ///
			- gtbu[.,tvj]:*gtbu[.,tvi]) ///
			+(1:-y1):*(	-exp(ptco):*thetach12[.,ct]:/Lc + gtbu[.,tvj]:*gtbu[.,tvi]:*exp(plva)	)
		H[(cols(Gb)+tvi),(cols(Gb)+tvj)] = moptimize_util_matsum(M, 2, tvi, h32, mis)
		}
	}
}
H=lowertriangle(H)+lowertriangle(H,0)'
}


void pgmhaz_gf2mataviewl(transmorphic M,
					real scalar todo, 
					real rowvector b, 
					real colvector fv, 
					real matrix S,
					real matrix H)
{

// y1=1 if an individual reports failure, 0 if right censored

// PARAMETERS TO BE ESTIMATED
// plva: log variance scalar
// ptco: time constance covariates column vector with
//       each row equal to x_i * \hat B_tc
// ptv[k]: the kth time varying covariate parameter estimate
//         which is a scalar   

y1   = moptimize_util_depvar(M, 1)
ptco = moptimize_util_xb(M, b, 1)

// The following three lines are in the pgmhazgf2.ado file
// that runs this mata file.
// mata: st_view(v1=.,.,"`id'", "`touse'") 
// mata: st_view(v2=.,.,"`3'", "`touse'")
// mata: idinfo=panelsetup(v1,1)

external v1
external v2

moptimize_init_by(M,v1)
idinfo=panelsetup((*moptimize_util_by(M)),1)

// moptimize for each time varying covariate, put these in a matrix
ptv=J(cols(v2),1,NULL)
aptv=J(1,cols(v2),.)

for (i=1; i<=cols(v2); i++) {
j=i+1
ptv[i] = &moptimize_util_xb(M, b, j) 
aptv[1,i]=(*ptv[i])
}

plva = moptimize_util_xb(M, b, cols(v2)+2)

thetac=J(panelstats(idinfo)[1],1,.)
thetal=J(panelstats(idinfo)[1],1,.)
thetacij=exp(v2*aptv')

for (k=1; k<=panelstats(idinfo)[1]; k++) {
	// aptv is a vector of the time varying covariate parameter estimates
    tv=panelsubmatrix(thetacij,k,idinfo)
	thetac[k,1]=sum(tv, 0)
	thetal[k,1]=thetac[k,1]-tv[rows(tv),1]
}

Lc=(1:+exp(plva:+ptco):*thetac)
Ll=(1:+exp(plva:+ptco):*thetal)

fv  = y1:*ln(Ll:^(-exp(-plva))-Lc:^(-exp(-plva)))-(1:-y1):*exp(-plva):*ln(Lc)
// quadcolsum(fv)

if (todo==0) return
// gradient
mis = 0 
// covariates that are constant over time
gbu  = y1:*exp(ptco):*(	      Lc:^(-exp(-plva)-1) :*thetac ///
							 -Ll:^(-exp(-plva)-1) :*thetal	)  ///
						:/(	  Ll:^(-exp(-plva)) ///
							 -Lc:^(-exp(-plva))	) ///
-(1:-y1):*exp(ptco):*thetac:/Lc
Gb = gbu:*moptimize_util_indepvars(M, 1)

// gamma variance
gvu = y1:* ///
 		( ///
		 (exp(-plva):*ln(Ll)-exp(ptco):*thetal:/Ll):*Ll:^(-exp(-plva)) ///
		-(exp(-plva):*ln(Lc)-exp(ptco):*thetac:/Lc):*Lc:^(-exp(-plva)) ///
		) ///
		:/ (	   Ll:^(-exp(-plva)) ///
				  -Lc:^(-exp(-plva))	) ///
	+(1:-y1):*(exp(-plva):*ln(Lc) ///
	          -exp(ptco):*thetac:/Lc)
Gv = gvu:*moptimize_util_indepvars(M, 2)

// each time varying covariate
thetacgij=v2:*thetacij
thetacg=J(panelstats(idinfo)[1],cols(v2),.)
thetalg=J(panelstats(idinfo)[1],cols(v2),.)
Gtva=J(panelstats(idinfo)[1],cols(v2),.)
	for (i=1; i<=panelstats(idinfo)[1]; i++) {
		tv=panelsubmatrix(thetacgij,i,idinfo)
			thetacg[i,.]=quadcolsum(tv)
			thetalg[i,.]=thetacg[i,.]-tv[rows(tv),.]
	}	   

gtbu = y1:*exp(ptco):*(	      Lc:^(-exp(-plva)-1) :*thetacg ///
							 -Ll:^(-exp(-plva)-1) :*thetalg	)  ///
						:/(	  Ll:^(-exp(-plva)) ///
							 -Lc:^(-exp(-plva))	) ///
-(1:-y1):*exp(ptco):*thetacg:/Lc

for (tvi=1; tvi<=cols(v2); tvi++) {
j=tvi+2
Gtva[.,tvi] = gtbu[.,tvi]:*moptimize_util_indepvars(M,j)
}

S = Gb, Gtva, Gv
// quadcolsum(S) // This displays the gradient vector

if (todo==1) return
// Hessian
// time constant
h11=y1:*(gbu /// 
	     +(Lc:^(-exp(-plva)-2) :*(-exp(2:*ptco)-exp(plva:+2:*ptco)):*thetac:^2 ///
	      -Ll:^(-exp(-plva)-2) :*(-exp(2:*ptco)-exp(plva:+2:*ptco)):*thetal:^2  )  ///
	 :/(   Ll:^(-exp(-plva)) ///
		  -Lc:^(-exp(-plva))	   ) /// 
		- gbu:^2) ///  
 +(1:-y1):*(gbu - exp(plva:+ptco):*thetac:*gbu:/Lc  )
H[(1..cols(Gb)),(1..cols(Gb))] = moptimize_util_matsum(M, 1, 1, h11, mis)

// gamma variance
h22 =	y1:* ( ///
				-gvu ///
				+(  exp(plva:+2:*ptco):*thetal:^2 :*Ll:^(-exp(-plva)) :/ Ll:^2 ///
				   -exp(plva:+2:*ptco):*thetac:^2 :*Lc:^(-exp(-plva)) :/ Lc:^2 ) ///
			       :/(   Ll:^(-exp(-plva)) ///
				        -Lc:^(-exp(-plva))	   ) ///
			    + (	  (		exp(-plva):*ln(Ll)-exp(ptco):*thetal:/Ll ):^2 :*Ll:^(-exp(-plva)) ///
				     -(		exp(-plva):*ln(Lc)-exp(ptco):*thetac:/Lc ):^2 :*Lc:^(-exp(-plva))		) ///
			       :/(   Ll:^(-exp(-plva)) ///
				        -Lc:^(-exp(-plva))	   ) ///
				-gvu:^2
			  ) ///
		+(1:-y1):*(-gvu:+  thetac:^2 :*exp(plva:+2:*ptco):/Lc:^2)
H[(cols(Gb)+cols(v2)+1),(cols(Gb)+cols(v2)+1)] = moptimize_util_matsum(M, cols(v2)+1, cols(v2)+1, h22, mis)

// time varying
thetachij=v2:^2:*thetacij
thetach=J(panelstats(idinfo)[1],cols(v2),.)
thetalh=J(panelstats(idinfo)[1],cols(v2),.)
	for (i=1; i<=panelstats(idinfo)[1]; i++) {
		tv=panelsubmatrix(thetachij,i,idinfo)
		thetach[i,.]=quadcolsum(tv)
		thetalh[i,.]=thetach[i,.]-tv[rows(tv),.]
		}
htb = y1:*(((exp(2:*ptco)+exp(plva:+2:*ptco)):*Ll:^(-exp(-plva)-2) :*thetalg:^2-exp(ptco):*Ll:^(-exp(-plva)-1) :*thetalh ///
		   -(exp(2:*ptco)+exp(plva:+2:*ptco)):*Lc:^(-exp(-plva)-2) :*thetacg:^2+exp(ptco):*Lc:^(-exp(-plva)-1) :*thetach  ) ///
			:/ (Ll:^(-exp(-plva)) ///
			   -Lc:^(-exp(-plva)) ) ///
		  -gtbu:^2) ///
  +(1:-y1):*(-exp(ptco):*thetach:/Lc+gtbu:^2 :*exp(plva))
  for (tvi=1; tvi<=cols(v2); tvi++) {
j=tvi+1
k=tvi+cols(Gb)
H[k,k] = moptimize_util_matsum(M, j, j, htb[.,tvi], mis)
}

// gamma variance, time constant
h21 = y1:*(		///
		   ( ///
			 exp(plva:+2:*ptco):*thetal:^2 :*Ll:^(-exp(-plva)-2) - exp(-plva:+ptco):*ln(Ll):*Ll:^(-exp(-plva)-1) :*thetal + exp(2:*ptco):*thetal:^2 :*Ll:^(-exp(-plva)-2) ///
			-exp(plva:+2:*ptco):*thetac:^2 :*Lc:^(-exp(-plva)-2) + exp(-plva:+ptco):*ln(Lc):*Lc:^(-exp(-plva)-1) :*thetac - exp(2:*ptco):*thetac:^2 :*Lc:^(-exp(-plva)-2) ///
			) ///
			:/ (Ll:^(-exp(-plva)) ///
			   -Lc:^(-exp(-plva)) ) ///
		   -gvu:*gbu ///
		) ///
	+(1:-y1):*gbu:^2 :*exp(plva)
H[cols(Gb)+cols(v2)+1,(1..cols(Gb))] = moptimize_util_matsum(M, cols(v2)+2, 1, h21, mis) 

// time varying, time constant
h31 = y1:*(gtbu+exp(ptco):^2 :*((1+exp(plva)):*Ll:^(-exp(-plva)-2) :*thetal:*thetalg ///
							   -(1+exp(plva)):*Lc:^(-exp(-plva)-2) :*thetac:*thetacg ) ///
							:/ (Ll:^(-exp(-plva)) ///
							   -Lc:^(-exp(-plva)) ) ///
		  -gtbu:*gbu) ///
	+(1:-y1):*(gtbu+gbu:*gtbu:*exp(plva))
for (tvi=1; tvi<=cols(v2); tvi++) {
H[(cols(Gb)+tvi),(1..cols(Gb))] = moptimize_util_matsum(M, tvi+1, 1, h31[.,tvi], mis)
}

// gamma, time varying
h32  = y1:*(( (exp(plva)+1):*exp(2:*ptco):*thetal:*thetalg:*Ll:^(-exp(-plva)-2) - exp(-plva:+ptco):*ln(Ll):*Ll:^(-exp(-plva)-1) :*thetalg ///
			:-(exp(plva)+1):*exp(2:*ptco):*thetac:*thetacg:*Lc:^(-exp(-plva)-2) + exp(-plva:+ptco):*ln(Lc):*Lc:^(-exp(-plva)-1) :*thetacg  ) ///
			:/ (Ll:^(-exp(-plva)) ///
			   -Lc:^(-exp(-plva)) ) ///
			:-gvu:*gtbu) ///
	+(1:-y1):*(gtbu:*gbu:*exp(plva))
for (tvi=1; tvi<=cols(v2); tvi++) {
H[(cols(Gb)+cols(v2)+1),(cols(Gb)+tvi)] = moptimize_util_matsum(M, cols(v2)+2, tvi+1, h32[.,tvi], mis) 
}

// combinations of time varying
for (tvj=1; tvj<cols(v2); tvj++) {
	
		thetach12ij=(v2)[.,tvj]:*(v2)[.,tvj+1::cols(v2)]:*thetacij

		thetach12=J(panelstats(idinfo)[1],cols(thetach12ij),.)
		thetalh12=J(panelstats(idinfo)[1],cols(thetach12ij),.)
		for (i=1; i<=panelstats(idinfo)[1]; i++) {
				thetach12[i,.]=quadcolsum(panelsubmatrix(thetach12ij,i,idinfo))
				thetalh12[i,.]=thetach12[i,.]-panelsubmatrix(thetach12ij,i,idinfo)[rows(panelsubmatrix(thetach12ij,i,idinfo)),.]
		}
	
	ct=0
	for (tvi=tvj+1; tvi<=cols(v2); tvi++) {
			ct=ct+1
			h32=y1:*(exp(ptco):*((1+exp(plva)):*exp(ptco):*Ll:^(-exp(-plva)-2) :*thetalg[.,tvj]:*thetalg[.,tvi] - Ll:^(-exp(-plva)-1) :*thetalh12[.,ct]  ///
								-(1+exp(plva)):*exp(ptco):*Lc:^(-exp(-plva)-2) :*thetacg[.,tvj]:*thetacg[.,tvi] + Lc:^(-exp(-plva)-1) :*thetach12[.,ct]  ) ///
								:/(Ll:^(-exp(-plva)) ///
							      -Lc:^(-exp(-plva)) ) ///
			- gtbu[.,tvj]:*gtbu[.,tvi]) ///
			+(1:-y1):*(	-exp(ptco):*thetach12[.,ct]:/Lc + gtbu[.,tvj]:*gtbu[.,tvi]:*exp(plva)	)
		H[(cols(Gb)+tvi),(cols(Gb)+tvj)] = moptimize_util_matsum(M, 2, tvi, h32, mis)
	}
}
H=lowertriangle(H)+lowertriangle(H,0)'
}

mata mosave pgmhaz_gf2mataviewo(), dir(C:\ado\plus\p) replace 
mata mosave pgmhaz_gf2mataviewl(), dir(C:\ado\plus\p) replace 
end
