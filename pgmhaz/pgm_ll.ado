*! version 1.2.1 Stephen P. Jenkins July 1997  STB-39 sbe17

program define pgm_ll

	version 5.0
	local b "`1'" 	 /* parameter values */
	local lnf "`2'"  /* scalar name to contain lnL fn. value */

	tempvar I sum sum1 sum2 lnfi

	tempname bcoef v lnvarg

		/* create Index scores, i.e. X_itj*b_j */

	local k = colsof(`b') - 1
	matrix `bcoef' = `b'[1,1..`k']
	matrix score double `I' = `bcoef' if $S_mlwgt

	scalar `lnvarg' = `b'[1,`k'+1]

	capture assert `I' + `lnvarg' < 200 if $S_mlwgt
	if _rc {
		di in red "exponential overflow"
		di in green "I_ij + ln(v) >= 200"
		exit 1400
	}

	qui by $S_E_id: gen double `sum' = sum(exp(`I' + `lnvarg')) /*
		*/ if $S_mlwgt

	qui by $S_E_id: gen double `sum1' = `sum'[_N] if _n==_N
	qui by $S_E_id: gen double `sum2' = `sum'[_N-1] if _n==_N
	qui by $S_E_id: replace `sum2' = 0 if _N == 1 & $S_mlwgt

	scalar `v' = exp(`lnvarg')
	
	ge double `lnfi' = sum( cond( $S_E_dd, /*
	  */  ln( (1+`sum2')^(-1/`v') - (1+`sum1')^(-1/`v') ), /*
	  */  -ln(1+`sum1')/`v'   ))

	scalar `lnf' = `lnfi'[_N]

end

