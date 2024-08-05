*! version 1.2.1 Stephen P. Jenkins, July 1997  STB-39 sbe17
*! Prentice-Gloeckler-Meyer hazard model with Gamma unobserved heterogeneity
*! Syntax: pgmhaz <varlist> [if <exp>] [in <range>], id(idvar)
*!      dead(deadvar) seq(seqvar) lnvar0(#) [trace eform level(#) nolog nocons]

program define pgmhaz
	version 5.0
/*------------------------------------------------ playback request */
	local options "LEvel(integer $S_level) EForm"
	if "`*'" == "" | substr("`1'",1,1) == "," {
		if "$S_E_cmd" ~= "pgmhaz" {
			error 301
		}
		parse "`*'"
	}
/*------------------------------------------------ estimation */
	else {
		local varlist "req ex"
		local if "opt"
		local in "opt"
		local options "Id(string) Dead(string) Seq(string) TRace D1 D2"
		local options "`options' LNVar0(real -1) noLOG NOCONS"
		local options "`options'  LEvel(integer $S_level) EForm"
		parse "`*'"

		if "`id'" == "" {
			di in red "Variable identifying each person must "
			di in red " be specified in id(idvar) "
			exit 198
		}
		confirm variable `id'
		local n : word count `id'
		if `n' ~= 1 {
			di in red "id(idvar) must contain only one variable"
			exit 198
		}

		if "`dead'" == "" {
			di in red "Per-interval censoring indicator variable must "
			di in red "be specified in dead(deadvar):"
			di in green "see -pgmhaz- help file"
			exit 198
		}
		confirm variable `dead'
		local n : word count `dead'
		if `n' ~= 1 {
			di in red "dead(.) must contain only one variable"
			exit 198
		}

		capture assert `dead'==1 | `dead'==0
		if _rc~=0 {
			di in red "Per-interval censoring indicator variable must "
			di in red "equal one or zero: see -pgmhaz- help file"
			exit 198	
		}

		if "`seq'" == "" {
			di in red "Integer-valued variable identifying "
			di in red "spell time interval be specified " _c
			di in red "in seq(seqvar) "
			exit 198
		}
		confirm variable `seq'
		local n : word count `seq'
		if `n' ~= 1 {
			di in red "seq(seqvar) must contain only one variable"
			exit 198
		}

		if "`log'" != "" {
			local prel "quietly "
		}

		if "`nocons'" != "" {
			local noconst "constant(01)"
		}

		parse "`varlist'", parse(" ")

			/* identify estimation sample */
		tempvar touse
		mark `touse' `if' `in'
		markout `touse' `varlist' `id' `dead' `seq'

		tempvar mysamp 
		tempname b b0 b1 b2 lnf V lnvar

		set more 1

			/* estimate intercept-only model and save logL */

		qui glm `dead'  if `touse', f(b) l(c) 
		local LL0 = -0.5 * $S_3
		local rdf0 = $S_2
		di " "

			/* get estimates of no-het model */

		di in gr "(1) PGM hazard model without " _c
		di in gr "unobserved heterogeneity"
		glm `dead' `varlist' if `touse', f(b) l(c) /*
		*/ level(`level') `eform' `log' `nocons'

		local LL1 = -0.5 * $S_3
		local df = `rdf0' - $S_2
		local chisq = -2 * (`LL0' - `LL1')
		di in gr "Log likelihood (-0.5*Deviance) = " in ye `LL1' 
		di in gr "   Cf. log likelihood for intercept-only model (Model 0) = " /*
		*/ in ye `LL0' 
		di in gr "   Chi-squared statistic for Model (1) vs. Model (0) = " in ye `chisq'
		di in gr "   Prob. > chi2(`df') = " in ye chiprob(`df',`chisq')
		di " "

		di in gr "(2) PGM hazard model with " _c
		di in gr "Gamma distributed unobserved heterogeneity"
		di " "
			/* get coeff starting values */
		matrix `b0' = get(_b)
		matrix coleq `b0' = hazard

			/* now pack out -b0- with starting
			  value for lnvarg  */

		matrix `lnvar' = (`lnvar0')
		matrix colnames `lnvar' = ln_varg:_cons
		matrix `b1' = `b0',`lnvar'

			/* estimate full model */
		
		sort `id' `touse'  `seq'

		eq hazard : `dead' `varlist'
		eq ln_varg : 

	/* stuff below in terms of lnvarg = ln(varg) rather than
	   varg itself, in order to constrain varg values to be >0  */

		global S_E_id "`id'"
		global S_E_dd "`dead'"
		
		ml begin
		ml function pgm_ll

		ml method deriv0
		ml model `b' = hazard ln_varg, depv(10) from(`b1') `noconst'
			/* starting values for regressor coeffs are
			   taken from -glm- estimates, and 
			   lnvarg is  set equal to minus one by default
			   i.e. varg = exp(lnvarg) = exp(-1) ~= .37 */		
		ml sample `mysamp' if `touse'
		`prel' ml maximize `lnf' `V' , `trace' `d1' `d2'

		ml post pgmhaz,  /*
		*/ title("PGM hazard model with Gamma heterogeneity")

		global S_1 = $S_E_ll
		global S_2 = `LL1'
		global S_E_cmd = "pgmhaz"

	}

	if `level' < 10 | `level' > 99 {
		local level = 95
	}

	if "`eform'" ~= "" {
		local eform eform(`dead':)
	}

	ml mlout pgmhaz, level(`level') `eform'

		/* v = exp(lnvarg); s.e.(v) = v * s.e.(lnvarg) */
	local v = exp( [ln_varg]_coef[_cons] )
	local sev = `v' * [ln_varg]_se[_cons] 
	di in gr "Gamma variance, exp(ln_varg) = " in ye `v'  _c
	di in gr "; Std. Err. = " in ye `sev' "; z = " in ye `v'/`sev'

	tempname lltest
	scalar `lltest' = -2 * ($S_2 - $S_1)
	di " "
	di in gr  "Likelihood ratio statistic for testing " _c
	di in gr "models (1) vs (2) = " in ye `lltest'
	di in gr "Prob. test statistic > chi2(1) = " in ye /*
		*/ chiprob(1,`lltest')

end



