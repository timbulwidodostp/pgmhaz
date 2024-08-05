*! pgmhazgf2 v1.0.3 Michael G. Farnworth, 16 January 2012

capture program drop pgmhazgf2
program define pgmhazgf2, sortpreserve byable(onecall) properties(svylb)

// -svylb- is and undocumented property that is used with -gf- programs.
// A similar document property is -svyr- and information on this is in the
// Stata Programing Reference Manual, version 11.
//
// The svyr property indicates to the svy prefix command that the requirements
// to use the vce(linearized) variance estimation method has been satisfied.

version 11.1
tempname pre_est
capture _return hold `pre_est'

if _by() {
	by `_byvars' `_byrc0': Estimate `0'
} 
else {
	if replay() {
			if "`e(cmd)'" != "pgmhazgf2" {
				error 301
				}
			if _by() {
				 error 190
				} 
			Replay `0'	
		}
	else {
		Estimate `0'
		}
}

_return restore `pre_est'
end

capture program drop Estimate
program define Estimate, eclass byable(recall)

	syntax anything [if] [in] [fweight pweight iweight] [,  ///
		Id(varname) Seq(varname numeric) LNVar0(real -1) noLOg ///
		noCONStant noBeta0 Level(cilevel) EForm  USEMYB USEMYBC LARGENTV * ]
		
	ereturn clear
		
	local title1 "Complementary log-log model"

	local title "Prentice-Gloeckler-Meyer hazard"
		
	mlopts mlopts, `options'	

	if "`weight'" != "" {
		local wgt "[`weight'`exp']"
	}
	
	if "`2'" == "" | "`2'" == " " | "`2'" == "  " {
		di as error ""
		di as error "At least one variable must be listed inside the second set of quotes."
		di as error "This can be a variable that is always equal to one along with the"
		di as error "noconstant option."
		di as error ""
		exit 198
	}	
	
	tempname check3
	fvexpand `3'
	local counttv=0
	foreach i in `r(varlist)' {
	 local counttv=`counttv'+1
	 }
	 if `counttv'<2 {
	 	di as error ""
		di as error "At least two variables must be listed inside the third set of quotes."
		di as error "Note: These can be time constant variables."
		di as error ""
		exit 198
	 }

		if "`id'" == "" {
		di as error ""
		di as error "A variable identifying each person must "
		di as error " be specified in id(idvar) "
		exit 198
			}
			
		if "`seq'" == "" {
		di as error ""
		di as error "An integer-valued variable identifying "
		di as error "each discret time period must be specified "
		di as error "in seq(seqvar) "
		exit 198
	}
		
	capture assert "`3'"=="" | "`3'"==" " | "`3'"=="   "
			if _rc==0 {
				di as error ""
		di as error "Three sets of variables must be listed."
		di as error ""
		di as error "The first is the dependent variable within double quotes,"
		di as error "the second is a list of time constant covariates within"
		di as error "double quotes, and the third is a list of time varying covariates"
		di as error "within double quotes."
		di as error ""
		di as error "Note: If there are no time varying covariates then one time constant"
		di as error "covariate can be listed within the third pair of quotes."
		di as error ""
		di as error "Note: This program requires at least two covariates and at least one"
		di as error "of them must be constant over time each for individual. When the noconstant"
		di as error "option is used a variable that is always equal to one can be included as a"
		di as error "covariate."
			}
		
	marksample touse
	markout `touse' `1' `2' `3'
	local log2 = cond("`beta0'" == "", "noisily", "quietly")
	
	sort  `id' `touse' `seq'
	
	mata: st_view(v1=.,.,"`id'", "`touse'") 
	mata: st_view(v2=.,.,"`3'", "`touse'")
	mata: idinfo=panelsetup(v1,1)
	
	display ""
	tempname number_of_subjects
	mata: st_numscalar("`number_of_subjects'",panelstats(idinfo)[1])
	display "number of subjects = " scalar(`number_of_subjects')
	quietly summ `1' if `1'>0 & `touse'
	display "number of subjects who report failure = " r(N)
	tempname minimum_timeperiods
	mata: st_numscalar("`minimum_timeperiods'",panelstats(idinfo)[3])
	display "minimum number of time periods that a subject is observed = " scalar(`minimum_timeperiods')
	tempname maximum_timeperiods
	mata: st_numscalar("`maximum_timeperiods'",panelstats(idinfo)[4])
	display "maximum number of time periods that a subject is observed = " scalar(`maximum_timeperiods')
	tempname total_timeperiods
	mata: st_numscalar("`total_timeperiods'",panelstats(idinfo)[2])
	display "number of time periods, summed across subjects = " scalar(`total_timeperiods')
	
	tempname b0z b0c
	mata: st_view(tcvars=.,1,"`2'", "`touse'")
	if "`constant'" ~= "" {
		mata: b0z=J(1,cols(tcvars)+cols(v2)+1,0)
	}
	else {
		mata: b0z=J(1,cols(tcvars)+cols(v2),0)
	}
	mata: st_matrix("`b0z'",b0z)
	mata: mata drop b0z
	matrix `b0c'=`b0z'
	
	local usemyb0c=cond("`usemybc'"=="","*","")
	`usemyb0c' capture matrix `b0c'=myb0c
	
	tempname checkdiffc 
	scalar checkdiffc=mreldif(`b0c',`b0z')
	
// Put a time varying covariate name on
// each time varying equation. Also identify
// the number of time varying covariates for below.
fvexpand `3'
tempname mltv
scalar `mltv'=r(varlist)
mata: mltv=st_strscalar("`mltv'")
mata: mltv=tokens(mltv)
mata: ctvc=cols(mltv)
mata: mltv=subinstr(mltv, ".", "_")
mata: mltv="(":+mltv:+": )"
mata: mltv=invtokens(mltv)
mata: st_local("mltv", mltv)
mata: st_numscalar("ctvc",ctvc)
mata: mata drop ctvc
mata: mata drop mltv 
local ctvc=scalar(ctvc)
	
	if checkdiffc~=0 {	
	display as text ""
	display as text "compelmentary log-log starting values: matrix myb0c"
	`log2' cloglog `1' `2' `3' `wgt' if `touse', `constant' level(`level') `eform' `log' from(`b0c', copy)
	}
	else {
	`log2' cloglog `1' `2' `3' `wgt' if `touse', `constant' level(`level') `eform' `log' 
	}
	

	local LL1 = e(ll)
	
	capture assert `1'==1 | `1'==0 if e(sample)
		if _rc~=0 {
		di as error ""
		di as error "The dependent variable must equal 0 or 1. Zero identifies"
		di as error "descrite time periods that an individual survives and a 1"
		di as error "identifies a discrete time period within which a person fails."
		exit 198
		}

	* starting values
	tempname b0
	mata: b0=st_matrix("e(b)")
	if "`constant'" == "" {
		mata: b0=b0[|1\cols(tcvars)|],b0[cols(tcvars)+cols(v2)+1],b0[|cols(tcvars)+1\cols(tcvars)+cols(v2)|]
	}
	mata: st_matrix("`b0'",b0)
	mata: mata drop tcvars b0
	matrix `b0'=`b0',`lnvar0'
	tempname b0check
	matrix `b0check'=`b0'
	
local usemyb0=cond("`usemyb'"=="","*","")	
capture `usemyb0' matrix `b0'=myb0

tempname checkdiff 
scalar checkdiff=mreldif(`b0',`b0check')

if checkdiff==0 {
display as text ""
display as text "pgmhaz starting values: complementary log-log parameter estimates"
display as text "                        and /ln_gvar=`lnvar0'"
}
else {
display as text ""
display as text "pgmhaz starting values: matrix myb0"
}

	global S_E_id "`id'"
	
	global tvv "`3'"
	
	tempname last_obs
	quietly generate byte `last_obs'=0
	quietly by `id': replace `last_obs'=1 if _n==_N

tempname pgmhazo
scalar pgmhazo=cond("`largentv'"=="",1,2)
if pgmhazo==1 {
		ml model gf2 pgmhaz_gf2mataviewo() (`1': `1'= `2', `constant' ) `mltv'  /ln_gvar `wgt' ///
		if `touse' & `last_obs'==1, `log'  ml maximize search(off) title(`title') ///
		nopreserve missing waldtest(0) init(`b0', copy) `mlopts'
		}
else {
		display ""
		display "Note: The largentv option has been employed"
		display ""
		ml model gf2 pgmhaz_gf2mataviewl() (`1': `1'= `2', `constant' ) `mltv'  /ln_gvar `wgt' ///
		if `touse' & `last_obs'==1, `log'  ml maximize search(off) title(`title') ///
		nopreserve missing waldtest(0) init(`b0', copy) `mlopts'
}
	
	ereturn local cmd pgmhazgf2

		tempname tag
		egen `tag' = tag(`id') if e(sample) 
		quietly count if `tag' & e(sample)
		ereturn scalar N_spell = r(N)
		ereturn scalar N = `total_timeperiods'

		ereturn local cmd "pgmhazgf2"		
		ereturn local depvar "`1'"

		ereturn local idvar "`id'"
		ereturn local seqvar "`seq'"

	* display ancillary parameters all together
	local ancl=`ctvc'+1
	ereturn scalar k_aux = `ancl' 

			ereturn scalar ll_nofr = `LL1'

                // v = exp(lnvarg); s.e.(v) = v * s.e.(lnvarg) 
	        ereturn scalar gammav = exp( [ln_gvar]_cons ) 
	        ereturn scalar se_gammav = `e(gammav)' * [ln_gvar]_se[_cons] 
			ereturn scalar lltest = -2 * ( `e(ll_nofr)' - `e(ll)' )
			ereturn scalar lltest_p = .5*chiprob(1,`e(lltest)')
			
	Replay, level(`level') `eform'
end

capture program drop Replay
program define Replay,  

	syntax [, Level(integer `c(level)') EForm * ]

	local diopts "`options'"
	ml display, level(`level') `eform' plus
	if `level' < 10 | `level' > 99 {
		local level = 95
		}
	 
			
		if "`eform'" == "eform" {
			foreach i of global tvv  {
				_diparm `i', exp level(`level') label(exp(_b[`i':_cons])) prob 
			}
			
		}
		_diparm ln_gvar, exp level(`level') label("Gamma var.") prob
		display as text "{hline 13}{c BT}{hline 64}"
		display as text "LR test of Gamma var. = 0: {help j_chibar:chibar2(01) = } " as res %8.6g `e(lltest)' _c
		display as text "  Prob.>=chibar2 = " as res  %8.6g  `e(lltest_p)'
		display as text " "

end
