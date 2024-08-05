{smcl}
{* *! version 1.1.0 19jan2011}{...}
{cmd:help pgmhazgf2}{right: ({browse "http://www.stata-journal.com/article.html?article=st0256":SJ12-2: st0256})}
{hline}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{hi:pgmhazgf2} {hline 2}}Discrete time (grouped data) proportional hazards models{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 17 2}
{cmd:pgmhazgf2} {cmd:"}{it:censoring_status_variable}{cmd:"}
{cmd:"}{it:time_constant_covariates}{cmd:"}
{cmd:"}{it:time_varying_covariates}{cmd:"}
{ifin}
{weight}{cmd:,}
	{cmdab:i:d(}{it:idvar}{cmd:)}  
	{cmdab:s:eq(}{it:seqvar}{cmd:)}
	[{it:options}]

{p 8 17 2}{cmd:pgmhazgf2} [{cmd:,} {cmdab:ef:orm}
{cmdab:l:evel(}{it:#}{cmd:)}] 

{pstd}
When {cmd:pgmhazgf2} is the last estimation command used, typing the
{cmd:pgmhazgf2} command without arguments redisplays the previous estimates.

{pstd}
Note: {cmd:pgmhazgf2} requires the file called {cmd:pgmhazgf2.ado}, as
well as two files called {cmd:pgmhaz_gf2mataviewo.mo} and
{cmd:pgmhaz_gf2mataviewl.mo}.  Typing {cmd:do pgmhaz_gf2mataview.mata} creates
these two {cmd:.mo} files in the path "C:\ado\plus\p".  If this path does not
exist, it should be created before running this Mata file.  If you prefer 
a different path, you may change the {it:path} in the command line
{cmd:mata mosave, dir(}{it:path}{cmd:)} within the
{cmd:pgmhaz_gf2mataview.mata} file.

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{p2coldent:* {cmdab:i:d(}{it:idvar}{cmd:)}}specify the variable uniquely identifying each subject{p_end}
{p2coldent:* {cmdab:s:eq(}{it:seqvar}{cmd:)}}specify the variable uniquely identifying each time interval at risk for each subject{p_end}
{synopt:{opt largentv}}change how the Hessian matrix is constructed {p_end}
{synopt:{cmdab:lnv:ar0(}{it:#}{cmd:)}}log-gamma variance starting value {p_end}
{synopt:{cmd:usemyb}}starting values for the Prentice-Gloeckler-Meyer hazard
regression with gamma frailty {p_end}
{synopt:{cmd:usemybc}}starting values for the {helpb cloglog} model {p_end}
{synopt:{cmdab:ef:orm}}report coefficients transformed to hazard ratio format {p_end}
{synopt:{cmdab:nocons:tant}}suppress constant term{p_end}
{synopt:{cmdab:nob:eta0}}suppress reporting {helpb cloglog} estimates {p_end}
{synopt:{cmdab:l:evel(}{it:#}{cmd:)}}significance level for confidence intervals {p_end}
{synopt:{it:maximize_options}}control the maximization process {p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* {opt id(idvar)} and {opt seq(seqvar)} are required.{p_end}
{p 4 4 2}{cmd:svy} and {cmd:by} prefixes are allowed; see {help prefix}.{p_end}
{p 4 6 2}
{cmd:fweight}s, {cmd:iweight}s, and {cmd:pweight}s  are allowed; see {help weight}.{p_end}
{p 4 6 2}{helpb predict} can be used after estimation.  With the
{cmd:score} option, an {helpb if} condition is required so that just one
prediction is made for each subject.{p_end}
{p 4 6 2}Currently, when the {helpb suest} command is used after
estimation, the following error messages appear: "unable to generate
scores" and "{cmd:suest} requires that {cmd:predict} allow the {cmd:score}
option".{p_end}
{p 4 6 2}Indicator variables can be applied with factor variables;
see {findalias frfvvarlists}.  Not all the factor-variable features
documented can be used under the {cmd:pgmhazgf2} command.


{title:Description}

{pstd}{cmd:pgmhazgf2} fits by maximum likelihood two discrete-time
(grouped data) proportional hazards regression models that are used to
examine single-spell data with no left-truncated survival times.  One of
these models incorporates a gamma mixture distribution to summarize
unobserved individual heterogeneity (or frailty).  Covariates may
include regressor variables summarizing observed differences between
persons (either fixed or time varying) and variables summarizing the
duration dependence of the hazard rate.  With suitable definition of
covariates, models with a fully nonparametric specification for duration
dependence may be estimated; so too may parametric specifications.
{cmd:pgmhazgf2} thus provides a useful complement to {helpb stcox} (for
continuous survival time data) and related programs.  For further
discussion of hazard regression models with unobserved heterogeneity,
see Jenkins (2008), especially {it:Lesson 7}.  This lesson also shows
how to derive predicted hazard and survivor functions from the estimates
of the model.

{pstd}{cmd:pgmhazgf2} fits two models by maximum likelihood: 1) the
Prentice and Gloeckler (1978) {helpb cloglog} model; and 2) the Prentice
and Gloeckler (1978) model incorporating a gamma mixture distribution to
summarize unobserved individual heterogeneity, as proposed by Meyer
(1990).  Specifically, the Prentice-Gloeckler-Meyer models fit are
those described by Meyer (1990, eq. 5 and 7).

{pstd}For suitably reorganized data, model 1 can be fit with the
{helpb cloglog} command; see Allison (1982) or Jenkins (1995).  To
fit the discrete-time proportional hazards model with normally
(rather than gamma) distributed heterogeneity, use {helpb xtcloglog}
with the data organized as here.  To fit a model with a discrete-mixture
distribution that summarizes subject heterogeneity as proposed by Heckman and
Singer (1984), see Jenkins (2004a).

{pstd}Model 2 is fit using {helpb ml} and a {helpb Mata}
{cmd:gf2}, with starting values for the coefficients taken from Model
1's estimates.  The program estimates the log of the gamma variance,
with a default value equal to -1, that is, a variance of approximately
0.37.  Different starting values may be set optionally; see below.


{title:Important note about data organization and variables}

{pstd}The dataset must be organized beforehand so that a data row
corresponds to each time interval at risk of the event for each subject.
This corresponds to "episode splitting" at the end of every time
interval, and {helpb expand} is useful for putting the data in this
form.  Also see the "data step" discussion in Jenkins (1995).

{pstd}With regard to this command, the censoring status variable
indicates censoring during each time interval at risk.  If a person is
censored, then the censoring status variable equals 0 for all
j=1,2,...,k_i.  If a person is not censored then this variable equals
zero for all j=1,2,3,...,k_i-1, and it is equal to 1 for j=k_i.

{pstd}For this command, three sets of variables must be separately
listed.  The first is the censoring status variable within double
quotes; the second is a list of time-constant covariates within double
quotes; and the third is a list of time-varying covariates within double
quotes.  This program requires at least two time-varying covariates.  If
there are less than two, then time-constant covariates can be listed
within the third quotes.  This program requires at least three
covariates, and at least one of them must be constant over time for each
individual.  When the {cmd:noconstant} option is used, a variable that
is always equal to 1 can be included.


{title:Options}

{phang}{cmd:id(}{it:idvar}{cmd:)} specifies the variable uniquely
identifying each subject, i.  {cmd:id()} is required.

{phang}{cmd:seq(}{it:seqvar}{cmd:)} is the variable uniquely identifying
each time interval at risk for each subject.  For each i, the variable
is the integer sequence 1,2,...,j.  {cmd:seq()} is required.

{phang}{opt largentv} changes how the Hessian matrix is constructed and
under some circumstances speeds up the estimation process.  Without this
option, a matrix is constructed that has a row for each observation
examined and a column for each possible combination of time-varying
covariates.  For example, with two time-varying covariates, there will
be 1 column and with five time-varying covariates, there will be 10
(4+3+2+1).  With a small number of time-varying covariates, this option
slows the maximization process down.  With a large number of
time-varying covariates, this option may speed up the maximization
process by working with smaller matrices rather than one large one.

{phang}{cmd:lnvar0(}{it:#}{cmd:)} specifies the value for the log of the
gamma variance that is used as the starting value in the maximization.
The default is {cmd:lnvar0(-1)}.

{phang}{opt usemyb} implies that a matrix called {cmd:myb0} serves as
starting values for the Prentice-Gloeckler-Meyer hazard regression with
gamma frailty.  If this matrix does not exist, then {cmd:pgmhazgf2}
proceeds with the {helpb cloglog} parameter estimates and one proposed
log-gamma variance starting value.  To see a properly constructed
starting vector, run {cmd:pgmhazgf2} and then submit the following two
commands:

{phang2}{cmd:. matrix myb0=e(b)}{p_end}
{phang2}{cmd:. matrix list myb0}{p_end}

{pmore}The order within this row vector is important for identifying the
equation and variable that each starting value is associated with.
These initial values are set with the {cmd:ml init, copy} option; see
{helpb ml init}.

{phang}{opt usemybc} implies that a matrix called {cmd:myb0c} serves as
starting values for the {helpb cloglog} regression.  If this matrix does
not exist, then the {cmd:cloglog} regression proceeds with the starting
values that the {cmd:cloglog} command uses.  To see a properly
constructed starting vector, estimate a {cmd:cloglog} model and then
submit the following two commands:

{phang2}{cmd:. matrix myb0c=e(b)}{p_end}
{phang2}{cmd:. matrix list myb0c}{p_end}

{pmore}The order within this row vector is important for identifying the
variable that each starting value is associated with.  These initial
values are set with the {cmd:ml init, copy} option; see {helpb ml init}.

{phang}{opt eform} reports the coefficients transformed to hazard ratio
format, that is, {cmdab:exp(_b[}{it:equation:varname}{cmd:])} rather
than {cmdab:_b[}{it:equation:varname}{cmd:]}.  Standard errors and
confidence intervals are similarly transformed.  {opt eform} may be
specified at estimation or when redisplaying results; see 
{manhelpi eform_option R}.  For the gamma parameter and time-varying
covariates, both the parameter estimates and the associated hazard
ratios are reported.  For time constant covariates, just the hazard
ratios are reported.

{phang}{opt noconstant} suppresses the constant term (intercept) in the model.

{phang}{opt nobeta0} suppresses reporting of the estimates from Model
(1).

{phang}{cmd:level(}{it:#}{cmd:)} species the significance level, as a
percentage, for confidence intervals. 
The default is {cmd:level(95)} or as set by {helpb set level}.

{phang}{it:maximize_options} controls the maximization process; see
{manhelp maximize R}.  For difficult maximization problems, using the
{cmd:difficult} or {cmd:technique()} option may help convergence.  So
too might different starting values.


{title:Examples}

{phang}
{cmd:. use http://www.stata-press.com/data/r10/bc.dta}{p_end}
{phang}
{cmd:. generate td = ceil(t) // construct a discrete-time variable}{p_end}
{phang}
{cmd:. isid t age dietfat}{p_end}
{phang}
{cmd:. sort t age dietfat}{p_end}
{phang}
{cmd:. generate id = _n // id variable for each person}{p_end}
{phang}
{cmd:. expand td}{p_end}
{phang}
{cmd:. sort id}{p_end}
{phang}
{cmd:. by id: generate newt = _n // new time variable}{p_end}
{phang}
{cmd:. sort id newt}{p_end}
{phang}
{cmd:. by id: generate died = dead==1 & _n==_N // identify censoring status during each time interval at risk}{p_end}
{phang}
{cmd:. generate newt_sq=newt^2 // time varying covariate}{p_end}
{phang}
{cmd:. pgmhazgf2 "died" "age smoking" "newt newt_sq", id(id) seq(newt)}{p_end}
{phang}
{cmd:. pgmhazgf2, eform}{p_end}

{pstd}
Identify how many time categories exist with at least one failure in each one.{p_end}

{phang}
{cmd:. tabulate newt died}{p_end}
{phang}
{cmd:. recode newt (8=9) (16/17=18) (20=21) (25=26) (28/29=30) (31/34=35), generate(newtd)}{p_end}
{phang}
{cmd:. tabulate newtd died}{p_end}
{phang}
{cmd:. recode age (26/35=35) (36/45=45) (46/55=55) (56/63=63) , generate(ageg)}{p_end}
{phang}
{cmd:. pgmhazgf2 "died" "i.ageg smoking" "i.newtd", id(id) seq(newt)}{p_end}


{title:Saved results}

{pstd}In addition to the usual results saved after {helpb ml},
{cmd:pgmhazgf2} also saves the following in {cmd:e()}:

{synoptset 18 tabbed}{...}
{p2col 5 18 29 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_spell)}}number of spells in the data{p_end}
{synopt:{cmd:e(gammav)}}estimated gamma variance{p_end}
{synopt:{cmd:e(se_gammav)}}standard error of the estimated gamma variance{p_end}
{synopt:{cmd:e(ll_nofr)}}log-likelihood value from model 1{p_end}
{synopt:{cmd:e(lltest)}}likelihood-ratio statistic for the test of model 1
versus model 2, that is, for the hypothesis that the gamma variance is equal
to 0{p_end}
{synopt:{cmd:e(lltest_p)}}associated p-value{p_end}

{p2col 5 25 29 2: Macros}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(idvar)}}contain the names specified in {opt id()}{p_end}
{synopt:{cmd:e(seqvar)}}contain the names specified in {opt seq()}{p_end}

{pstd}Note: After estimation, {cmd:mata: mata describe} lists the names
of the matrices and functions in memory.  The real column vector
{cmd:v1} is a Mata view of the ID variable, and the real matrix {cmd:v2}
is a Mata view of the time-varying covariates.  The real matrix ID info
is constructed with the command {cmd:idinfo=panelsetup(v1,1)}; see
{manhelp mf_panelsetup M-5:panelsetup()}.


{title:Acknowledgments and contributions of this program}

{pstd}This program is an extension of the {cmd:d0} maximum likelihood
{cmd:pgmhaz8} program written by Jenkins (2004b).  Furthermore, the
examples and several parts of this document are directly from the
{cmd:pgmhaz8} help document (Jenkins 2004b).  The following command
provides a description of this program:

{phang}
{cmd:. ssc describe pgmhaz8}{p_end}

{pstd}The main contribution is expression of the gradient vector and
Hessian matrix under a {helpb cloglog} model that accounts for gamma
distributed unobserved heterogeneity.  This increases the speed at which
a likelihood function is maximized and it allows one to use computer
programs that require an expression of the gradient vector.

{pstd}While developing this program, advice from Stephen P. Jenkins was
essential, and his comments were very helpful.  Isabel Canette, who has
a PhD in Statistics/Data Analysis and is a Senior Statistician with
Stata, was also very helpful in answering questions I had while working
on this.


{title:References}

{phang}
Allison, P. D.  1982.  Discrete-time methods for the analysis of
event histories.  {it:Sociological Methodology} 13: 61-98.

{phang}
Heckman, J. J., and B. Singer.  1984.  A method for minimizing the
impact of distributional assumption in econometric models for duration
data.  {it:Econometrica} 52: 271-320.

{phang}
Jenkins, S. P.  1995.  Easy estimation methods for discrete-time
duration models.  {it:Oxford Bulletin of Economics and Statistics} 57: 129-138.

{phang}
------.  2004a.  hshaz: Stata module to estimate discrete time (grouped data)
proportional hazards models.  Statistical Software Components S438501,
Department of Economics, Boston College.
{browse "http://ideas.repec.org/c/boc/bocode/s438501.html":http://ideas.repec.org/c/boc/bocode/s438501.html}.

{phang}
------.  2004b.  pgmhaz8: Stata module to estimate discrete time
(grouped data) proportional hazards models.  Statistical Software Components
S438501, Department of Economics, Boston College.
{browse "http://ideas.repec.org/c/boc/bocode/s438501.html":http://ideas.repec.org/c/boc/bocode/s438501.html}.

{phang}
------.  2008.  Survival Analysis with Stata.  University of
Essex module EC968.
{browse "http://www.iser.essex.ac.uk/study/resources/module-ec968":http://www.iser.essex.ac.uk/study/resources/module-ec968}.

{phang}
Meyer, B. D.  1990.  Unemployment insurance and unemployment spells.
{it:Econometrica} 58: 757-782.

{phang}
Prentice, R. L., and L. A. Gloeckler.  1978.  Regression analysis
of grouped survival data with application to breast cancer data.
{it:Biometrics} 34: 57-67.


{title:Author}

{pstd}Michael G. Farnworth{p_end}
{pstd}Department of Economics{p_end}
{pstd}University of New Brunswick{p_end}
{pstd}Fredericton, New Brunswick, Canada{p_end}
{pstd}{browse "mailto:MikeFarn@unb.ca":MikeFarn@unb.ca}{p_end}


{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 12, number 2: {browse "http://www.stata-journal.com/article.html?article=st0256":st0256}

{p 7 14 2}Note: 
A PDF help file has also been downloaded with this software.
The PDF help file provides more details about the survival and maximum
likelihood functions.  Furthermore, an appendix derives the relationship
between the gamma-distributed term and the survival function.  The appendix
also reports the likelihood function, gradient vector, and Hessian matrix.

{p 7 14 2}Help:  {helpb pgmhaz8} (if installed), {helpb stcox},
{helpb cloglog}, {helpb xtcloglog}, {helpb expand},
{helpb hshaz} (if installed){p_end}
