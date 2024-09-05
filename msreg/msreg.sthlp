{smcl}
{* *! version 1.0.0  23mar2020}{...}
{viewerjumpto "Title" "msreg##title"}{...}
{viewerjumpto "Syntax" "msreg##syntax"}{...}
{viewerjumpto "Description" "msreg##description"}{...}
{viewerjumpto "Options" "msreg##options"}{...}
{viewerjumpto "Examples" "msreg##examples"}{...}
{viewerjumpto "Stored results" "msreg##results"}{...}
{viewerjumpto "Reference" "msreg##reference"}{...}
{viewerjumpto "Authors" "msreg##authors"}{...}
{viewerjumpto "Also see" "msreg##alsosee"}{...}
{cmd:help msreg}{right: ({browse "https://doi.org/10.1177/1536867X211000008":SJ21-1: st0630})}
{hline}

{marker title}{...}
{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{cmd:msreg} {hline 2}}Consistent estimation of linear regression models
using matched data{p_end}


{marker syntax}{...}
{title:Syntax}

{p 8 14 2}
{cmd:msreg}
{depvar}
[{it:varlist_X1}]
{cmd:(}{it:varlist_X2} {cmd:=} {it:varlist_Z}{cmd:)}
{cmd:using} {help filename:{it:filename}} 
{ifin}
[{cmd:,} {it:options}]

{synoptset 19}{...}
{synopthdr}
{synoptline}
{synopt:{cmd:vce(}{help msreg##vce_spec:{it:vce_spec}}{cmd:)}}specify the type of variance and covariance
matrix used in computation{p_end}
{synopt:{cmd:estimator(}{help msreg##est_spec:{it:est_spec}}{cmd:)}}specify the type of estimator{p_end}
{synopt:{opt nneighbor(#)}}specify the number of matches per observation{p_end}
{synopt:{cmd:metric(}{help msreg##metric_spec:{it:metric_spec}}{cmd:)}}specify the distance matrix used as
the weight matrix in a quadratic form that transforms the multiple distances
into a single distance measure{p_end}
{synopt:{opt order(#)}}specify the order of polynomial in the power-series
approximation{p_end}
{synopt:{opt nocon:stant}}suppress the constant term{p_end}
{synopt:{opt l:evel(#)}}specify the level of significance for output
table{p_end}
{synopt :{help msreg##display_options:{it:display_options}}}control
INCLUDE help shortdes-displayoptall
INCLUDE help shortdes-coeflegend
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{it:varlist_X1}, {it:varlist_X2}, and {it:varlist_Z} may contain factor
variables; see {help fvvarlist}.{p_end}

{synoptset 19}{...}
{marker vce_spec}{...}
{synopthdr:vce_spec}
{synoptline}
{synopt :{opt vi}}VI-type variance-covariance assuming sample-size ratio
converges to a nonzero constant; the default{p_end}
{synopt :{opt vii}}VII-type variance-covariance assuming sample-size ratio
converges to zero{p_end}
{synopt :{opt viii}}VIII-type variance-covariance assuming sample-size ratio
converges to infinity{p_end}
{synoptline}

{marker est_spec}{...}
{synopthdr:est_spec}
{synoptline}
{synopt :{opt two:step}}specify to use the two-step bias-corrected estimator;
the default{p_end}
{synopt :{opt one:step}}specify to use the one-step bias-corrected
estimator{p_end}
{synoptline}

{marker metric_spec}{...}
{synopthdr:metric_spec}
{synoptline}
{synopt :{opt maha:lanobis}}specify to use the inverse of the sample covariance
matrix of matching variables; the default{p_end}
{synopt :{opt eucl:idean}}specify to use the inverse of only diagonal elements
of the sample covariance matrix of matching variables{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:msreg} fits a linear regression model when some of the covariates are
missing in one sample and we need to use a matched sample for estimation.  The
ordinary least-squares (OLS) estimator using the matched sample is
inconsistent in this setting.  {cmd:msreg} can consistently estimate the model
coefficients using matched sample indirect inference estimators developed in
Hirukawa and Prokhorov (2018).


{marker options}{...}
{title:Options}

{phang}
{opt vce(vce_spec)} specifies the type of variance-covariance matrix used in
computation.  {it:vce_spec} can be one of {cmd:vi}, {cmd:vii}, or {cmd:viii}.
The default is {cmd:vce(vi)}.

{phang}
{opt estimator(est_spec)} specifies the type of estimator.  {it:est_spec} can
be either {cmd:onestep} or {cmd:twostep}.  {cmd:onestep} specifies to use the
one-step bias-corrected estimator.  {cmd:twostep} specifies to use the
two-step bias-corrected estimator.  The default is {cmd:estimator(twostep)}.

{phang}
{opt nneighbor(#)} specifies the number of matches per observation.  The
default is {cmd:nneighbor(1)}.  The maximum allowed number of matches is 10.
Each observation is matched with the mean of the specified number of
observations from the other dataset.

{phang}
{opt metric(metric_spec)} specifies the distance matrix used as the weight
matrix in a quadratic form that transforms the multiple distances into a
single distance measure.  {it:metric_spec} can be either {cmd:mahalanobis} or
{cmd:euclidean}.  {cmd:metric(mahalanobis)} specifies to use the inverse of
the sample covariance matrix of matching variables, which is the default.
{cmd:metric(euclidean)} specifies to use the inverse of only diagonal elements
of the sample covariance matrix of matching variables.

{phang}
{opt order(#)} specifies the order of polynomials in the power-series
approximation for fully modified matched-sample indirect inference.  The
default is {cmd:order(2)}.  The maximum allowed number of order is 5.

{phang}
{cmd:noconstant} suppresses the constant term.

{phang}
{opt level(#)} specifies the level of significance for the output table.

{marker display_options}{...}
INCLUDE help displayopts_list

{phang}
{cmd:coeflegend} specifies that the legend of the coefficients and how to
specify them in an expression be displayed rather than displaying the
statistics for the coefficients.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. import excel "https://raw.githubusercontent.com/timbulwidodostp/msreg/main/msreg/msreg.xlsx", sheet("Sheet2") firstrow clear}{p_end}
{phang2}{cmd:. save "C:\msreg_2.dta"}{p_end}
{phang2}{cmd:. import excel "https://raw.githubusercontent.com/timbulwidodostp/msreg/main/msreg/msreg.xlsx", sheet("Sheet1") firstrow clear}{p_end}

{pstd}
Use two-step bias-corrected estimator{p_end}
{phang2}{cmd:. msreg y x11 x12 (x21 x22 = z1 z2) using "C:\msreg_2.dta", vce(vi)}{p_end}

{pstd}
Specify 3 matches per observation{p_end}
{phang2}{cmd:. msreg y x11 x12 (x21 x22 = z1 z2) using "C:\msreg_2.dta", vce(vi) nneighbor(3)}{p_end}

{pstd}
Use the third order of polynomial in the power-series approximation{p_end}
{phang2}{cmd:. msreg y x11 x12 (x21 x22 = z1 z2) using "C:\msreg_2.dta", vce(vi) order(3)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:msreg} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt: {cmd:e(N_1)}}number of observations in the first sample{p_end}
{synopt: {cmd:e(N_2)}}number of observations in the second sample{p_end}
{synopt: {cmd:e(nneighbor)}}number of nearest neighbors matched{p_end}
{synopt: {cmd:e(chi2)}}chi-squared{p_end}
{synopt: {cmd:e(p)}}p-value for test of variables{p_end}
{synopt: {cmd:e(df_m)}}degree of freedom in the model{p_end}
{synopt: {cmd:e(order)}}order of polynomial in the power-series estimation{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt: {cmd:e(cmd)}}{cmd:msreg}{p_end}
{synopt: {cmd:e(title)}}title in estimation output{p_end}
{synopt: {cmd:e(chi2type)}}{cmd:Wald}{p_end}
{synopt: {cmd:e(vce)}}{it:vce_spec} specified in {cmd:vce()}{p_end}
{synopt: {cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt: {cmd:e(properties)}}{cmd:b V}{p_end}
{synopt: {cmd:e(footnote)}}footnote displayed under the output table{p_end}
{synopt: {cmd:e(common)}}common variables that exist in both samples{p_end}
{synopt: {cmd:e(matched}}matched variables{p_end}
{synopt: {cmd:e(metric)}}type of distance matrix{p_end}
{synopt: {cmd:e(estimator)}}type of estimator{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt: {cmd:e(b)}}coefficient vector{p_end}
{synopt: {cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{marker reference}{...}
{title:Reference}

{marker HP2018}{...}
{phang}
Hirukawa, M., and A. Prokhorov. 2018. Consistent estimation of linear
regression models using matched data. {it:Journal of Econometrics} 203:
344-358.  {browse "https://doi.org/10.1016/j.jeconom.2017.07.006"}.


{marker authors}{...}
{title:Authors}

{pstd}
Masayuki Hirukawa{break}
Ryukoku University{break}
Kyoto, Japan{break}
hirukawa@econ.ryukoku.ac.jp

{pstd}
Di Liu{break}
StataCorp{break}
College Station, TX{break}
dliu@stata.com

{pstd}
Artem Prokhorov{break}
University of Sydney Business School{break}
Sydney, Australia, and{break}
CEBA, St. Petersburg State University{break}
St. Petersburg, Russia, and{break}
CIREQ, University of Montreal{break}
Montreal, Canada{break}
artem.prokhorov@sydney.edu.au


{marker alsosee}{...}
{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 21, number 1: {browse "https://doi.org/10.1177/1536867X211000008":st0630}{p_end}
