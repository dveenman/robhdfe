{smcl}
{* *! version 1.0.0 2026XXXX David Veenman}{...}
{vieweralsosee "[R] regress" "help regress"}{...}
{vieweralsosee "[R] reghdfe" "help reghdfe"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[R] robreg" "help robreg"}{...}
{viewerjumpto "Syntax" "robhdfe##syntax"}{...}
{viewerjumpto "Description" "robhdfe##description"}{...}
{viewerjumpto "Options" "robhdfe##options"}{...}
{viewerjumpto "Examples" "robhdfe##examples"}{...}
{viewerjumpto "Stored results" "robhdfe##results"}{...}
{viewerjumpto "References" "robhdfe##references"}{...}
{viewerjumpto "Author" "robhdfe##author"}{...}
{title:Title}

{pstd}{hi:robhdfe} {hline 2} Robust high-dimensional fixed effects estimation with clustered standard errors

{marker syntax}{...}
{title:Syntax}

{p 4 17 2}
{cmd:robhdfe} {it:subcommand} {depvar} {indepvars} {ifin}
[{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt :{opt absorb(varlist)}}specify the fixed effects to absorb{p_end}
{synopt :{opt efficiency(real)}}normal efficiency of the estimator (must be between 63.7 and 99.9){p_end}

{syntab:Optional}
{synopt :{opt cluster(varlist)}}cluster standard errors by specified variables{p_end}
{synopt :{opt tolerance(real)}}tolerance for convergence (default 1e-10){p_end}
{synopt :{opt weightvar(varname)}}store robust weights in specified variable{p_end}
{synopt :{opt keepsin}}do not drop singleton observations{p_end}
{synoptline}
{p2colreset}{...}

{phang}
{it:subcommand} must be {cmd:m}, for Huber M-estimation. 


{marker description}{...}
{title:Description}

{phang}
    {cmd:robhdfe} performs a robust regression estimation with high-dimensional fixed effects using the Huber 
    objective function, as implemented in Gassen & Veenman (2026). It is designed to handle estimation with many fixed 
    effects and provides resistance to outliers and heavy tails in the dependent variable. The command uses an iterative 
    reweighted least squares (IRWLS) algorithm with initial estimates determined by a scale parameter obtained using 
    residuals from a quantile regression with fixed effects following methods of moments approach (MM-QR) of Machado and 
    Santos Silva (2019) and the extension to multiple fixed effects dimensions based on Rios-Avila, Siles, and 
    Canavire-Bacarreza (2024).

{phang}
    The estimator achieves the specified normal efficiency while being robust to deviations from normality and outliers 
    in the regression errors. Internally, fixed effects are absorbed using the {help reghdfe} package, allowing for efficient 
    estimation with large fixed effects dimensions. For Stata versions 19 and later, the program relies on the built-in 
    {help areg} function for the IRWLS step, which similarly allows for the absorption of multiple fixed effect dimensions.

{phang}
    Standard errors can be clustered to account for dependence within (up to two) groups. Degrees-of-freedom and finite-sample 
    corrections follow the implementation in {cmd:reghdfe}.

{phang}
    This command requires the {cmd:moremata}, {cmd:reghdfe}, and {cmd:hdfe} packages.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt absorb(varlist)} specifies the variables representing the fixed effects to be absorbed. All variables in {it:varlist} must be different and not nested within each other.

{phang}
{opt efficiency(real)} sets the normal efficiency of the M-estimator, expressed as a percentage. Must be between 63.7 and 99.9. Higher values provide more efficient estimates under normality, but less robustness.

{dlgtab:Optional}

{phang}
{opt cluster(varlist)} specifies variables to cluster the standard errors on. Up to two clustering dimensions are supported. If not specified, heteroskedasticity-robust standard errors are computed.

{phang}
{opt tolerance(real)} sets the convergence tolerance for the IRWLS algorithm. Default is 1e-10. 

{phang}
{opt weightvar(varname)} stores the final robust weights in the specified new or existing variable. If the variable already exists, it will be replaced.

{phang}
    {opt keepsin} specifies that singleton observations are not dropped in the estimation. The default is to drop singletons. Different 
    from OLS estimation, singletons not only affect standard error calculations through degrees-of-freedom adjustment, but also 
    affect the initial MM-QR scale estimate. Singleton observations do not contribute to the estimation of slope coefficients, so 
    their residuals are effectively zero and not informative about the dispersion of the true error distribution. As a result, retaining 
    these observations leads to a downward bias in the scale estimate. 


{marker examples}{...}
{title:Examples}

{pstd}Basic usage with single fixed effect:{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. reghdfe price weight length, absorb(rep78)}{p_end}
{phang2}{cmd:. robhdfe m price weight length, absorb(rep78) eff(95)}{p_end}

{pstd}With multiple fixed effects and clustering:{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. reghdfe ln_w grade age ttl_exp tenure not_smsa south, absorb(idcode year) cluster(idcode)}{p_end}
{phang2}{cmd:. robhdfe m ln_w grade age ttl_exp tenure not_smsa south, absorb(idcode year) cluster(idcode) eff(95)}{p_end}

{pstd}Storing robust regression weights:{p_end}
{phang2}{cmd:. robhdfe m ln_w grade age ttl_exp tenure not_smsa south, absorb(idcode year) cluster(idcode) eff(95) weightvar(w95)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:robhdfe} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(r2_p)}}pseudo R-squared{p_end}
{synopt:{cmd:e(scale)}}robust scale estimate{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters (if clustered){p_end}
{synopt:{cmd:e(N_clust1)}}number of clusters in first dimension (if double-clustered){p_end}
{synopt:{cmd:e(N_clust2)}}number of clusters in second dimension (if double-clustered){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}
Correia, S. (2017). Linear Models with High-Dimensional Fixed Effects: An Efficient and Feasible Estimator. {browse "http://scorreia.com/research/hdfe.pdf"}.

{phang}
Gassen, J., Veenman, D. (2026) Estimation Precision and Robust Estimation in Archival Research. {browse "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4975569"}.

{phang}
Huber, P.J. (1964). Robust Estimation of a Location Parameter. {it:Annals of Mathematical Statistics}, 35(1), 73-101.

{phang}
Machado, J.A.F., Santos Silva, J.M.C. (2019). Quantiles Via Moments. {it:Journal of Econometrics}, 213(1), 145-173.

{phang}
Rios Avila, F., Siles, L., Canavire-Bacarreza, G.J. (2024). Estimating Quantile Regressions with Multiple Fixed Effects Through Method of Moments. {browse "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4944894"}.

{marker author}{...}
{title:Author}

{pstd}
David Veenman{p_end}
{pstd}
Amsterdam Business School, University of Amsterdam{p_end}
{pstd}
d.veenman@uva.nl{p_end}


