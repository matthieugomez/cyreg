{smcl}
{* *! version 0.1 7feb2021}{...}
{viewerjumpto "Syntax" "cyreg##syntax"}{...}
{viewerjumpto "Description" "cyreg##description"}{...}
{viewerjumpto "Options" "cyreg##options"}{...}
{viewerjumpto "Examples" "cyreg##examples"}{...}
{viewerjumpto "References" "cyreg##references"}{...}
{viewerjumpto "Author" "cyreg##contact"}{...}



{title:Title}

{p2colset 4 24 24 8}{...}
{p2col :{cmd:cyreg} {hline 2}}Implements "Efficient tests of stock return predictability" Campbell Yogo (2006){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:cyreg} {depvar} {indepvars}{cmd:,} [ {help cyreg##options:options}]{p_end}

{marker description}{...}
{title:Description}

{pstd}
The program computes the 90% confidence band using the CY method.

{marker options}{...}
{title:Options}

{synoptset 25 tabbed}{...}
{synoptline}

{synopt:{opt nog:graph}} Does not display the graph corresponding to the 90% CI {p_end}
{synopt:{opth lag(numlist)}} Impose number of lags. When left unspecified, the optimal number of lag is chosen via BIC criterion.


{marker examples}{...}
{title:Examples}
{pstd} Prepare dataset{p_end}
{phang2}{cmd:. set obs 100}{p_end}
{phang2}{cmd:. gen time = _n}{p_end}
{pstd} Generate serially correlated predictor and returns{p_end}
{phang2}{cmd:. gen predictor = 1}{p_end}
{phang2}{cmd:. replace predictor = 0.8 * predictor[_n-1] + rnormal() if _n > 1}{p_end}
{phang2}{cmd:. gen return = 0.3 * predictor[_n-1] + rnormal() if _n > 1}{p_end}
{pstd} Run regression {p_end}
{phang2}{cmd:. tsset time}{p_end}
{phang2}{cmd:. cyreg return L.predictor}{p_end}

{phang}
Please report issues on Github
{browse "https://github.com/matthieugomez/cyreg":https://github.com/matthieugomez/cyreg}
{p_end}


