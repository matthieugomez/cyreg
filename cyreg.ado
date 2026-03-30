* Compute Campbell and Yogo(2006) confidence band
* the syntax is such that y = r_t and x = x_{t-1}
* Should be testing for a positive predictive coefficient


* users should We set the maximum lag length p to four for annual, six for quarterly, and eight for monthly data

program define cyreg, eclass
  version 12.1
  syntax varlist(min=2 numeric ts) [if] [in],[ Lag(integer 0) Nlag(integer -1) MAXlag(integer 4) NOGraph * ]
  if `nlag' >= 0 {
    display as text "Warning: nlag() is deprecated, use lag() instead"
    if `lag' > 0 {
      display as error "Cannot specify both lag() and nlag()"
      exit 198
    }
    local lag = `nlag'
  }
  ereturn clear
  cap tsset
  if _rc{
    display as error "dataset must be tsset"
    exit 198
  }

  marksample touse
  tokenize `varlist'
  local originaly `1'
  local originalx `2'
  fvrevar `varlist' if `touse'
  local varlist = r(varlist)
  tokenize `varlist'
  local y `1'
  local x `2'
  tempvar tousex
  qui gen `tousex' = !missing(`x') `if' `in'
  qui replace `tousex' = 0 if `tousex' == .
  if `lag' == 0 {
    bic F.`x', maxlag(`maxlag')
    local biclag = r(optlag)
  }
  else{
    local biclag = `lag'
  }
  if `biclag' > 1{
        * I *think* T equals number of periods used in the estimation
        markout `touse' L(1/`=`biclag'-1').`x'
      } 
      qui count if `touse'
      local T = r(N)


      /* compute b, V and check delta < 0 */
      tempvar resu
      qui reg `y' `originalx' if `touse'
      tempname b V
      mat `b' = e(b)
      mat `V' = e(V)
      qui predict `resu', res

      tempvar rese dx
      qui gen `dx' = F.`x' - `x'
      if `biclag' == 1{
        qui reg `dx' `x' if `touse'
      }
      else{
        qui reg `dx' `x' L(1/`=`biclag'-1').`dx' if `touse'
      }
      qui predict `rese', res
      tempvar resue
      qui gen `resue' = `rese' * `resu'
      qui sum `resue' if `touse'
      local invx = r(mean) > 0
      if `invx'{
        local inv
        tempvar newx
        qui gen `newx' = -`x'
      }
      else{
        local newx `x'
      }


    * ----------------------------
    * Step 1: run OLS regressions
    * ----------------------------
    * run the predictive regression
    tempvar resu resusq
    qui reg `y' `newx' if `touse', `options'
    local SEbeta = _se[`newx']
    local r2 = e(r2)
    local r2_a = e(r2_a)

    qui predict `resu', res
    qui gen `resusq' = `resu'^2
    qui sum `resusq' if `touse'
    local sigma_u = sqrt(r(sum) / (`T' - 2))

    * run the AR(p) estimation
    tempvar rese resesq
    tempvar dx
    qui gen `dx' = F.`newx' - `newx'
    local sumpsi = 0
    if `biclag' == 1{
      qui reg `dx' `newx' if `touse'
    }
    else{
      qui reg `dx' `newx' L(1/`=`biclag'-1').`dx' if `touse'
      foreach i of numlist 1/`=`biclag'-1'{
        local sumpsi = `sumpsi' + _b[L`i'.`dx']
      }
    }
    qui predict `rese', res
    qui gen `resesq' = `rese'^2
    qui sum `resesq' if `touse'
    local sigma_e = sqrt(r(sum) / (`T' - 2))
    local omega = sqrt(`sigma_e'^2/(1-`sumpsi')^2)
    * compute covariance of errors
    tempvar resue
    qui gen `resue' = `rese' * `resu'
    qui sum `resue' if `touse'
    local sigma_ue = r(sum) /  (`T' - 2)

    local delta = `sigma_ue' / (`sigma_u' * `sigma_e')

    * ----------------------------
    * Step 2: Run the AR(1) regression
    * ----------------------------
    tempvar tempv tempvsq
    qui reg F.`newx' `newx' if `touse'
    qui predict `tempv', res
    qui gen `tempvsq' = `tempv'^2
    qui sum `tempvsq' if `touse'
    local sigma_v = sqrt(r(sum) / (`T' - 2))
    local SErho = _se[`newx']

    * ----------------------------
    * Step 3: Compute the DF-GLS statistic for rho
    * ----------------------------
    /* use tousex to mimic cambell yogo table*/
    qui reg F.`newx' `newx' if `tousex'
    local rho = _b[`newx']
    qui dfgls `newx' if `tousex', maxlag(`=`biclag'-1') notrend
    if `biclag' == 1{
      local tstat = r(dft0)
    }
    else if `biclag' > 1{
      matrix r = r(results)
      local tstat = r[`=`biclag'-1', 5]
    }
    * check whether values are in table
    if  !inrange(`tstat', -5, 1) {
      display as error "DF-GLS tstat `tstat' is outside the range -5/1"
      exit
    }


    * ----------------------------
    * Step 4: Compute the confidence interval for c
    * ----------------------------
    * compute ci for rho
    preserve
    qui findfile cyreg_table1_appendix.dta
    qui use "`r(fn)'", clear
    qui keep if size == "95%"
    qui gen tstat_dist = abs(`tstat' - tstat)
    sort tstat_dist
    local cmin = cmin[1]
    local cmax = cmax[1]
    local minrho = 1 + `cmin' / `T'
    local maxrho = 1 + `cmax' / `T'
    restore

    * pretest. Look if confidence interval reject unit root clearly enough to have valid statistic
    preserve
    qui findfile cyreg_table1.dta
    qui use "`r(fn)'", clear
    rename (v1 v2 v3) (delta cmin cmax)
    qui gen delta_dist = abs(`delta' - delta)
    sort delta_dist
    if `cmin' >= `=cmax[1]' | `cmax' <= `=cmin[1]'{
      local pretest "No"
    }
    else{
      local pretest "Yes"
    }
    restore

    /* now we're picking the  CI in the Bonferroni q test */
    preserve
    qui findfile cyreg_table2_appendix.dta
    qui use "`r(fn)'", clear
    qui gen tstat_dist = abs(`tstat' - tstat)
    qui gen delta_dist = abs(`delta' - delta)
    sort delta_dist tstat_dist
    local Qminrho = 1 + `=cmin[1]' / `T'
    local Qmaxrho = 1 + `=cmax[1]' / `T'
    local Qminc = cmin[1]
    local Qmaxc = cmax[1]
    restore

    * ----------------------------
    * Step 4: Compute the Bonferroni interval
    * note that omega = sigma_e if p == 1 so the second term disappears in the summation
    * ----------------------------
    if `biclag' == 1{
      assert abs(`sigma_e' - `sigma_v') <= 1e-5
    }
    foreach suffix in min max{
      tempvar y`suffix'
      qui gen `y`suffix'' = `y' - (`sigma_ue') / (`sigma_e' * `omega') * ///
      (F.`newx' - (`Q`suffix'rho') * `newx')
      qui reg `y`suffix'' `newx' if `touse'
      local Q`suffix'bmin = _b[`newx'] ///
      + (`T'-2)/2 * (`sigma_ue') / ((`sigma_e') * (`omega')) * ((`omega')^2/(`sigma_v')^2 - 1) * (`SErho')^2 ///
      - 1.645 *  sqrt((1 - (`delta')^2))  * (`SEbeta')
      local Q`suffix'bmax = _b[`newx'] ///
      + (`T'-2)/2 * (`sigma_ue') / ((`sigma_e') * (`omega')) * ((`omega')^2/(`sigma_v')^2 - 1) * (`SErho')^2 ///
      + 1.645 * sqrt((1 - (`delta')^2))  * (`SEbeta')
    }
    * output values




    ereturn clear
    ereturn post `b' `V', obs(`T') esample(`touse') dof(2) depname("`originaly'") 
    ereturn scalar r2 = `r2'
    ereturn scalar r2_a = `r2_a'
    /* non corrected regression */
    /* correlation */
    if `invx'{
      ereturn scalar delta = -`delta'
    }
    else{
      ereturn scalar delta = `delta'
    }

    /* ar(p) structure */
    ereturn local pretest = "`pretest'"
    ereturn scalar N = `T'
    ereturn scalar nlag = `biclag'
    ereturn scalar DFGLS = `tstat'
    ereturn scalar minrho = `minrho'
    ereturn scalar rho = `rho'
    ereturn scalar maxrho = `maxrho'
    ereturn scalar Qminc = `Qminc'
    ereturn scalar Qmaxc = `Qmaxc'
    ereturn scalar Qminrho = `Qminrho'
    ereturn scalar Qmaxrho = `Qmaxrho'

    /* beta */
    if `invx'{
      ereturn scalar Qminbmax = -`Qminbmin'
      ereturn scalar Qminbmin = -`Qminbmax'
      ereturn scalar Qmaxbmax = -`Qmaxbmin'
      ereturn scalar Qmaxbmin = -`Qmaxbmax'
    }
    else{
      ereturn scalar Qminbmin = `Qminbmin'
      ereturn scalar Qminbmax = `Qminbmax'
      ereturn scalar Qmaxbmin = `Qmaxbmin'
      ereturn scalar Qmaxbmax = `Qmaxbmax'
    }
    ereturn scalar bscale = (`sigma_e')/(`sigma_u')

    * ----------------------------
    * Step 5-6: Produce the chart if needed
    * ----------------------------
    ereturn display


    tempname left right
    .`left' = {}
    .`right' = {}
    local width 78
    local colwidths 1 30 41 67
    local i 0
    foreach c of local colwidths {
      local ++i
      local c`i' `c'
      local C`i' _col(`c')
    }
    local max_len_title = `c3' - 2
    local c4wfmt1 = `c4wfmt' + 1
    .`right'.Arrpush `C3' "Number of obs" `C4' "= " as res %10.3g e(N)
    .`right'.Arrpush `C3' "Number of lags" `C4' "= " as res %10.3g e(nlag)
    .`right'.Arrpush `C3' "Correlation (delta)" `C4' "= " as res %10.3f e(delta)
    .`right'.Arrpush `C3' "Persistence (rho) 95 % CI" `C4' "= " "[`: display %4.3f e(minrho)', `: display %4.3f e(maxrho)']"
    .`right'.Arrpush `C3' "Pretest : does the usual t-test exceeds 7.5%?" `C4' "= " "`pretest'"

    if `invx'{
      .`right'.Arrpush `C3' "Coefficient (beta) 90% CI" `C4' "= "  "[`: display %4.3f e(Qminbmin)', `: display %4.3f e(Qmaxbmax)']"
    }
    else{
      .`right'.Arrpush `C3' "Coefficient (beta) 90% CI" `C4' "= "  "[`: display %4.3f e(Qmaxbmin)', `: display %4.3f e(Qminbmax)']"
    }
    HeaderDisplay `left' `right' `"`title'"' `"`title2'"' `"`title3'"' `"`title4'"' `"`title5'"'
    if "`nograph'" == "" {
      local ylabelmin = floor(10 * min(e(Qminbmin), e(Qmaxbmin), -0.1)) / 10
      local ylabelmax = ceil(10 * max(e(Qminbmax), e(Qmaxbmax), 0.1)) / 10
      local xlabelmin = ceil((`Qminrho' - 0.5 * (`Qmaxrho' - `Qminrho')) * 100) / 100
      local xlabelmax = ceil((max(1,`Qmaxrho') + 0.5 * (`Qmaxrho' -`Qminrho')) * 100) / 100
      twoway (scatteri `=e(Qminbmax)' `Qminrho' `=e(Qmaxbmax)' `Qmaxrho' ///
      , connect(l)) ///
      (scatteri `=e(Qminbmin)' `Qminrho' `=e(Qmaxbmin)' `Qmaxrho' ///
      , connect(l)) ///
      , xline(1, lpattern(solid) lcolor(black)) ///
      yline(0, lpattern(solid) lcolor(black)) ///
      ylabel(`ylabelmin'(0.1)`ylabelmax', angle(0)) ///
      xlabel(`xlabelmin'(0.02)`xlabelmax') ///
      legend(order(1 "upper" 2 "lower")) ///
      xtitle("Persistence {&rho}") ytitle("Predictive coefficient {&beta}") ///
      graphregion(fcolor(white))  bgcolor(white)
    }
  end

  program define bic, rclass
    syntax varlist(min=1 numeric ts) [if] [in], [maxlag(integer 10)]
    marksample touse
    fvrevar `varlist' if `touse'
    local x = r(varlist)
    cap assert `maxlag' >= 1
    if _rc{
      di as error "maxlag must be >= 1"
      exit
    }
    qui varsoc `x' if `touse', maxlag(`maxlag')
    mat s = r(stats)
    local optlag = 0
    local bic = .
    foreach nlag of numlist 1/`maxlag'{
      if s[1 + `nlag',  9] < `bic'{
        local optlag = `nlag'
        local bic = s[1 + `nlag',  9]
      }
    }
    return scalar optlag = `optlag'
  end


  program define HeaderDisplay
    args left right title1 title2 title3 title4 title5

    local nl = `.`left'.arrnels'
    local nr = `.`right'.arrnels'
    local K = max(`nl',`nr')

    di
    if `"`title1'"' != "" {
      di as txt `"`title'"'
      forval i = 2/5 {
        if `"`title`i''"' != "" {
          di as txt `"`title`i''"'
        }
      }
      if `K' {
        di
      }
    }

    local c _c
    forval i = 1/`K' {
      di as txt `.`left'[`i']' as txt `.`right'[`i']'
    }
  end
