* Compute Campbell and Yogo(2006) confidence band
* the syntax is such that y = r_t and x = x_{t-1}
* Should be testing for a positive predictive coefficient


cap program drop cyreg
program define cyreg, eclass
        syntax varlist(min=2 numeric ts) [if] [in], tabledir(string) [nlag(string) NOGraph]

        cap tsset
        if _rc{
                display as error "dataset must be tsset"
        }
        marksample touse
        fvrevar `varlist' if `touse'
        local varlist = r(varlist)
        tokenize `varlist'
        local y `1'
        local x `2'
        qui count if `touse'
        * T equals number of periods
        * should it equal the number of observations in the regression (ie with L(1/p).x)?
        local T = r(N)

        if "`nlag'" == ""{
                bic F.`x', maxlag(10)
                local nlag = r(optlag)
        }

        * ----------------------------
        * Step 1: run OLS regressions
        * ----------------------------

        * run the predictive regression
        tempvar tempu tempusq
        qui reg `y' `x' if `touse'
        local bet = _b[`x']
        local SEbet = _se[`x']
        qui predict `tempu', res
        qui gen `tempusq' = `tempu'^2
        qui sum `tempusq'
        local sigma_u = sqrt(r(sum) / (`T' - 2))

        * run the AR(p) estimation
        tempvar tempe tempesq
        tempvar dx
        qui gen `dx' = F.`x' - `x'
        local sumpsi = 0
        if `nlag' == 1{
                qui reg `dx' `x' if `touse'
        }
        else{
                qui reg `dx' `x' L(1/`=`nlag'-1').`dx' if `touse'
                foreach i of numlist 1/`=`nlag'-1'{
                        local sumpsi = `sumpsi' + _b[L`i'.`dx']
                }
        }
        qui predict `tempe', res
        qui gen `tempesq' = `tempe'^2
        qui sum `tempesq'
        local sigma_e = sqrt(r(sum) / (`T' - 2))
        local omega = sqrt(`sigma_e'^2/(1-`sumpsi')^2)

        * compute covariance of errors
        tempvar tempue
        qui gen `tempue' = `tempe' * `tempu'
        qui sum `tempue' if `touse'
        local sigma_ue = r(sum) /  (`T' - 2)

        local delta = `sigma_ue' / (`sigma_u' * `sigma_e')

        * ----------------------------
        * Step 2: Run the AR(1) regression
        * ----------------------------

        tempvar tempv tempvsq
        qui reg F.`x' `x' if `touse'
        qui predict `tempv', res
        qui gen `tempvsq' = `tempv'^2
        qui sum `tempvsq'
        local sigma_v = sqrt(r(sum) / (`T' - 2))
        local SErho = _se[`x']

        * ----------------------------
        * Step 3: Compute the DF-GLS statistic for rho
        * ----------------------------
        qui dfgls F.`x' if `touse', maxlag(`=`nlag'-1') notrend
        if `nlag' == 1{
                local tstat = r(dft0)
        }
        else if `nlag' > 1{
                matrix r = r(results)
                local tstat = r[`=`nlag'-1', 5]
        }

        * ----------------------------
        * Step 3: Compute the confidence interval for c
        * ----------------------------

        * check whether values are in table
        cap assert inrange(`delta', -1 , 0)
        if _rc{
                display as error "delta `delta' is outside the range -1/-0. Use the opposite of your predictor."
                exit
        }
        cap assert inrange(`tstat', -5, 1) 
        if _rc{
                display as error "tstat `tstat' is outside the range -5/1"
                exit
        }

        * if values are in table, pick up range
        preserve
        use "`tabledir'/table2_appendix", clear
        qui gen tstat_dist = abs(`tstat' - tstat)
        qui gen delta_dist = abs(`delta' - delta)
        sort delta_dist tstat_dist
        local rhomin = 1 + `=cmin[1]' / `T'
        local rhomax = 1 + `=cmax[1]' / `T'
        local cmin = cmin[1]
        local cmax = cmax[1]
        restore

        * ----------------------------
        * Step 4: Compute the Bonferroni interval
        * note that omega = sigma_v if p == 1 so the second term disappears in the summation
        * ----------------------------


        foreach suffix in min max{
                tempvar y`suffix'
                qui gen `y`suffix'' = `y' - (`sigma_ue') / (`sigma_e'^2) * ///
                (F.`x' - (`rho`suffix'') * `x')
                qui reg `y`suffix'' `x' if `touse'
                local b`suffix'min = _b[`x'] ///
                + (`T'-2)/2 * `sigma_ue' / (`sigma_e' * `omega') * (`omega'^2/`sigma_v'^2 - 1) * `SErho' ///
                - 1.645 *  sqrt((1 - (`delta')^2))  * `SEbet'
                local b`suffix'max = _b[`x'] ///
                + (`T'-2)/2 * `sigma_ue' / (`sigma_e' * `omega') * (`omega'^2/`sigma_v'^2 - 1) * `SErho' ///
                + 1.645 * sqrt((1 - (`delta')^2))  * `SEbet'
        }
        * output values
        ereturn clear
        ereturn scalar nlag = `nlag'
        ereturn scalar bet = `bet'
        ereturn scalar delta = `delta'
        ereturn scalar DFGLS = `tstat'
        ereturn scalar cmin = `cmin'
        ereturn scalar cmax = `cmax'
        ereturn scalar rhomin = `rhomin'
        ereturn scalar rhomax = `rhomax'
        ereturn scalar bminmin = `bminmin'
        ereturn scalar bminmax = `bminmax'
        ereturn scalar bmaxmin = `bmaxmin'
        ereturn scalar bmaxmax = `bmaxmax'

        * ----------------------------
        * Step 5-6: Produce the chart if needed
        * ----------------------------
        di "Number of lags : `nlag'"
        di "Coefficient:" %9.3g `bet',  
        di "Confidence interval: ["  %9.3g `bmaxmin' "," %9.3g `bminmax' "]"

        if "`nograph'" == "" {
          local yrange = `=e(bminmax)'/2
          local ylabelmin = -ceil(0.5 * e(bmaxmax) * 10) / 10
          local ylabelmax = ceil(1.5 * e(bmaxmax) * 10) / 10
          local xlabelmin = ceil((`rhomin' - 0.5 * (`rhomax' - `rhomin')) * 100) / 100
          local xlabelmax = ceil((max(1,`rhomax') + 0.5 * (`rhomax' -`rhomin')) * 100) / 100
          twoway (scatteri `=e(bminmax)' `rhomin' `=e(bmaxmax)' `rhomax' ///
            , connect(l)) ///
          (scatteri `=e(bminmin)' `rhomin' `=e(bmaxmin)' `rhomax' ///
            , connect(l)) ///
          , xline(1, lpattern(solid) lcolor(black)) ///
          yline(0, lpattern(solid) lcolor(black)) ///
          ylabel(`ylabelmin'(0.1)`ylabelmax') ///
          xlabel(`xlabelmin'(0.02)`xlabelmax') ///
          legend(order(1 "upper" 2 "lower")) ///
          xtitle("Persistence rho") ytitle("Predictive coefficient beta") ///
          graphregion(fcolor(white)) 
  }
end



discard
cap program drop bic
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
        foreach nlag of numlist 1/`maxlag'{
                qui reg `x' L(1/`nlag').`x' if `touse'
                qui estat ic
                mat s = r(S)
                if `nlag' == 1{
                        local optlag = 1
                        local bic = s[1, 6]
                }
                else{
                        local newbic = s[1, 6]
                        if `newbic' < `bic'{
                                local optlag = `nlag'
                                local bic = `newbic'
                        }
                }
        }
        return scalar optlag = `optlag'
end

/***************************************************************************************************
* side function computing the DF-GLS statistic
* computes the same value as " dfgls `var' , maxlag(p-1) notrend " (I checked)

cap program drop aux_dfgls
program define aux_dfgls, eclass
        syntax varlist(min=1 numeric ts) [if] [in], [lag(integer 1)]
        marksample touse
        fvrevar `varlist' if `touse'
        local x = r(varlist)
        tempvar newx newrho xbar dxbar
        qui tsset
        local time = r(timevar)
        qui sum `time' if `touse'
        local T = r(max) - r(min) + 1
        local mintime = r(min)
        scalar rho_gls = 1 - 7 / `T'
        qui gen `newx' = `x' if `time' == `mintime' & `touse'
        qui replace `newx' = `x' - `=rho_gls' * L.`x' if `time' > `mintime'

        qui gen `newrho' = 1 if `time' == `mintime'  & `touse'
        qui replace `newrho' = 1 -`=rho_gls' if `time' > `mintime'  & `touse'

        qui reg `newx' `newrho' if `touse', nocons

        qui gen `xbar' = `x' - _b[`newrho'] if `touse'

        qui gen `dxbar' = `xbar' - L.`xbar' if `touse'
        if `lag' == 1{
                qui reg `dxbar' L.`xbar' if `touse', nocons
        }
        else{
                qui reg `dxbar' L.`xbar' L(1/`=`lag'-1').`dxbar' if `touse', nocons
        }
        local tstat = `=_b[L.`xbar']' / `=_se[L.`xbar']'
        ereturn clear
        ereturn scalar tstat = `tstat'
end

***************************************************************************************************/

