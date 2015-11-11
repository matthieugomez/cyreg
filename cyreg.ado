cap program drop cyreg
program define cyreg, eclass
        syntax varlist(min=2 numeric ts) [if], tablepath(string)

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
        replace `touse' = 0 if missing(F.`x')
        qui count if `touse'
        local T = r(N)
        tempvar tempe tempesq
        qui reg F.`x' `x' if `touse'
        qui predict `tempe', res
        qui gen `tempesq' = `tempe'^2
        qui sum `tempesq' if `touse'
        local sigma_e = sqrt(r(sum) / (`T' - 2))
        tempvar tempu tempusq
        qui reg `y' `x' if `touse'
        qui predict `tempu', res
        qui gen `tempusq' = `tempu'^2
        qui sum `tem2sq' if `touse'
        local sigma_u = sqrt(r(sum) / (`T' - 2))

        tempvar tempue
        qui gen `tempue' = `tempe' * `tempu'
        qui sum `tempue' if `touse'
        local sigma_ue = r(sum) /  (`T' - 2)

        local delta = `sigma_ue' / (`sigma_u' * `sigma_e')
        aux_dfgls F.`x' if `touse'
        local tstat = e(tstat)
        preserve
        use "`tablepath'", clear
        qui gen delta_dist = abs(`delta' - delta)
        qui gen tstat_dist = abs(`tstat' - tstat)
        sort delta_dist tstat_dist
        cap assert delta_dist[1] <= 0.025
        if _rc{
                display as error "delta `delta' is outside the range -1/-0.025"
                exit
        }
        cap assert tstat_dist[1] <= 0.1
        if _rc{
                display as error "tstat `tstat' is outside the range -5/1"
                exit
        }
        local rhomin = 1 + `=cmin[1]' / `T'
        local rhomax = 1 + `=cmax[1]' / `T'
        restore

        foreach suffix in min max{
                tempvar y`suffix'
                qui gen `y`suffix'' = `y' - `sigma_ue' / `sigma_e'^2 * (`x' - `rho`suffix'' * L.`x')
                qui reg `y`suffix'' `x' if `touse'
                local b`suffix'min = _b[`x'] - 1.645 * sqrt((1- `delta'^2)) * _se[`x']
                local b`suffix'max = _b[`x'] + 1.645 * sqrt((1- `delta'^2)) * _se[`x']
        }
        ereturn clear
        ereturn scalar rhomin = `rhomin'
        ereturn scalar bminmin = `bminmin'
        ereturn scalar bminmax = `bminmax'
        ereturn scalar rhomax = `rhomax'
        ereturn scalar bmaxmin = `bmaxmin'
        ereturn scalar bmaxmax = `bmaxmax'
        ereturn local ci = "[`bminmax', `bmaxmin']"

        local yrange = `=e(bminmax)'/2
        local ylabelmin = -ceil(0.5 * e(bmaxmax) * 10) / 10
        local ylabelmax = ceil(1.5 * e(bmaxmax) * 10) / 10
        twoway (scatteri `=e(bminmax)' `rhomin' `=e(bmaxmax)' `rhomax', connect(l)) (scatteri `=e(bminmin)' `rhomin' `=e(bmaxmin)' `rhomax', connect(l)), xline(1, lpattern(solid) lcolor(black)) yline(0, lpattern(solid) lcolor(black)) ylabel(`ylabelmin'(0.1)`ylabelmax') legend(order(1 "upper" 2 "lower"))

end


cap program drop aux_dfgls
program define aux_dfgls, eclass
        syntax varlist(min=1 numeric ts) [if] [in]
        marksample touse
        fvrevar `varlist' if `touse'
        local x = r(varlist)
        tempvar newx newrho xbar dxbar
        qui tsset
        local time = r(timevar)
        qui sum `time' if `touse'
        local mintime = r(min)
        local maxtime = r(max)
        scalar rho_gls = 1 - 7 / (`maxtime' - `mintime' + 1)
        qui gen `newx' = `x' if `time' == `mintime' & `touse'
        qui replace `newx' = `x' - `=rho_gls' * L.`x' if `time' > `mintime'

        qui gen `newrho' = 1 if `time' == `mintime'  & `touse'
        qui replace `newrho' = 1 -`=rho_gls' if `time' > `mintime'  & `touse'

        qui reg `newx' `newrho' if `touse', nocons

        qui gen `xbar' = `x' - _b[`newrho'] if `touse'

        qui gen `dxbar' = `xbar' - L.`xbar' if `touse'
        qui reg `dxbar' L.`xbar' L.`dxbar' if `touse', nocons
        local tstat = `=_b[L1.`xbar']'
        ereturn clear
        ereturn scalar tstat = `tstat'
end
