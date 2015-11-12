adopath + "/Users/Matthieu/Dropbox/Github/stata-cyreg"
cd "/Users/Matthieu/Dropbox/Github/stata-cyreg"
do cyreg.ado
insheet using "test/CRSP_A.txt", clear
tsset time
gen ep = ret - rf
reg ep L.ldp
cyreg ep L.ldp, tablepath("table.dta")
assert e(nlag) ==  1
assert e(bet) ==  .1756569861580354
assert e(delta) ==  -.7001412489946537
assert e(DFGLS) ==  -.8822743570558439
assert e(cmin) ==  -4.920000076293945
assert e(cmax) ==  .9110000133514404
assert e(rhomin) ==  .9352631568908691
assert e(rhomax) ==  1.01198684228094
assert e(bminmin) ==  .0989143626194303
assert e(bminmax) ==  .2477809241752174
assert e(bmaxmin) ==  .0296727139488352
assert e(bmaxmax) ==  .1785392755046223


cyreg ep L.lep, tablepath("table.dta")
assert e(nlag) ==  1
assert e(bet) ==  .1746210803117575
assert e(delta) ==  -.9308960406398238
assert e(DFGLS) ==  -2.066228580582374
assert e(cmin) ==  -15.73299980163574
assert e(cmax) ==  -2.834000110626221
assert e(rhomin) ==  .7929868447153192
assert e(rhomax) ==  .9627105248601813
assert e(bminmin) ==  .1992346435669487
assert e(bminmax) ==  .2710889497909728
assert e(bmaxmin) ==  .0435807108426722
assert e(bmaxmax) ==  .1154350170666962

*tests are here only to check that changes in code does not introduce unexpected changes in results
*values have not been checked agaoinst campbell-yogo (except the nlag)