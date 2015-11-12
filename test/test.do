adopath + "/Users/Matthieu/Dropbox/Github/stata-cyreg"
cd "/Users/Matthieu/Dropbox/Github/stata-cyreg"
do cyreg.ado
insheet using "test/CRSP_A.txt", clear
tsset time
gen ep = ret - rf
cyreg ep L.ldp, tablepath("table.dta")
cyreg ep L.ldp, tablepath("table.dta") p(3)
cyreg ret L.ldp, tablepath("table.dta")

*??? does not give Fig4