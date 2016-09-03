/***************************************************************************************************

***************************************************************************************************/
cd "/Users/Matthieu/Dropbox/Github/stata-cyreg"
cap program drop cyreg
cap program drop bic
cap program drop make_table1_appendix 
cap program drop make_table2_appendix 

do cyreg.ado
insheet using "test/CRSP_A.txt", clear
tsset time
/* Test with table in Campblel 
Compare graph with Fig 4 A*/
gen ep = ret - rf
cyreg F.ret ldp

assert e(nlag) ==  1
assert (e(delta) + 0.7213) <= 0.01
assert (e(DFGLS) + 1.0) <= 0.1
assert abs(e(minrho) - .903) <= 0.01
assert abs(e(maxrho) - 1.050) <= 0.01

assert abs(e(minrho) - .903) <= 0.01
assert abs(e(maxrho) - 1.050) <= 0.01


cyreg F.ret lep

