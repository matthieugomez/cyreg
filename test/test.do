/***************************************************************************************************
ldp lep means log.
***************************************************************************************************/
cd "/Users/Matthieu/Dropbox/Github/stata-cyreg"
cap program drop cyreg
cap program drop bic
cap program drop make_table1_appendix 
cap program drop make_table2_appendix 

qui do cyreg.ado
insheet using "test/CRSP_A.txt", clear
tsset time

/* mimic values from Table 4 Panel A Annual*/
* first test ldp not using L. 
gen ep = ret - rf
cyreg F.ret ldp
assert e(nlag) ==  1
assert (e(delta) + 0.721) <= 0.01 
assert (e(DFGLS) + 1.033) <= 0.1 
assert abs(e(minrho) - 0.903) <= 0.01 
assert abs(e(maxrho) - 1.050) <= 0.01 
/* eyeballing from figure 4*/
assert abs(e(Qminrho) - .925) <= 0.01 
assert abs(e(Qmaxrho) - 1.01) <= 0.01 
/*Table 5 Panel A Annual */
assert abs(e(Qmaxbmin) * e(bscale) - 0.014) <= 0.01 
assert abs(e(Qminbmax) * e(bscale) - 0.188) <= 0.01 

/* mimic values from Table 4 Panel C Annual */
cyreg F.ret ldp if inrange(time, 1952, 2002)
assert e(nlag) ==  1
assert (e(delta) + 0.749) <= 0.01 
assert (e(DFGLS) + 0.462) <= 0.1 
assert abs(e(minrho) - 0.917) <= 0.01 
assert abs(e(maxrho) - 1.087) <= 0.01 
assert abs(e(Qmaxbmin) * e(bscale) + 0.007) <= 0.01 
assert abs(e(Qminbmax) * e(bscale) - 0.183) <= 0.01 
/* now the weird thing ius that rho overline is 1.01 in their Fig4 but checking the value for delta and DFGLS in Table4 I obtain much higher Qmaxrho */


/* check that i obtain other fic with this */
gen temp = - ldp
cyreg F.ret temp
cyreg F.ret temp if inrange(time, 1952, 2002)








