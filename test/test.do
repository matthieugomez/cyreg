/***************************************************************************************************
test.do
Validate cyreg against published Campbell-Yogo (2006) tables.
ldp means log dividend-price ratio.
***************************************************************************************************/

cap program drop cyreg
cap program drop bic
cap program drop HeaderDisplay

adopath + "table"
qui do "cyreg.ado"
import delimited using "test/CRSP_A.txt", clear
tsset time

/* Table 4 Panel A Annual (full sample 1926-2010) */
gen ep = ret - rf
cyreg F.ret ldp, nograph
assert e(nlag) == 1
assert (e(delta) + 0.721) <= 0.01
assert (e(DFGLS) + 1.033) <= 0.1
assert abs(e(minrho) - 0.903) <= 0.01
assert abs(e(maxrho) - 1.050) <= 0.01
/* eyeballing from figure 4 */
assert abs(e(Qminrho) - .925) <= 0.01
assert abs(e(Qmaxrho) - 1.01) <= 0.01
/* Table 5 Panel A Annual (reported value by CY is scaled by bscale) */
assert abs(e(Qmaxbmin) * e(bscale) - 0.014) <= 0.01
assert abs(e(Qminbmax) * e(bscale) - 0.188) <= 0.01

/* Table 4 Panel C Annual (subsample 1952-2002) */
cyreg F.ret ldp if inrange(time, 1952, 2002), nograph
assert e(nlag) == 1
assert (e(delta) + 0.749) <= 0.01
assert (e(DFGLS) + 0.462) <= 0.1
assert abs(e(minrho) - 0.917) <= 0.01
assert abs(e(maxrho) - 1.087) <= 0.01
assert abs(e(Qmaxbmin) * e(bscale) + 0.007) <= 0.01
assert abs(e(Qminbmax) * e(bscale) - 0.183) <= 0.01

/* Sign inversion */
gen temp = - ldp
cyreg F.ret temp, nograph
cyreg F.ret temp if inrange(time, 1952, 2002), nograph

/* lag() option */
cyreg F.ret ldp, lag(1) nograph
assert e(nlag) == 1
