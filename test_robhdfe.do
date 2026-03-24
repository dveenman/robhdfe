clear all
set obs 1000
set seed 1234567
local years 10

gen n=_n
gen firm=ceil(_n/`years') 
bys firm: gen year=_n
gen group=ceil(_n/(_N/50)) 

// Create program to generate sample from data-generating process:
capture program drop panelframe
program define panelframe
	qui sort firm year
	qui gen fe=rnormal() if year==1
	qui replace fe=fe[_n-1] if fe==.
	qui gen x=rnormal()+sqrt(.5)*fe
	qui gen x2=rnormal()
    qui gen u=exp(rnormal())+sqrt(.5)*fe
	qui gen y=x+u
	qui replace y=10000 if _n==_N
	qui sort firm year
	drop u fe
end

panelframe

reghdfe y x x2, absorb(firm year) cluster(group)
robreg m y i.year x x2, ivar(firm) cluster(group) eff(95)
robhdfe m y x x2, absorb(firm year) cluster(group) eff(95) 
robhdfe m y x x2, absorb(firm year) cluster(group year) eff(95)

drop if _n<10
robhdfe m y x x2, absorb(firm year) cluster(group) eff(95) 
robhdfe m y x x2, absorb(firm year) cluster(group) eff(95) keepsin



clear all
sysuse auto, clear
reghdfe price weight length, absorb(rep78) 
robhdfe m price weight length, absorb(rep78) eff(95) 

robreg m price weight length, ivar(rep78) cluster(rep78) eff(95) 
robhdfe m price weight length, absorb(rep78) cluster(rep78) eff(95) 



clear all
webuse nlswork, clear
reghdfe ln_w grade age ttl_exp tenure not_smsa south , absorb(idcode year) cluster(idcode) 
robhdfe m ln_w grade age ttl_exp tenure not_smsa south, absorb(idcode year) cluster(idcode) eff(95) 
