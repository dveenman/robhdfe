//////////////////////////////////////////////////
// Test based on simulated data:
//////////////////////////////////////////////////
clear all
set seed 1234567

local firms 50000
local periods 100
local groups 50
local obs=`firms'*`periods'
set obs `obs'

gen n=_n
gen firm=ceil(_n/`periods') 
bys firm: gen period=_n
gen group=ceil(_n/(_N/`groups')) 

// Create program to generate sample from data-generating process:
capture program drop panelframe
program define panelframe
	qui sort firm period
	qui gen fe=rnormal() if period==1
	qui replace fe=fe[_n-1] if fe==.
	qui gen x=rnormal()+sqrt(.5)*fe
	qui gen x2=rnormal()
    qui gen u=exp(rnormal())+sqrt(.5)*fe // long-tailed error distibution
	qui gen y=x+u
	qui sort firm period
	drop u fe
end

panelframe

timer clear
timer on 1
robreg m y i.period x x2, ivar(firm) cluster(group) eff(95)
timer off 1

timer on 2
robhdfe m y x x2, absorb(firm period) cluster(group) eff(95) 
timer off 2

timer on 3
robhdfe m y x x2, absorb(firm period) cluster(group) eff(95) julia
timer off 3

// Time for robreg:
timer list 1

// Time for robhdfe:
timer list 2

// Time for robhdfe/julia:
timer list 3

// Effect of singleton observation:
sum period
drop if _n<r(max) 
robhdfe m y x x2, absorb(firm period) cluster(group) eff(95) 
robhdfe m y x x2, absorb(firm period) cluster(group) eff(95) keepsin


//////////////////////////////////////////////////
// Test based on external data:
//////////////////////////////////////////////////
clear all
sysuse auto, clear
reghdfe price weight length, absorb(rep78) 
robhdfe m price weight length, absorb(rep78) eff(95) 

robreg m price weight length, ivar(rep78) cluster(rep78) eff(95) 
robhdfe m price weight length, absorb(rep78) cluster(rep78) eff(95) 

clear all
webuse nlswork, clear
reghdfe ln_w grade age ttl_exp tenure not_smsa south, absorb(idcode year) cluster(idcode) 
robhdfe m ln_w grade age ttl_exp tenure not_smsa south, absorb(idcode year) cluster(idcode) eff(95) 

reghdfe ln_w grade age ttl_exp tenure not_smsa south, absorb(idcode year) cluster(idcode) resid
robhdfe m ln_w grade age ttl_exp tenure not_smsa south if _reghdfe_resid!=., absorb(idcode year) cluster(idcode) eff(95) 
robreg m ln_w i.year grade age ttl_exp tenure not_smsa south if _reghdfe_resid!=., ivar(idcode) cluster(idcode) eff(95) 
