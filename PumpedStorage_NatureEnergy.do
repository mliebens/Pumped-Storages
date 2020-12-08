clear
use "PumpedStorage_NatureEnergy.dta" // Import dataset

** Create variables for the common market DE/AT	(MWh):
** Electricity production by source
	gen bio_deat = bio_de+bio_at
	gen fosgas_deat = fosgas_de + fosgas_at
	gen foshardcoal_deat = foshardcoal_de + foshardcoal_at
	gen lignite_deat = lignite_de // only available for DE, AT has no lignite electricity
	gen fosoil_deat = fosoil_de + fosoil_at
	gen nuclear_deat = nuclear_de // only available for DE, AT has no nuclear electricity
	gen geo_deat = geo_de + geo_at
	gen hrorap_deat = hrorap_de + hrorap_at
	gen hwr_deat = hwr_de + hwr_at
	gen other_deat = other_de + other_at
	gen otherren_deat = otherren_de + otherren_at
	gen solar_deat = solar_de + solar_at
	gen waste_deat = waste_de + waste_at
	gen windon_deat = windon_de + windon_at
	gen windoff_deat = windoff_de // only available for DE, not for AT

** Load (MWh)
	egen load = rowtotal(*_deat*) // Load in DE/AT
	egen load_de = rowtotal(*_de) // Load in DE
	
** Variable: feed-in from intermittent renewables (MWh)
	gen res = windoff_deat + windon_deat + solar_deat

** Wind (MWh)
	gen wind = windon_deat + windoff_deat
	gen solar = solar_deat

** Temperature^2
	gen tempavg2 = tempavg^2
	
** Generate seasonal indicator variables	
	gen year = year(date)
	gen month = month(date)
	gen week = week(date)
	gen dow = dow(date)
	gen day = day(date)

** Create seasonal fixed effects	
	tab year, gen(y_fe)
	tab month, gen(m_fe)
	tab week, gen(w_fe)
	tab dow, gen(d_fe)
	tab hour, gen(h_fe)
	
** Profits	(industry profit from spot market)
	gen profit = (hps_gen_at - hps_con_at)* price 
	sum profit*
	scalar minprof = r(min)*(-1)
	gen logPi = log(profit + minprof + 1.1)
	sum logPi

** Peak and off-peak prices
	gen ppeak = price if hour==8 | hour==9 | hour==10 | hour==11 | hour==12 | hour==13 | hour==14 | hour==15 | hour==16 | hour==17 | hour==18 | hour==19 | hour==20
	gen poffpeak = price if hour==1 | hour==2 | hour==3 | hour==4 | hour==5 | hour==6 | hour==7 | hour==21 | hour==22 | hour==23 | hour==24

** Price variance	
	egen avgpriceperday = mean(price), by(date)
	sort date hour
	gen pricevariance = price-avgpriceperday	
		
** Declare data to be time-series data	
	sort date hour
	egen hourid = group(date hour)
	tsset hourid
	
** Save sample statistics for later use
	local var "wind solar price res p_co2"
	foreach x in `var' {
		sum `x',d
		scalar `x'mean = r(mean)
		scalar `x'min = r(min)
		scalar `x'p25 = r(p25)
		scalar `x'p50 = r(p50)
		scalar `x'p75 = r(p75)
		scalar `x'max = r(max)
		}
				
** Label variables
	label variable logPi "$log(\pi)$"
	label variable price "Price"
	label variable res "RE forecast"
	label variable load "Load"
	label variable tempavg "Temp"
	label variable tempavg2 "Temp^2"
	label variable p_coal "P_coal"
	label variable p_gas "P_gas" 
	label variable p_co2 "P_co2"
	label variable sunshineavg "Sunshine"
	label variable sunshineavg2 "Sunshine^2"
	label variable windspavg "Windspeed"
	label variable windspavg2 "Windspeed^2"	
	
** Supplementary Table 1: Sample statistics
	sum profit price res load tempavg p_coal p_gas p_co2 windspavg sunshineavg if e(sample)
	
** Supplementary Table 2: Correlations
	pwcorr logPi res load tempavg p_coal p_gas p_co2

** Supplementary Figure 1: Wholesale electricity price & RE feed-in	
	sort res
	graph twoway (scatter price res, msize(tiny) msymbol(smcircle) mcolor(gray)) (lfit price res, lwidth(thin) lcolor(black)  graphregion(color(white)) bgcolor(white) ytitle("Price (€/MWh)") xtitle("RE feed-in (MWh)") legend(off)) 

** Supplementary Figure 2: Feed-in profiles of wind and solar in DE/AT	
	sort hour
	egen windbyhour = mean(wind), by(hour)
	egen solarbyhour = mean(solar), by(hour)
	graph twoway (line windbyhour hour, lcolor(black) lpattern(dash) yaxis(1)) (line solarbyhour hour, lcolor(black) lpattern(longdash) yaxis(1) graphregion(color(white)) bgcolor(white) ytitle("MWh") legend(on)) 
	
** Table 1 (part 1): causal chain via price
	** Normal 2SLS regression (SE not adjusted for first-order serial correlation)
	ivreg2 logPi (price = res) load tempavg tempavg2 p_coal p_gas p_co2 *_fe*, robust first savefirst savefprefix(alpha)
		scalar b_price = _b[price]
		estimates restore alphaprice
		scalar b_res = _b[res]
		di "causal chain = " b_res * b_price * resmean
	** 2SLS regression using Newey-West standard errors (SE robust to heteroskedasticity and allow for first-order serial correlation)
	newey2 logPi (price = res) load tempavg tempavg2 p_coal p_gas p_co2 *_fe*, lag(1) first force
		di "causal chain = " .0127658 * -.0013314 * resmean

** Table 1 (part 2): causal chain via price variance
	** Normal 2SLS regression (SE not adjusted for first-order serial correlation)
	eststo regVariance: ivreg2 logPi (pricevariance = res) load tempavg tempavg2 p_coal p_gas p_co2 *_fe*, robust first savefirst savefprefix(alpha)
		scalar b_pricevar = _b[pricevariance]
		estimates restore alphapricevariance
		scalar b_res = _b[res]
		di "causal chain = " b_res * b_pricevar * resmean
	** 2SLS regression using Newey-West standard errors (SE robust to heteroskedasticity and allow for first-order serial correlation)
	newey2 logPi (pricevariance = res) load tempavg tempavg2 p_coal p_gas p_co2 *_fe*, lag(1) first force	
		di "causal chain = " -.0003912 * .0434444 * resmean
	
** Supplementary Table 3: Direct effect of RE on pump-storage profitability
	tsset hourid
	** Normal OLS fixed-effects regression (SE not adjusted for first-order serial correlation)
	reg logPi res load tempavg tempavg2 p_coal p_gas p_co2 *_fe*, robust
	** regression using Newey-West standard errors (SE robust to heteroskedasticity and allow for first-order serial correlation)
	newey logPi res load tempavg tempavg2 p_coal p_gas p_co2 *_fe*, lag(1) force
		di "causal chain = " _b[res]*resmean

** Table 2: Two-stage estimates of the effect of carbon pricing		
	** Normal 2SLS regression (SE not adjusted for first-order serial correlation)
	ivreg2 logPi (price = res p_co2) load tempavg tempavg2 p_coal p_gas *_fe*, robust first savefirst savefprefix(beta)
	** 2SLS regression using Newey-West standard errors (SE robust to heteroskedasticity and allow for first-order serial correlation)
	newey2 logPi (price = res p_co2) load tempavg tempavg2 p_coal p_gas *_fe*, lag(1) first force
		** Effect of mean RES on profits
		di "causal chain = " -.0013314 * .0125818 * resmean
		** Effect of mean P_CO2 on profits
		di "causal chain = " 1.00729 * .0125818 * p_co2mean

** Supplementary Table 4: Non-constantmarginal effect of RE
	gen res2 = res^2
	gen solar2 = solar^2
	gen wind2 = wind^2
	// Non-Linear with RES
	** First stage:
	reg price c.res##c.res load tempavg tempavg2 p_coal p_gas p_co2 *_fe*, robust
	margins, dydx(res) atmeans
	scalar b_res_nl = el(r(b),1,1)
	** Second stage:
	eststo sq_reg3: ivreg2 logPi (price = res res2) load tempavg tempavg2 p_coal p_gas p_co2 *_fe*, robust first savefirst savefprefix(alpha)
	scalar b_price_nl1 = _b[price]
		** Effect of mean RES on price
		di b_res_nl * resmean
		** Effect of mean RES on profits
		di "causal chain = " b_price_nl1 * b_res_nl * resmean
		
** Supplementary Table 5: Separate effects of wind and solar power	
	ivreg2 logPi (price = wind solar) load tempavg tempavg2 p_coal p_gas p_co2 *_fe*, robust first savefirst savefprefix(alpha)
	newey2 logPi (price = wind solar) load tempavg tempavg2 p_coal p_gas p_co2 *_fe*, lag(1) first  force
	 di -.0013075*windmean
	 di -.0015632*solarmean
	 di .0133445*-.0013075*windmean
	 di .0133445*-.0015632*solarmean

* Supplementary Table 6: Effects for different values of RE	
	di  .0127658 * -.0013314 * resmin
	di  .0127658 * -.0013314 * resp25
	di  .0127658 * -.0013314 * resmean
	di  .0127658 * -.0013314 * resp50
	di  .0127658 * -.0013314 * resp75
	di  .0127658 * -.0013314 * resmax
					
** Figure 2: Production & consumption of pump storages by hour of day		
clear
use "PumpedStorage_NatureEnergy.dta" // Import dataset
	gen  hps_netprod_at = hps_gen_at-hps_con_at
	collapse (mean) hps_gen_at hps_con_at price, by(hour)
	twoway line hps_gen_at hps_con_at hour, yaxis(1)   || line price hour, yaxis(2) legend(order(1 "Production" 2 "Consumption" 3 "Spot price")) graphregion(color(white)) bgcolor(white) ytitle("MWh", axis(1)) ytitle("€/MWh", axis(2)) xlabel(1(1)24) xtitle("")

** Figure 3: Carbon price (EUR/tCO2)
clear
use "PumpedStorage_NatureEnergy.dta" // Import dataset
** Collapse data to daily frequency
	collapse (mean) p_co2, by(date)
	twoway line p_co2 date, yline(8.7, lpattern(dash) lcolor(gray)) graphregion(color(white)) bgcolor(white)
	
