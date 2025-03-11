clear all

* Set parameters
local N = 5 // number of categories
local B = 2^`N' - 1 //number of baskets
local S = 4 //number of segments, this can be determined with k-median analysis



/* create all the possible baskets combinations */
clear
set obs `B'
gen num = _n
forval j = 1/`N' {
    gen V`j' = mod(floor(num/(2^(`j'-1))),2)
}
gen sortkey = ""
forval j=1/`N' {
    replace sortkey = sortkey + string(V`j') 
}

egen elements = rowtotal(V*)
gen group = 0
replace group = 1 if elements == 1
replace group = 2 if elements == 2
replace group = 3 if elements == 3
replace group = 4 if elements == 4
replace group = 5 if elements == 5

gsort group -sortkey
mkmat V1-V5, matrix(Z)



import delimited using "Nth.csv", clear
save "Nth.dta", replace

/*segment analysis*/
//long-run data to uncover segments structure
collapse (sum) bn, by(h b)
tempfile Nh
save `Nh'

clear 
svmat Z
forvalues i = 1/`N' {
    rename Z`i' a`i'
}
gen b = _n
save "Zb.dta", replace

use `Nh', clear
save "Nh.dta", replace    
use "Nh.dta", clear
merge m:1 b using "Zb.dta"
drop _merge
save lnhb, replace
forvalues i = 1/`N' {
    gen new_a`i' = bn * a`i'
}


collapse (sum) new_a*, by(h)

foreach var of varlist new_a* {
    egen std_`var' = std(`var')  
    drop `var'                   
    rename std_`var' `var'       
}

forvalues i = 1/`N' {
    rename new_a`i' a`i'
}
/// the above step is used to compute category counts, which we will use to run cluster analysis to identify customer- segment membership


/*cluster analysis */

cluster kmedians a1 a2 a3 a4 a5, k(4) measure(L1) start(krandom)

rename _clus_1 segment //this is the segment-customer membership
save "membership.dta", replace



gen byte count_seg = 1
collapse (sum) count_seg, by(segment)
list

matrix pi0 = J(1, `S', 0)  


forvalues i = 1/`S' {
sum count_seg if segment == `i'
    matrix pi0[1,`i'] = r(sum)
}


mata: st_matrix("pi0", st_matrix("pi0") :/ sum(st_matrix("pi0")))

matrix list pi0

use "Nth.dta", clear
merge m:1 h using "membership.dta", nogen
drop a1 a2 a3 a4 a5
save "Nth.dta",replace
// a1 a2 ... represents alpha of 5 categories

import delimited using "price.csv", clear
mkmat v1-v5, matrix(price)


import delimited using "prom.csv", clear
mkmat v1-v5, matrix(prom)

mata:

Z = st_matrix("Z")  
B = strtoreal(st_local("B"))
N = strtoreal(st_local("N"))

thetab = J(B, (N*(N-1))/2, 0)
for(i = 1; i <= B; i++) {
    Zi = J(N, 2, 0)
    Zi[.,1] = Z[i,.]'  
    
    Zi_full = Zi * Zi'
    
    
    idx = 1
    for(r = 2; r <= N; r++) {
        for(c = 1; c < r; c++) {
            thetab[i,idx++] = Zi_full[r,c]
        }
    }
}


price = st_matrix("price")
prom = st_matrix("prom")
priceb = Z * price'
promb = Z * prom'


st_matrix("thetab", thetab)
st_matrix("priceb", priceb)
st_matrix("promb", promb)
end


clear
svmat thetab, names(v)
gen b = _n
save "thetab.dta", replace


clear
svmat priceb, names(price)
gen b = _n
reshape long price, i(b) j(t)
save "priceb.dta", replace


clear
svmat promb, names(prom)
gen b = _n
reshape long prom, i(b) j(t)
save "promb.dta", replace


use "Nth.dta", clear
merge m:1 b using "Zb.dta", nogen
merge m:1 b t using "priceb.dta", nogen
merge m:1 b t using "promb.dta", nogen
merge m:1 b using "thetab.dta", nogen
save "snhb_new.dta", replace


duplicates drop h t, force
gen ht = _n
keep h t ht
save "ht.dta",replace
use "snhb_new.dta",clear
merge m:1 h t using "ht.dta", nogen
save "snhb_new.dta", replace



* Define variables to collect
local var_list a1 a2 a3 a4 a5 v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 price prom _cons

* Count total variables 
local total_vars : word count `var_list'

* Create results matrix for all segments plus full sample
matrix results = J(`total_vars', `S', .)
matrix rownames results = `var_list'
matrix colnames results = "Segment1" "Segment2" "Segment3" "Segment4" 

* Create SE matrix
matrix se_results = J(`total_vars', `S', .)
matrix rownames se_results = `var_list'
matrix colnames se_results = "Segment1" "Segment2" "Segment3" "Segment4" 

* First phase: For each segment
forvalues i = 1/`S' {
    preserve 
    keep if segment == `i'
    
    * First-stage model
    xtset ht
    xtpoisson bn a1 a2 a3 a4 a5 v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 price prom, fe
    
    * Store coefficients from first stage
    foreach var of local var_list {
        if "`var'" != "_cons" {
            matrix results[rownumb(results,"`var'"), `i'] = _b[`var']
            matrix se_results[rownumb(se_results,"`var'"), `i'] = _se[`var']
        }
    }
    
    * Second-stage model
    predict fitted_values, nu0
    generate ln_fitted = ln(fitted_values)
    poisson bn ln_fitted
    
    * Store constant from second stage
    matrix results[rownumb(results,"_cons"), `i'] = _b[_cons]
    matrix se_results[rownumb(se_results,"_cons"), `i'] = _se[_cons]
    
    restore
}

* Display results in Stata
matrix list results
matrix list se_results


*** code below is one way to export the results

* Export to Excel using a simple approach. 
putexcel set "Segment_Results.xlsx", replace

* Add headers
putexcel A1 = "Variable"
putexcel B1 = "Segment1"
putexcel C1 = "Segment2"
putexcel D1 = "Segment3"
putexcel E1 = "Segment4"


* Add results row by row
local row = 2
foreach var of local var_list {
    putexcel A`row' = "`var'"
    
    * Add coefficient for Segment 1
    local coef1 = results[rownumb(results,"`var'"), 1]
    local se1 = se_results[rownumb(se_results,"`var'"), 1]
    local t1 = abs(`coef1'/`se1')
    local stars1 = ""
    if `t1' > 2.58 local stars1 = "***"
    else if `t1' > 1.96 local stars1 = "**" 
    else if `t1' > 1.65 local stars1 = "*"
    putexcel B`row' = "`=string(`coef1', "%9.3f")''`stars1''"
    putexcel B`=`row'+1' = "(`=string(`se1', "%9.3f")')"
    
    * Add coefficient for Segment 2
    local coef2 = results[rownumb(results,"`var'"), 2]
    local se2 = se_results[rownumb(se_results,"`var'"), 2]
    local t2 = abs(`coef2'/`se2')
    local stars2 = ""
    if `t2' > 2.58 local stars2 = "***"
    else if `t2' > 1.96 local stars2 = "**" 
    else if `t2' > 1.65 local stars2 = "*"
    putexcel C`row' = "`=string(`coef2', "%9.3f")''`stars2''"
    putexcel C`=`row'+1' = "(`=string(`se2', "%9.3f")')"
    
    * Add coefficient for Segment 3
    local coef3 = results[rownumb(results,"`var'"), 3]
    local se3 = se_results[rownumb(se_results,"`var'"), 3]
    local t3 = abs(`coef3'/`se3')
    local stars3 = ""
    if `t3' > 2.58 local stars3 = "***"
    else if `t3' > 1.96 local stars3 = "**" 
    else if `t3' > 1.65 local stars3 = "*"
    putexcel D`row' = "`=string(`coef3', "%9.3f")''`stars3''"
    putexcel D`=`row'+1' = "(`=string(`se3', "%9.3f")')"
    
    * Add coefficient for Segment 4
    local coef4 = results[rownumb(results,"`var'"), 4]
    local se4 = se_results[rownumb(se_results,"`var'"), 4]
    local t4 = abs(`coef4'/`se4')
    local stars4 = ""
    if `t4' > 2.58 local stars4 = "***"
    else if `t4' > 1.96 local stars4 = "**" 
    else if `t4' > 1.65 local stars4 = "*"
    putexcel E`row' = "`=string(`coef4', "%9.3f")''`stars4''"
    putexcel E`=`row'+1' = "(`=string(`se4', "%9.3f")')"
    
      * Move to next variable (skipping the SE row)
    local row = `row' + 2
}

* Add note about significance
putexcel A`row' = "* p<0.1, ** p<0.05, *** p<0.01"

di "Results exported to Segment_Results.xlsx with all 4 segments"

/// in Segment_Results.xlsx, a1,..,a5 represents alpha of 5 categories, V1...V10 represents the lower_tri elements of theta matrix: theta_21, theta_31, theta_32, theta_41, theta_42, theta_43, theta_51,theta_52, theta_53, theta_54. price, prom represents the price and prom coefficients (beta). _cons is the frequency scale R
