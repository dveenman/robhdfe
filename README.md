# robhdfe

A Stata program to combine robust regression estimation (Huber M) with high-dimensional fixed effects and clustered standard errors. The program accompanies the Gassen & Veenman (2026) study on "Estimation Precision and Robust Inference in Archival Research'' (https://ssrn.com/abstract=4975569). The accompanying R package `ferols` can be found at https://github.com/joachim-gassen/ferols.


---

Installation:
```
net install robhdfe, replace from(https://raw.githubusercontent.com/dveenman/robhdfe/main/)
```

The program requires `moremata`, `reghdfe`, and `hdfe` to be installed in Stata:
```
ssc inst moremata, replace
ssc inst hdfe, replace
ssc inst reghdfe, replace 
```

---

Syntax:
