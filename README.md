# Better_Summaries
Better summary() functions for lm, glm, gee, and geese objects in R

Description:<br />
The summary() function may not give you a satisfactory summary table. Use these functions instead to obtain
the summary table portion of the summary() function that includes confidence intervals (by normal approximation 
or, optionally for glms, profile likelihood), p-values whose number of significant digits is adustable, and 
appropriately exponentiated and labeled results for logistic and poisson regressions.

Usage:<br />
lm.mysummary( fit, alpha=.05, dig=3, p.dig=4, ci="profile.lik")<br />
glm.mysummary(fit, alpha=.05, dig=3, p.dig=4, ci="profile.lik")<br />
geese.mysummary (fit, alpha=.05, dig=3, p.dig=4)<br />
gee.mysummary(fit, alpha=.05, dig=3, p.dig=4)<br />
<br />

Arguments:<br />
```
fit   - lm model object for lm.mysummary, glm model object for glm.mysummary, geese model object for 
        geese.mysummary, and gee model object for gee.mysummary
alpha - significance level
ci    - for use in glm.mysummary(), either "profile.lik" for a profile likelihood- based confidence interval 
        (preferred) or "normal.approx" for a normal approximation- based confidence interval.
dig   - number of significant digits to report for all values except the p-value
p.dig - number of significant digits to report for the p-value
```

Aside: Note that, unless family is independence or exchangeable, gee() and geese() both need the data set 
to be sorted by time within subject, and missing measurements for a subject should be still be given a row 
in the data set (with NAs (but no NAs in the id or time columns)). Also, the geese function assumes that no 
values are missing. Hence if one measurement is missing on a subject, then this subject has to be excluded 
from the analysis.
