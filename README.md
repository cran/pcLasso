Bug fixes:  

- `predict.pcLasso` now works when `family = “binomial”` (previously, the intercept term was being added in an incorrect manner).
- Previously, `standardize = TRUE` scaled the `beta` coefficients and intercept `a0` incorrectly. This has been fixed.
- `pcLasso` now generates lambda values for the objective function RSS/(2n) + penalty, instead of that for RSS/2 + penalty.