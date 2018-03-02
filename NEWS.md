# chemmodlab 1.0.1

## Bug fixes

* `CombineSplits()` when ppv was used as the model performance measure, it at times 
  evaluated to NA because there were no predicted positives for a model. `Performance()`
  may be used to find these cases.  The ppv for these splits is now imputed using the mean
  of the other splits that are not NA.  If there are no splits that are not NA, 
  `CombineSplits()` returns an error

# chemmodlab 1.0.0

First feature complete release of chemmodlab on CRAN