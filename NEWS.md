# chemmodlab 1.1.0.9000

# chemmodlab 1.1.0

## Major changes

* This version of `chemmodlab` has all of the functionality mentioned in the following paper: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0309-4

* Added function `ApplicabilityDomain` which uses a Hotelling T2 control chart to identify outliers in an external data set for which predictions are desired, hence identifying observations whose model predictions may be considered extrapolations.

* `ModelTrain` now has a S3 method which accepts molecule objects that are created by the package rcdk. rcdk supports most of the widely used chemical file formats (SMILES, SDF, InChI, Mol2, CML, etc). When molecules are provided to `ModelTrain`, the names of predefined descriptor sets and/or fingerprints must also be provided.y 

# chemmodlab 1.0.1

## Bug fixes

* `CombineSplits()` when ppv was used as the model performance measure, it at times 
  evaluated to NA because there were no predicted positives for a model. `Performance()`
  may be used to find these cases.  The ppv for these splits is now imputed using the mean
  of the other splits that are not NA.  If there are no splits that are not NA, 
  `CombineSplits()` returns an error

# chemmodlab 1.0.0

First feature complete release of chemmodlab on CRAN
