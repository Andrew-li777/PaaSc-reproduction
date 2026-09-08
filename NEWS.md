# PaaSc 1.0.1

Bug fix release.

## Fixed

* `computeMCA()` aborted with

  ```
  Error in (function (cond) : error in evaluating the argument 'x' in selecting a
  method for function 'as.matrix': The `slot` argument of `GetAssayData()` was
  deprecated in SeuratObject 5.0.0 and is now defunct.
  ```

  on any installation with SeuratObject >= 5.4.0. `computeMCA()` delegated to
  `CelliD::RunMCA()`'s Seurat method, which still fetches the expression matrix with
  `GetAssayData(X, slot = slot)`. PaaSc now reads the matrix itself - with `layer` on
  SeuratObject >= 5.0.0 and `slot` on earlier versions - and passes it to CelliD's matrix
  method. The reduction is filled in exactly as CelliD does (dimension names `mca_1...`,
  `sqrt()` applied to `stdev`, `mca.flag` set), so results and downstream functions,
  including `computeScore(weight = TRUE)`, are unchanged. Seurat 4 installations are
  unaffected.

## Changed

* Declared the dependencies the code actually uses. `NormalizeData()`,
  `DefaultAssay()` and `CreateDimReducObject()` are exported by Seurat/SeuratObject, not
  by CelliD, so they were previously resolved from the search path and `computeMCA()`
  only worked if the user had attached Seurat themselves. They are now imported in
  NAMESPACE, and Seurat, SeuratObject, stats and utils are listed under Imports.
