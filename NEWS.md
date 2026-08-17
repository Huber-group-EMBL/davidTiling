# davidTiling 1.53.0

* The pre-built ExpressionSet object, previously available via 
  `data(davidTiling)` has been removed, alongside the bundled CEL files in
  `inst/celfiles`.
  To get the data from now on, please use the new `getDavidTilingData()` 
  function, which will download the raw CEL files from ArrayExpress (or use a
  local cache if they were downloaded previously), and convert them to the
  ExpressionSet object via `tilingArray::readCel2eSet()`.
