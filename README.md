# stata-cyreg
Implements "Efficient tests of stock return predictability" Campbell Yogo (2006)

The program computes the 90% confidence band using the CY method for the regression that corresponds for ``reg r L.x``.
```
tsset date
cyreg r L.x, tablepath(table.dta)
```

The option ``nograph`` (or just ``nog``) allows to not display the graph.
