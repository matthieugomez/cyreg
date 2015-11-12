# stata-cyreg
Implement "Efficient tests of stock return predictability" Campbell Yogo (2006)

```
tsset date
cyreg r L.x, tablepath(table.dta)
```

The option ``nograph`` (or just ``nog``) allows to not display the graph.
