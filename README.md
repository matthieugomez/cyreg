# stata-cyreg
Implements "Efficient tests of stock return predictability" Campbell Yogo (2006)

The program computes the 90% confidence band using the CY method.

Syntax
```
tsset date
cyreg r L.x
```

Options:
- (optional) `nograph` (or just `nog`) to not display the graph
- (optional) `lag` to impose number of lags (otherwise chosen via BIC criterion)
