# cyreg

The program implements "Efficient tests of stock return predictability" Campbell Yogo (2006) in Stata

It returns the 90% confidence band using the Campbell Yogo method.

The syntax is:
```
tsset date
cyreg r L.x
```

Options:
- (optional) `nograph` to suppress the graph
- (optional) `lag` to impose number of lags (otherwise chosen via BIC criterion)
