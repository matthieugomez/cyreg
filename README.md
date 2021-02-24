Implements "Efficient tests of stock return predictability" Campbell Yogo (2006) in Stata

The program computes the 90% confidence band using the Campbell Yogo method.

The syntax is as follows:
```
tsset date
cyreg r L.x
```

Options:
- (optional) `nograph` (or just `nog`) to not display the graph
- (optional) `lag` to impose number of lags (otherwise chosen via BIC criterion)
