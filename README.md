# cyreg

This Stata package implements "Efficient tests of stock return predictability" (Campbell and Yogo, 2006). It computes 90% Bonferroni confidence bands for predictive regressions where the predictor may be persistent (near unit root).

## Syntax
```stata
tsset date
cyreg r L.x [, options]
```

## Options
- `lag(#)` — impose the number of lags for the AR(p) specification. When left unspecified, the optimal lag length is chosen via BIC criterion.
- `maxlag(#)` — set the maximum lag length considered by BIC (default: 4). Campbell and Yogo (2006) recommend 4 for annual data, 6 for quarterly, and 8 for monthly.
- `nograph` — suppress the confidence region graph.

## Installation
```stata
net install cyreg, from("https://raw.githubusercontent.com/matthieugomez/cyreg/main/")
```

## Reference
Campbell, John Y., and Motohiro Yogo. 2006. "Efficient tests of stock return predictability." *Journal of Finance* 61(2): 715-747.
