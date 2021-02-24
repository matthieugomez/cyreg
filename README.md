# cyreg

The program implements "Efficient tests of stock return predictability" Campbell Yogo (2006) in Stata

It returns the 90% confidence band using the Campbell Yogo method.

The syntax is:
```
tsset date
cyreg r L.x
```

Options:
- `nograph` to suppress the graph
- `lag` to impose number of lags (otherwise chosen via BIC criterion)


# Installation

```
net install cyreg from("https://raw.githubusercontent.com/matthieugomez/cyreg/master/")
```
If you have a version of Stata < 13, you need to install it manually

Click the "Download ZIP" button in the right column to download a zipfile.

Extract it into a folder (e.g. ~/SOMEFOLDER)

Run
```
cap ado uninstall cyreg
net install cyreg, from("~/SOMEFOLDER")
```