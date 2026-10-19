# About
The code in this repository was used to analyze TESS observations of [CALSPEC](https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec) stars.

## Who is this for?
People trying to reproduce the results in my paper [LINK WILL BE HERE, eventually]. If you're trying to do an actual independent analysis, you should probably use [lightkurve](https://lightkurve.github.io/lightkurve/).

## Usage Instructions
You can either run your analysis in the TJCS folder or pip install that directory.

Stars are intended to be analyzed on an individual basis. You input a TIC ID and it automagically runs. There are a few helper commands to make plots and see periodograms.
```
from star import Star

s = Star(352817378)
s.lcs[0].make_plots()
```
