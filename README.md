# NMRQuanTR


**NMRQuanTR** is an R package for calibration-curve-based metabolite
quantification in NMR data. It uses data from 2D HSQC NMR spectra processed with
CCPN (ref.) to first derive calibration curves and then use those curves to
quantify the molecules in the mixture.

### Installation

To get the code either clone this repository:

```git clone http://github.com/ecaamano/NMRQuanTR```

Or download it in a zip file by clicking download button  on the right side of the page.

For running this package you will need to have R programming environment v.3.2.2 or newer.
It might work with older versions but it has only been tested with R-3.2.2.

Package requirements include:

* ggplot2 - for plotting
* gridExtra - for plotting
* testthat - for running tests

All these packages will be installed if you run tests first as explained below.

### Using the package

The package contains three R scripts.

* lib.R - functions that perform calculations and plotting.
* main.R - pipeline functions.
* pipelinesCL.R - an R script for running the pipelines from the command line without opening R.

The easiest way to start using the code is to run the test suite that will 
install the required packages and make sure that everything is working.

To run the test suite navigate to the NMRquanTR/tests folder and run the following line
in the command line:

```Rscript run_tests.R```

Wait for the package installations to finish. If no errors occured NMRQuanTR is ready
to be used.

...
