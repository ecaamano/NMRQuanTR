
library(testthat)

#source('tests/test_lib.R')
source('lib.R')
source('main.R')

test_results = test_dir('tests', reporter = 'summary')
