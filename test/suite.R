# Test Suite
library(testthat)
library(lme4)
library(mvtnorm)
library(MCMCpack)

# Set up
  # setwd("/Users/Edoardo/DriveUni/MasterThesis/BayesianThesis")
  # ! WD: should be project location (main folder)
  code_path <- "./R"    #path to folder containing all the package scripts
  test_path <- "./test" #path to folder containing different test-name.R scripts

# Autotest
  auto_test(code_path, test_path, reporter = TapReporter)

# Reporter choices:
# - LocationReporter: prints the location of every expectation and error
# - TapReporter: "ok 1 TEST: different comparison sings" (good for mismatch report)

# Output Interpretation:
# - Mismatch: if a mismatch between "tocheck" and "benchmark" appears, it means that the changes 
#   applied to the function of interest (e.g. statcheck) are providing different results from the 
#   stable function used in the package at the moment. When a mismatch ooccurs, testthat output
#   shows a difference. The first number is always the number of extraction obtained with the
#   function you are working on and the second one is always the number of extractions obtained
#   with the stable function. Using the LocationReporter option allows to see directly in which
#   script and line is the test that is failing.
# - Running tests: if the last line in the console is "Rerunning tests:  test-name.R"
#   it means that R is still (re)running a test.