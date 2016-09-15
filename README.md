# MSS-MLE-Fluctuation-Test-Python2.7
About: MSS-Maximum Likelihood test for fluctuation analysis implemented for python2.7. Also calculates 95% confidence intervals.

Inputs: 
  -f: Colony Counts in a text file. For each plate counted, enter the number of counted colonies as a single line entry
      Example:
      23
      45
      1
      0
      0
      156
      1501
      .
      .
      .
  -n: Number of cells plated. Input an integer amount > 0. For the test to be valid, the number of cells plated must be the same (approximately) for each experiment.
  -f: Percentage of culture plated. Input a decimal number 0<x<=1. 
  
Outputs:
  Adjusted M = The number of mutants per culture as estimated by the MSS-Maximum Likelihood method
  Mutation Rate = The probability of mutation per cell per division or generation
  95% Confidence Interval = Confidence Interval as calculated in 
        W.A. Rosche, P.L. Foster
        Determining mutation rates in bacterial populations
        Methods, 20 (2000), pp. 4–17
        
        
Example Usage:

    python2.7 mss-mle_calc.py -f test_data.txt -n 100000 -p 0.15
    Optimization terminated successfully.
         Current function value: -0.000000
         Iterations: 12
         Function evaluations: 24

    Adjusted M =  19.3514954238 
    Mutation Rate =  0.000134134344926 
    95% Confidence Interval = [ 0.000156447851887  ,  0.000113118237358  ]
  
