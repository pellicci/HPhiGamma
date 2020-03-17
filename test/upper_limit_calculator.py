import ROOT
import math
import os

debug = True

fInput1 = ROOT.TFile(histos/latest_production/MCbackgroundHistogram.root)
fInput2 = ROOT.TFile(histos/latest_production/histo_Signal.root)

alpha = 0.05
backgroundEvents = 20 # b
totalEvents = 21  # n = s + b
s_zero = 5*totalEvents

"""
def factorial(x):
    if x == 0:
       return 1
    else:
       return x * factorial(x-1)
"""


b = backgroundEvents
n = totalEvents
summation = 0.
print "alpha = ",alpha,"  b = ",b,"  n = ",n, "  s0 starting = ",s_zero 

def upper_limit(s):
   
    if debug:
        print ""
        print "s0 = ",s
   
    if s == 0:
        return s

    summation = 0.
    for k in range(1,n):
        exponential = math.exp(-b-s)
        power = (b+s)**k
        fac = math.factorial(k)
        summation += exponential * power * (1/fac) 
        
        if debug:
            print "exp = ",exponential
            print "pow = ",power
            print "factorial = ", fac
    print "summation = ",summation
            
    if summation <= alpha + alpha/1000 and summation >= alpha - alpha/1000:  
        return s

    else:
        s -= 0.5
        return upper_limit(s)

print "s0 = ", upper_limit(s_zero)



    
    


    
