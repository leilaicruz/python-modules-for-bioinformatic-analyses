# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 10:34:04 2020

@author: linigodelacruz
"""
## Testing ina nutshell 
## In software tests, expected results are compared with observed results in order to establish accuracy:
def fahrenheit_to_celsius(temp_f):
    temp_c = (temp_f - 32.0) * (5.0/9.0)
    return temp_c



# This is the test function: `assert` raises an error if something
# is wrong.
def test_fahrenheit_to_celsius():
    temp_c = fahrenheit_to_celsius(temp_f=100.0)
    expected_result = 37.777777
    assert abs(temp_c - expected_result) < 1.0e-6
    
    
## Defensive programming 
#Assume that mistakes will happen and introduce guards against them.

#Use assertions for things you believe will/should never happen.

#Use exceptions for anomalous or exceptional conditions requiring special processing.

def kelvin_to_celsius(temp_k):
    """
    Converts temperature in Kelvin
    to Celsius.
    """
    assert temp_k >= 0.0, "ERROR: negative T_K"
    temp_c = temp_k + 273.15
    return temp_c
