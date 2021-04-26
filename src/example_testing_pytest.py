# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 10:49:08 2020

@author: linigodelacruz
"""

def add(a, b):
    return a + b

## This code contains one genuine function and a test function. 
#pytest finds any functions beginning with test_ and treats them as tests.

def test_add():
    assert add(2, 3) == 5
    assert add('space', 'ship') == 'spaceship'