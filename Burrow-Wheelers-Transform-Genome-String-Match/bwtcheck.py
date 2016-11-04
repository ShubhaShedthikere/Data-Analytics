# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 19:20:11 2016

@author: shubha
"""

fhand=open("chrMapFragmentPart1328")
bwtref=fhand.read().split("\n")[:-1]
print bwtref
print type(bwtref)
print len(bwtref)