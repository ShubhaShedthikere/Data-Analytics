#Tester Code

import numpy as np
from Generate_Input_Files import Generate_Files 
import sequence_partial_match_detect as sd
#import Exons as ebin


#------------------------------------------------------------------------------   
#read the chromosome string file and generate bwt files
bwtLastColFileName, bwtRefFileName=Generate_Files("sampleChr.txt")
print "-------------------bwt col and reference files generated-----------------"

#------------------------------------------------------------------------------ 

bwtLastColFileName="sampleBWTCol.txt"
bwtRefFileName="sampleBWTRef.txt"
numOfAllowedMismatches=1
readMatchPositions=sd.PatternApproxMatch("TAG",bwtLastColFileName,bwtRefFileName,numOfAllowedMismatches)
print "position=",positions


#dctionary with start and end values of the bins
exonsBinPositions={10:15, 20:25,29:32,: 27:45}


  