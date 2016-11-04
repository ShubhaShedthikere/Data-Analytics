import numpy as np
import Band as bd
import re


#-------------------------------------------------------------------------------------
#Input Parameters and other input information

exonsBinPositions={10:90, 120:185,229:232}
readsFileName="reads.txt"
chromosomeFileName=""
bwtLastColFileName="chrX_last_col.txt"
bwtRefFileName="chrX_map.txt"
numOfAllowedMismatches=0


#-------------------------------------------------------------------------------------
# Global Variables

ACGTBin=[]
ReadMatrixSave=[]
totCountEach=np.zeros(4).astype(int)
firstColStartIndex=np.zeros(4).astype(int)
exonBinReadsCount={}
mapDict=dict({"A":0,"C":1,"G":2,"T":3})



def ReadBWT(bwtLastColFileName,bwtRefFileName):
    global bwtLastCol
    global ACGTBin
    global ReadMatrixSave
    global totCountEach
    global bwtRef
    global firstColStartIndex
    
    
    # Load the binary version of the BWT column
   
    
    bwtCol=[]
    bwtLastColFileName="chrX_last_col.txt"
    linecount=0
    n=20
    delta=100
    print "Loading BWT last column.............."
    with open(bwtLastColFileName) as fhand:
        for line in fhand:
              text=line
              totCountEach[mapDict["A"]]=totCountEach[mapDict["A"]]+line.count("A")
              totCountEach[mapDict["C"]]=totCountEach[mapDict["C"]]+line.count("C")
              totCountEach[mapDict["G"]]=totCountEach[mapDict["G"]]+line.count("G")
              totCountEach[mapDict["T"]]=totCountEach[mapDict["T"]]+line.count("T")
              if linecount%delta==0:
                  ReadMatrixSave.append(totCountEach)
                  
              bwtCol.append(text)
              linecount=linecount+1
              
    print "totCountEach",totCountEach          
    print "Num of elements=",len(bwtCol)
    
    
ReadBWT(bwtLastColFileName,bwtRefFileName)    
   
    
             
    
