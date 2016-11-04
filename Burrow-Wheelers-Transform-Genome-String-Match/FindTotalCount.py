
import numpy as np


totCountEach=np.zeros(4)
mapDict=dict({"A":0,"C":1,"G":2,"T":3})
Counts=np.zeros(4)
ReadMatrixSave=[]
linecount=0


delta=100

print "Start-----"
with open("chrX_last_col.txt") as fInp, open("ReadMatrixSave.txt",'w') as fOut1:
    for line in fInp:
        
        totCountEach[mapDict["A"]]=totCountEach[mapDict["A"]]+line.count("A")
        totCountEach[mapDict["C"]]=totCountEach[mapDict["C"]]+line.count("C")
        totCountEach[mapDict["G"]]=totCountEach[mapDict["G"]]+line.count("G")
        totCountEach[mapDict["T"]]=totCountEach[mapDict["T"]]+line.count("T")
        
        if linecount%delta==0:
            print>>fOut1,totCountEach
          
        linecount=linecount+1   
       
                
fhand=open("totEachCountFile.txt",'w')
print>>fhand,totCountEach   
fhand.close() 

print "totCountEach=",totCountEach

#print ReadMatrixSave