#Tester Code

import numpy as np
#import sequence_detect


#------------------------------------------------------------------------------   
#read the chromosome string file 
def ReadChr(filename):    
    fhand=open(filename)  
    fstring=fhand.read()
    chrString=fstring.replace("\n","")
    return chrString

#------------------------------------------------------------------------------
#Obtain BWT column and the index of the sorted rotated string in the original
# chromosome
def shiftRightOnce(inpString):
    return inpString[1:]+inpString[0]
    
        
    

def GenerateBWT(inpChr):
    posiDict={inpChr:0}
    shiftedStr=inpChr
    for i in range(len(inpChr)):
        posiDict[shiftedStr]=i
        shiftedStr=shiftRightOnce(shiftedStr)
        
    sortedListFinal=sorted(posiDict)
    #sortedListFinal=sortedList[1:]
    #sortedListFinal.append(sortedList[0])  
    #print sortedListFinal              
    bwtLastCol= [ key[-1]for (key, value) in sorted(posiDict.items())]
   # firstItem=bwtLastCol.pop(0)
   # bwtLastCol.append(firstItem)
    bwtRefIndex= [value for (key, value) in sorted(posiDict.items())]
    #firstItem=bwtRefIndex.pop(0)
    #bwtRefIndex.append(firstItem)
    bwtLastColString=''.join(bwtLastCol)
    text_file = open("sampleBWTCol.txt", "w")
    text_file.write(bwtLastColString)
    text_file.close()
    text_file = open("sampleBWTRef.txt", "w")
    for item in bwtRefIndex:
        print>>text_file, item
    
    text_file = open("sampleBWTMatrix.txt", "w")
    for item in sortedListFinal:
        print>>text_file, item    
  
    return   "sampleBWTCol.txt","sampleBWTRef.txt"
  
    
#------------------------------------------------------------------------------
def Generate_Files(ChrFileName):
    inpStr=ReadChr(ChrFileName)    
    f1,f2=GenerateBWT(inpStr)
    return f1,f2
   

bwtLastColFileName, bwtRefFileName=Generate_Files("sampleChr.txt")



  