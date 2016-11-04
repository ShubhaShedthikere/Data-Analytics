import numpy as np
import Band as bd
import re


#-------------------------------------------------------------------------------------
#Input Parameters and other input information

exonsBinPositions={149249758:149249869, \
                   149256128:149256424, \
                   149258413:149258581, \
                   149260049:149260214, \
                   149261769:149262008, \
                   149264291:149264401, \
                   149293259:149293555,\
                   149295543:149295711,\
                   149297179:149297344,\
                   149298899:149299138}
readsFileName="reads.txt"
chromosomeFileName=""
#bwtLastColFileName="sampleBWTCol.txt"
bwtLastColFileName="chrX_last_col.txt"
bwtRefFileName="chrX_map.txt"
#bwtRefFileName="sampleBWTRef.txt"
numOfAllowedMismatches=2
numOfCharPerLineBWT=100
numOfLinesSkippedReadSave=100


#-------------------------------------------------------------------------------------
# Global Variables
bwtCol=[]
ACGTBin=[]
ReadMatrixSave=[]
totCountEach=np.zeros(4).astype(int)
firstColStartIndex=np.zeros(4).astype(int)
exonBinReadsCount={}
mapDict=dict({"A":0,"C":1,"G":2,"T":3})

#-------------------------------------------------------------------------------------
# read the BWT last column and create a matrix of cumulative sum
def ReadBWT(bwtLastColFileName,bwtRefFileName):
    global bwtCol
    global ACGTBin
    global ReadMatrixSave
    global totCountEach
    global bwtRef
    global firstColStartIndex
    
    
    # Load the binary version of the BWT column
   
   
    linecount=0
    
    print "Loading BWT last column.............."
    with open(bwtLastColFileName) as fhand:
        for line in fhand:
              text=line
              #print text
              totCountEach[mapDict["A"]]=totCountEach[mapDict["A"]]+line.count("A")
              totCountEach[mapDict["C"]]=totCountEach[mapDict["C"]]+line.count("C")
              totCountEach[mapDict["G"]]=totCountEach[mapDict["G"]]+line.count("G")
              totCountEach[mapDict["T"]]=totCountEach[mapDict["T"]]+line.count("T")
              
              if linecount%numOfLinesSkippedReadSave==0:
                  ReadMatrixSave.append(list(totCountEach))
                  
                  
              bwtCol.append(text.replace("\n",""))
              linecount=linecount+1
              
             
              
              
              
   # print "totCountEach",totCountEach          
   # print "Num of elements=",len(bwtCol)
   # print "ReadMatrixSave=",ReadMatrixSave
    
    
             
    
#--------------------------------------------------------------------------------------
#It takes the start index, end index and the character of a band of bwt columns
# and returns a binary string of that band which has 1 where the character is
# is present, else 0    
    

def OccuranceFetch(startIndex,endIndex,character):
    #print "Inside OccuranceFetch....."
    ACGTFetchStartIndex=startIndex/numOfCharPerLineBWT
    #print "startIndex=",startIndex,"endIndex=",endIndex
    remStart=startIndex%numOfCharPerLineBWT
    #print "remStart",remStart
    ACGTFetchEndIndex=endIndex/numOfCharPerLineBWT
    remEnd=endIndex%numOfCharPerLineBWT
    numOfOccurances=0
    firstIndex=-1
    lastIndex=-1   
    #print "ACGTFetchStartIndex=",ACGTFetchStartIndex
    #print "ACGTFetchEndIndex=",ACGTFetchEndIndex
    
    for i in range(ACGTFetchStartIndex,ACGTFetchEndIndex+1):
        
        #print "i=",i,"bwtCol[i]",bwtCol[i]
        if (ACGTFetchStartIndex==ACGTFetchEndIndex):
            numOfOccurances= bwtCol[i].count(character,remStart,remEnd+1)
            if numOfOccurances!=0:
               firstIndex=i*numOfCharPerLineBWT+bwtCol[i].find(character,remStart,remEnd+1)
               lastIndex=i*numOfCharPerLineBWT+bwtCol[i].rfind(character,remStart,remEnd+1)
        
        else:
              
           if i==ACGTFetchStartIndex:
              numOfOccurances=numOfOccurances+\
                            bwtCol[i].count(character,remStart)
              if numOfOccurances!=0:
                 firstIndex=i*numOfCharPerLineBWT+bwtCol[i].find(character,remStart)
                     
           elif i==ACGTFetchEndIndex:  
                numOfOccurances=numOfOccurances+\
                            bwtCol[i].count(character,0,remEnd+1)
                if firstIndex==-1 and numOfOccurances!=0:
                    firstIndex=i*numOfCharPerLineBWT+bwtCol[i].find(character,0,remEnd+1)          
               
           else:
               numOfOccurances=numOfOccurances+\
                            bwtCol[i].count(character) 
                            
               if firstIndex==-1 and numOfOccurances!=0:
                  firstIndex=i*numOfCharPerLineBWT+bwtCol[i].find(character)                    
         
     
    for i in range(ACGTFetchEndIndex,ACGTFetchStartIndex-1,-1):
               if   (ACGTFetchEndIndex!=ACGTFetchStartIndex):              
                   if i==ACGTFetchEndIndex and bwtCol[i].count(character,0,remEnd+1)!=0 \
                      and lastIndex==-1: 
                         lastIndex=i*numOfCharPerLineBWT+bwtCol[i].rfind(character,0,remEnd+1)
                   elif i== ACGTFetchStartIndex and  bwtCol[i].count(character,remStart)!=0 \
                        and lastIndex==-1:                
                          lastIndex=i*numOfCharPerLineBWT+bwtCol[i].rfind(character,remStart)
                   else: 
                          if lastIndex==-1 and bwtCol[i].count(character)!=0:
                           lastIndex=i*numOfCharPerLineBWT+bwtCol[i].rfind(character)
               
     
    #print "firstIndex=",firstIndex,"lastIndex=",lastIndex,"numOfOccurances",numOfOccurances
    
    return firstIndex,lastIndex,numOfOccurances      
           
          
                 
    
    
#--------------------------------------------------------------------------------------
#Rank query
def RankQuery(character,rowIndex):
    #print "Inside rank query........"
    CharOffset=rowIndex%numOfCharPerLineBWT
    BwtColLineNum=rowIndex/numOfCharPerLineBWT
    LineOffset=BwtColLineNum%numOfLinesSkippedReadSave
    ReadSaveIndex=rowIndex/(numOfLinesSkippedReadSave*numOfCharPerLineBWT)
    #print "LineOffset",LineOffset,"ReadSaveIndex",ReadSaveIndex
    #print "BwtColLineNum",BwtColLineNum,"CharOffset",CharOffset
    #print ReadMatrixSave[ReadSaveIndex]
    if LineOffset==0 :
       if CharOffset==9:
           return ReadMatrixSave[ReadSaveIndex][mapDict[character]]-1
       else: 
           startIndex=rowIndex+1
           endIndex=ReadSaveIndex*numOfCharPerLineBWT*numOfLinesSkippedReadSave+\
                       (numOfCharPerLineBWT-1)
           firstIndex,lastIndex,numOfOccurances=OccuranceFetch(startIndex,endIndex,character) 
      # print "ReadMatrixSave[",ReadSaveIndex,"]",ReadMatrixSave[ReadSaveIndex]  
           return ReadMatrixSave[ReadSaveIndex][mapDict[character]]-1-numOfOccurances 
       
    else:
       startIndex= ReadSaveIndex*numOfCharPerLineBWT*numOfLinesSkippedReadSave+\
                   numOfCharPerLineBWT
       endIndex=rowIndex
       firstIndex,lastIndex,numOfOccurances=OccuranceFetch(startIndex,endIndex,character) 
      # print "ReadMatrixSave[",ReadSaveIndex,"]",ReadMatrixSave[ReadSaveIndex]  
       return ReadMatrixSave[ReadSaveIndex][mapDict[character]]-1+numOfOccurances 
   
    
  

#--------------------------------------------------------------------------------------    
#Returns the row Indices of the band of the corresponding characters in 
#first column    
def FindBand(character, minRank,maxRank):
    startBandIndex=firstColStartIndex[mapDict[character]]+minRank
    return startBandIndex, startBandIndex+(maxRank-minRank)
    
#--------------------------------------------------------------------------------------
# This function returns the a list of band objects which are the 
# valid sub-bands of the given band
def FindSubBands(bandObj,seqChar,maxMismatches):
    bandReturnList=[]
    indexOfSeqChar0=None
    #print "Inside Find Subbands............"
    for char in ['A','C','G','T']:
     #   print "--------------Finding the",char,"-subband of",bandObj.char," ------------------"
      #  print "bandObj.startIndex=",bandObj.startIndex,"bandObj.endIndex=",bandObj.endIndex
        firstIndex,lastIndex,numOfOccurances =OccuranceFetch(bandObj.startIndex,bandObj.endIndex,char)
        indexOfSeqChar0=firstIndex
        indexOfSeqChar1=lastIndex
       
        
        if (indexOfSeqChar0== -1):
            continue
        else:
            #for each if the entry in indexOfSeqChar obtain rank
           # print "The ranks of ",char," with the rowIndex",firstIndex," and ",lastIndex
            minRank=RankQuery(char,indexOfSeqChar0)
            maxRank=RankQuery(char,indexOfSeqChar1)
            #print "minRank=",minRank,"maxRank=",maxRank
            startBandIndex, endBandIndex = FindBand(char,minRank,maxRank)
            totBandMismatches = bandObj.mismatches + (char!=seqChar)
            if totBandMismatches<=maxMismatches:
                newBandObj=bd.Band(char,startBandIndex,endBandIndex,totBandMismatches) 
                bandReturnList.append(newBandObj)
            else:
                continue
    return bandReturnList

#--------------------------------------------------------------------------------------
# find the position in the reference chromosome
def  FindRefPosition(startRowIndex,endRowIndex):
     chrName="chrMapFragmentPart"
     partNumStart=startRowIndex/100000
     remStart=startRowIndex%100000
     remEnd=endRowIndex%100000
     partNumEnd=endRowIndex/100000
     position=[]
     if partNumStart==partNumEnd:
        #print "only one bwf file loaded"
        fileName=chrName+ str(partNumStart)
        fhand=open(fileName)
        bwtref=fhand.read().split("\n")[:-1]
        fhand.close()
        position = bwtref[remStart:remEnd+1]
     return position   
 
def FindRefPosition1(startRowIndex,endRowIndex):
    fhand=open(bwtRefFileName)
    bwtRefInp=fhand.read()
    fhand.close()
    bwtRef=bwtRefInp.split("\n")
   # print bwtRef
    posi=(np.asarray(bwtRef[startRowIndex:endRowIndex+1]).flatten()).astype(int)
    return posi  
           
        
#--------------------------------------------------------------------------------------   
# String Matching Algorithm
# Takes the string input and outputs the index in the original chromosome
    
def PatternApproxMatch(inpString,maxMismatches=0):    
    
    global firstColStartIndex
    firstColStartIndex=(np.insert(np.cumsum(totCountEach),0,0)).astype(int)
    #print "firstColStartIndex=",firstColStartIndex
    
    # Initializations
    seqChar=inpString[-1]
    bandList=[]
    for char in ['A','C','G','T']:
        minRank=0
        maxRank=totCountEach[mapDict[char]]-1
        startBandIndex, endBandIndex = FindBand(char,minRank,maxRank)
     #   print "startBandIndex=",startBandIndex, "endBandIndex=", endBandIndex
        bandObj=bd.Band(char,startBandIndex,endBandIndex,int(char!=seqChar))  
        if bandObj.mismatches<=maxMismatches:
           bandList.append(bandObj) 
      #  print "bandList=",bandList
        
       
    # Find the matching bands
    for i in range(1,len(inpString)):
       
        seqChar=inpString[-1-i]
       # print i,"-th character",seqChar
        numOfBandObj=len(bandList)
        #print "numOfBandObj",numOfBandObj
        bandObjCount=0
        while (bandObjCount <numOfBandObj):
         #   print "inside second alphabet while loop"
            bandObj=bandList.pop(0)
            bandList.extend(FindSubBands(bandObj,seqChar,maxMismatches))
            bandObjCount=bandObjCount+1   
          #  print bandList
            #x=raw_input("hit enter when done.....................................................")

    
    #Find the reference indices
    #x=raw_input("done bandlist... check bandlist")
    #print "bandlist=",bandList
    matchPositions=[]
    confidenceLevel=[]
     
    for i in range(len(bandList)):  
          numMisMatches=(np.ones(bandList[i].endIndex-bandList[i].startIndex+1)*bandList[i].mismatches).astype(float)
          numMisMatches[numMisMatches>0]=0.5
          numMisMatches[numMisMatches==0]=1
          posi=(np.asarray(FindRefPosition(bandList[i].startIndex,bandList[i].endIndex)).flatten()).astype(int)
          matchPositions.extend(posi)
          confidenceLevel.extend(numMisMatches)
         # positions=(np.asarray(matchPositions)).astype(int)        
         
    readMatchPositions=dict(zip(matchPositions,confidenceLevel))
    #print readMatchPositions
    return readMatchPositions    


#--------------------------------------------------------------------------------------
# Associates each read with the Red and Green Exons depending on the 
# matched position



def BinningReads(readMatchPositions,currentRead) :
    for (position,confidence) in readMatchPositions.items():
        for (startBinIndex,endBinIndex) in exonsBinPositions.items():
            if position>=startBinIndex and position+len(currentRead) <=endBinIndex:
                exonBinReadsCount[startBinIndex]=exonBinReadsCount[startBinIndex]+confidence
                break
    

#--------------------------------------------------------------------------------------


ReadBWT(bwtLastColFileName,bwtRefFileName)

newkeys=list(exonsBinPositions.keys())
exonBinReadsCount={key:0 for key in newkeys}
counter=0
with open(readsFileName) as readsFileHand:
     for currentRead in readsFileHand:
         currentRead=currentRead.replace("N","")
         currentRead=currentRead.replace("\n","")
        # print "Current Read:",currentRead
         #print "Read count:",counter
        # currentRead="ATTCTACAAATACCCCTCCAGTATATTTCTTTCCCTTTTTTTTGAAGTCTCAGTCTGCCACCCAGGCTGGAGTGCAGTGGTGTGATTTTGGCTCACTGCA"
         readMatchPositions=PatternApproxMatch(currentRead,numOfAllowedMismatches)
         #print "position=",readMatchPositions
         BinningReads(readMatchPositions,currentRead)
         print exonBinReadsCount
         counter=counter+1
         
         
         
#to determine the probability of colour Blindness
totalCount=0
for i in exonBinReadsCount.values():
    totalCount=totalCount+i
   

#--------------------------------------------------------------------------
#Version 1 : Base version
#Version 2 : Changed fhand= open(readsFileName) to with open(readsFileName)
#           : to read the file line by line, that is one read at a time
#Version 3: changed ReadBWT function
#         : Added ACGTFetch function
#         : changed rank query  
#         : added accessing file fragments for bwtrefcol
#         : using bwtcol directly


    
            
            
    
    
    
          
        
          
      
  
    
    
    
    
    
    
    
    









