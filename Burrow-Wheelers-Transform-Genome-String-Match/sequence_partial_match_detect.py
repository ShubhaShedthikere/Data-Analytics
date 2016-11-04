import numpy as np
import Band as bd
#import Exons as ex


#-------------------------------------------------------------------------------------
#Input Parameters and other input information

exonsBinPositions={10:90, 120:185,229:232}
readsFileName="reads.txt"
chromosomeFileName=""
bwtLastColFileName="chrX_last_col.txt"
bwtRefFileName="chrX_map.txt"
numOfAllowedMismatches=1
delta=10


#-------------------------------------------------------------------------------------
# Global Variables

bwtLastCol=[]
ACGTBin=[]
ReadMatrixSave=[]
totCountEach=[]
firstColStartIndex=[]
bwtRef=[]
exonBinReadsCount={}
mapDict=dict({"A":0,"C":1,"G":2,"T":3})

#-------------------------------------------------------------------------------------
# read the BWT last column and create a matrix of cumulative sum
def ReadBWT(bwtColfilename,bwtRefFileName):
    global bwtLastCol
    global ACGTBin
    global ReadMatrixSave
    global totCountEach
    global bwtRef
    
    
    fhand=open(bwtRefFileName)
    bwtRefInp=fhand.read()
    bwtRef=bwtRefInp.split("\n")
    #print bwtRef
    
    fhand=open(bwtColfilename)
    bwtLastColInp=fhand.read()
    bwtLastColList=list(bwtLastColInp)
    bwtLastCol=np.array(bwtLastColList)
    bwtLastCol=bwtLastCol[ bwtLastCol!="\n"]
   
    ABinCol=(bwtLastCol=="A").astype(int)
    CBinCol=(bwtLastCol=="C").astype(int)
    GBinCol=(bwtLastCol=="G").astype(int)
    TBinCol=(bwtLastCol=="T").astype(int)
    ACGTBin = np.matrix([ABinCol,CBinCol,GBinCol,TBinCol]).T
    

    #cumulative sum
    ABinColCount=np.cumsum(ABinCol)
    CBinColCount=np.cumsum(CBinCol)
    GBinColCount=np.cumsum(GBinCol)
    TBinColCount=np.cumsum(TBinCol)
    ReadMatrix = np.matrix([ABinColCount,CBinColCount,GBinColCount,TBinColCount]).T

    #sample every delta-th value to store on disk
    delta = 10
    ReadMatrixSave = ReadMatrix[::delta]
   # print ReadMatrixSave
    totCountEach=np.asarray(ReadMatrix[-1,:]).flatten()
    #print "totCountEach",totCountEach

#--------------------------------------------------------------------------------------
#Rank query
def RankQuery(character,rowIndex):
    
    result=rowIndex/delta
    startIndex=result*delta
    rem=rowIndex % delta
    endIndex= startIndex+rem
    sampled=ACGTBin[startIndex+1:endIndex+1, mapDict[character]]
    return ReadMatrixSave[result,mapDict[character]]-1+np.sum(sampled)
    
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
    for char in ['A','C','G','T']:
        sampleBinCol=np.asarray(ACGTBin[bandObj.startIndex:bandObj.endIndex+1,mapDict[char]]).flatten()
        indexOfSeqChar=(np.array(np.where(sampleBinCol==1))).flatten()
        if (len(indexOfSeqChar)== 0):
            continue
        else:
            #for each if the entry in indexOfSeqChar obtain rank
           # print "indexOfSeqChar=",indexOfSeqChar
            minRank=RankQuery(char,indexOfSeqChar[0]+bandObj.startIndex)
            maxRank=RankQuery(char,indexOfSeqChar[-1]+bandObj.startIndex) 
            startBandIndex, endBandIndex = FindBand(char,minRank,maxRank)
            totBandMismatches = bandObj.mismatches + (char!=seqChar)
            if totBandMismatches<=maxMismatches:
                newBandObj=bd.Band(char,startBandIndex,endBandIndex,totBandMismatches) 
                bandReturnList.append(newBandObj)
            else:
                continue
    return bandReturnList

#--------------------------------------------------------------------------------------   
# String Matching Algorithm
# Takes the string input and outputs the index in the original chromosome
    
def PatternApproxMatch(inpString,maxMismatches=0):    
    
    global firstColStartIndex
    firstColStartIndex=np.insert(np.cumsum(totCountEach),0,0)
   # print "firstColStartIndex=",firstColStartIndex
    
    # Initializations
    seqChar=inpString[-1]
    bandList=[]
    for char in ['A','C','G','T']:
        minRank=0
        maxRank=totCountEach[mapDict[char]]-1
        startBandIndex, endBandIndex = FindBand(char,minRank,maxRank)
        bandObj=bd.Band(char,startBandIndex,endBandIndex,int(char!=seqChar))  
        bandList.append(bandObj)  
    
    # Find the matching bands
    for i in range(1,len(inpString)):
        seqChar=inpString[-1-i]
        numOfBandObj=len(bandList)
        bandObjCount=0
        while (bandObjCount <numOfBandObj):
            bandObj=bandList.pop(0)
            bandList.extend(FindSubBands(bandObj,seqChar,maxMismatches))
            bandObjCount=bandObjCount+1            
            
    #Find the reference indices
    matchPositions=[]
    confidenceLevel=[]
     
    for i in range(len(bandList)):  
          numMisMatches=(np.ones(bandList[i].endIndex-bandList[i].startIndex+1)*bandList[i].mismatches).astype(float)
          numMisMatches[numMisMatches>0]=0.5
          numMisMatches[numMisMatches==0]=1
          posi=(np.asarray(bwtRef[bandList[i].startIndex:bandList[i].endIndex+1]).flatten()).astype(int)
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
         readMatchPositions=PatternApproxMatch(currentRead,numOfAllowedMismatches)
         #print "position=",readMatchPositions
         
         BinningReads(readMatchPositions,currentRead)
         print "Read Count :",counter
         print "ExonBinReadsCount:",exonBinReadsCount
         print "------------------------------------------------------------"
         counter=counter+1


#--------------------------------------------------------------------------
#Version 1 : Base version
#Version 2 : Changed fhand= open(readsFileName) to with open(readsFileName)



    
            
            
    
    
    
          
        
          
      
  
    
    
    
    
    
    
    
    









