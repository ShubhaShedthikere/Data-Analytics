import numpy as np

delta=10
bwtLastCol=[]
ACGTBin=[]
ReadMatrixSave=[]
totCountEach=[]
bwtRef=[]

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
    print bwtRef
    
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
    print ReadMatrixSave
    totCountEach=np.asarray(ReadMatrix[-1,:]).flatten()
    print "totCountEach",totCountEach






#Rank query
def RankQuery(character,rowIndex):
    mapDict=dict({"A":0,"C":1,"G":2,"T":3})
    result=rowIndex/delta
    startIndex=result*delta
    rem=rowIndex % delta
    endIndex= startIndex+rem
    sampled=ACGTBin[startIndex+1:endIndex+1, mapDict[character]]
    return ReadMatrixSave[result,mapDict[character]]-1+np.sum(sampled)
    
    
#Returns the row Indices of the band of the corresponding characters in 
#first column    
def FindBand(character, minRank,maxRank,firstColStartIndex):
    mapDict=dict({"A":0,"C":1,"G":2,"T":3})
    startBandIndex=firstColStartIndex[mapDict[character]]+minRank
    return startBandIndex, startBandIndex+(maxRank-minRank)
    
    
# String Matching Algorithm
# Takes the string input and outputs the index in the original chromosome
    
def PatternMatch(inpString,bwtLastColFileName,bwtRefFileName):    
    ReadBWT(bwtLastColFileName,bwtRefFileName)
    firstColStartIndex=np.insert(np.cumsum(totCountEach),0,0)
    print "firstColStartIndex=",firstColStartIndex
    isFirstIteration=True
    mapDict=dict({"A":0,"C":1,"G":2,"T":3})
    prevChar=None
    for i in range(len(inpString)):
        print "------------------------------------------------------------"
        seqChar=inpString[-1-i]
        if (isFirstIteration):
            minRank=0
            maxRank=totCountEach[mapDict[seqChar]]-1
            isFirstIteration=False
            prevChar=seqChar
            print "minRank=",minRank
            print "maxRank=",maxRank
            continue
        else:
            #Find Band
            startBandIndex, endBandIndex = FindBand(prevChar,minRank,maxRank,firstColStartIndex)
            print "startBandIndex=",startBandIndex
            print "endBandIndex=",endBandIndex
            sampledBwtLastCol=bwtLastCol[startBandIndex:endBandIndex+1]
            indexOfSeqChar=(np.array(np.where(sampledBwtLastCol==seqChar))).flatten()
        
        # if the sequence is not found then break,else continue matching
        if (indexOfSeqChar== ""):
            break
        else:
            #for each if the entry in indexOfSeqChar obtain rank
            print indexOfSeqChar
            minRank=RankQuery(seqChar,indexOfSeqChar[0]+startBandIndex)
            maxRank=RankQuery(seqChar,indexOfSeqChar[-1]+startBandIndex)
            
        prevChar=seqChar
        print "minRank=",minRank
        print "maxRank=",maxRank
            
    if (indexOfSeqChar==" "):
          return "not found"
    else:
           startBandIndex, endBandIndex = FindBand(prevChar,minRank,maxRank,firstColStartIndex)
           matchPositions=bwtRef[startBandIndex:endBandIndex+1]
           positions=(np.asarray(matchPositions)).astype(int)
    
    return positions       
          
        
          
      
  
    
    
    
    
    
    
    
    









