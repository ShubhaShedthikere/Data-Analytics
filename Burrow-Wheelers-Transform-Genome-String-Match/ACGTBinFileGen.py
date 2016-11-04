import numpy as np
 # form the binary version of the BWT column
ABinCol=[]
CBinCol=[]
GBinCol=[]
TBinCol=[]
bwtLastColFileName="chrX_last_col.txt"
linecount=0
n=20
print "Forming A C G T binary version of BWT last column.............."
with open(bwtLastColFileName) as fhand:
     for line in fhand:
              text=line
              textA= text.replace('A', '1').replace('C', '0').replace('G', '0').replace('T', '0').replace('\n','').replace('$','')         
              textA=[int(textA[i:i+n],2) for i in range(0, len(textA), n)]
              ABinCol.extend(textA)
              textC= text.replace('A', '0').replace('C', '1').replace('G', '0').replace('T', '0').replace('\n','').replace('$','')         
              textC=[int(textC[i:i+n],2) for i in range(0, len(textC), n)]
              CBinCol.extend(textC)
              textG= text.replace('A', '0').replace('C', '0').replace('G', '1').replace('T', '0').replace('\n','').replace('$','')         
              textG=[int(textG[i:i+n],2) for i in range(0, len(textG), n)]
              GBinCol.extend(textG)
              textT= text.replace('A', '0').replace('C', '0').replace('G', '0').replace('T', '1').replace('\n','').replace('$','')         
              textT=[int(textT[i:i+n],2) for i in range(0, len(textT), n)]
              TBinCol.extend(textT)
              linecount=linecount+1
              
     print "Num of elements=",len(ABinCol)
     ACGTBin = np.matrix([ABinCol,CBinCol,GBinCol,TBinCol]).T
     np.savetxt("ACGTBinMatrix.txt",ACGTBin)

