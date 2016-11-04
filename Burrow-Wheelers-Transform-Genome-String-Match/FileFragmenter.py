
#This function opens the file and reads N lines of the file and puts it into a
# new file

import numpy as np
lineCount=0
chr_Map_fragment=[]
chrName="chrMapFragmentPart"

with open("chrX_map.txt") as infile:
     for line in infile:
        
         if lineCount%100000==0:
            partNum=lineCount/100000
            fileName=chrName+ str(partNum)
            if partNum!=0:
              fhand.close()
            fhand=open(fileName,'w')
         
         print>>fhand, line.replace("\n","")
         lineCount=lineCount+1
         
         #if lineCount>=500000:
          #   break
            
         
 
print np.asarray(chr_Map_fragment )
print "done"       
        
        