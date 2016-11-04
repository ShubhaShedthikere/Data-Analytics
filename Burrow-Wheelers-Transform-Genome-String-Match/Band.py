class Band:
     def __init__(self,char, startIndex, endIndex, mismatches):
         self.char=char
         self.startIndex=startIndex
         self.endIndex=endIndex
         self.mismatches=mismatches
     
     def __repr__(self):
         return str(self)
     
     def __str__(self):
        return("Char: "+str(self.char)+"\n"+"startIndex:"+str(self.startIndex)+\
                "\n"+"endIndex:"+str(self.endIndex)+"\n"+"Mismatches:"\
                  +str(self.mismatches)) 
         
