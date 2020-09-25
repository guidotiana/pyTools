import pandas as pd
from pdb import Pdb
from os import path,mkdir
import sys
from math import log
from os import system

class Dali:

   def __init__(self, refPdb, refChain, tmpPath='.', localPdb='', verbose=False):

      self.analogs = pd.DataFrame(columns = [ 'pdb', 'chain', 'z', 'rmsd', 'lali', 'nres', 'id', 'title', 'sequence', 'aligned', 'organism_id', 'organism_scientific' ]) 
      self.equivalences = pd.DataFrame(columns = [ 'pdb2', 'chain2', 'from1', 'to1', 'from2', 'to2' ])
      self.refPdb = refPdb
      self.refChain = refChain
      self.refSequence = ''
      self.refLength = 0
      self.analogs['organism_id']=''

      if verbose: print('Creating Dali object for pdb '+refPdb+refChain)

      pdb=Pdb()

      if not path.isdir(tmpPath):
        try:
         mkdir(tmpPath)
        except:
         sys.exit('Cannot make dir '+tmpPath)

      if not path.isfile(tmpPath+'/'+refPdb+'.pdb'):
         if localPdb=='':
           pdb.Download(pdbCode=refPdb, verbose=verbose, path=tmpPath)
         else:
           pdb.TakeLocalFiles(pdbCode=refPdb, sourceDir=localPdb, targetDir=tmpPath )
      
      pdb.ReadFile(tmpPath+'/'+refPdb+'.pdb', skipNonAminoacid=True)
      self.refSequence = pdb.GetSequence(chain=refChain)
      del pdb
      self.refLength = len(self.refSequence)

      if verbose: print('Length of reference sequence is '+str(self.refLength))


   def __OpenOutput(self, filename, objName=''):

      if filename != "":
       try: 
         fp = open(filename,"w")
       except:
         sys.exit("ERROR: cannot open file to write "+objName)
      else:
       fp = sys.stdout

      return fp


   def ReadFile(self, filename, verbose=False):

     # open file
     try:
       fp = open(filename,"r")
     except: 
       sys.exit("ERROR: file "+filename+" not found.")

     # loop on file lines
     readAnalogs = False
     readEquivalences = False

     if verbose: print("Reading file "+filename)

     for i, line in enumerate(fp):
        
        line2 = line.replace('-',' ')
        word = line2.split()
                    
        if len(word)==3:
           if word[1]=='Query:' and self.refPdb=='': 
              self.refPdb = word[2][0:4]
              self.refChain = word[2][4]
              if verbose: print('Reference pdb is '+word[2])
            
        # read analogs
        if readAnalogs == True and len(word)>7 :
          self.analogs = self.analogs.append( { 'pdb': word[1],
                                                'chain': word[2],
                                                'z': word[3],
                                                'rmsd': word[4],
                                                'lali': word[5],
                                                'nres': word[6],
                                                'id': word[7],
                                                'title': ' '.join(word[8:]),
                                              }, ignore_index=True )
        #read equivalences
        if readEquivalences == True and len(word)>10 :
           self.equivalences = self.equivalences.append( {  'pdb2': word[3],
                                                            'chain2':word[4],
                                                            'from1': word[5],
                                                            'to1': word[6],
                                                            'from2': word[8],
                                                            'to2': word[9],
                                                         }, ignore_index=True )



        # decide what to read
        if len(word)>1:
           if word[0] == '#' and word[1]=='No:': 
              readAnalogs = True
           if word[0] == '#' and word[1]=='Structural': 
              readAnalogs = False
              readEquivalences = True
           if word[0] == '#' and word[1]=='Translation-rotation': 
              readAnalogs = False
              readEquivalences = False

     self.analogs['z'] = self.analogs['z'].astype(float)
     self.analogs['rmsd'] = self.analogs['rmsd'].astype(float)
     self.analogs['lali'] = self.analogs['lali'].astype(int)
     self.analogs['nres'] = self.analogs['nres'].astype(int)
     self.analogs['id'] = self.analogs['id'].astype(int)
     self.equivalences['from1'] = self.equivalences['from1'].astype(int)
     self.equivalences['from2'] = self.equivalences['from2'].astype(int)
     self.equivalences['to1'] = self.equivalences['to1'].astype(int)
     self.equivalences['to2'] = self.equivalences['to2'].astype(int)

     if verbose: print('Read '+str(len(self.analogs.index))+' proteins in Dali file')

   def RemoveDissimilar(self, z):

      self.analogs.drop( self.analogs[ self.analogs['z'] < z ].index, inplace=True )


   def RemoveInternalHomologs(self, id=-1, byOrganism=False, dropNoOrganism=False, orderBy=''):

      if id>-1:
        self.analogs.drop( self.analogs[ self.analogs['id'] > id ].index, inplace=True )      

      if byOrganism:
        if self.analogs["organism_id"].empty == True: sys.exit("Cannot remove internal homologs before reading pdb")
        if orderBy !='':
           self.analogs.sort_values(by=orderBy, ascending=0)
        self.analogs.drop_duplicates( subset="organism_id", keep=False, inplace=True)
        if dropNoOrganism: self.analogs.drop(  self.analogs[ self.analogs['organism_id'] =='' ].index, inplace=True  )


   def GetSequenceFromPdb(self, tmpPath='.', localPdb='', verbose=False):
 
      for i,prot in self.analogs.iterrows():

        pdb2name = prot['pdb']
        myChain = prot['chain']

        if verbose: print(" reading file "+tmpPath+'/'+pdb2name+'.pdb')

        pdb=Pdb()
        if not path.isfile(tmpPath+'/'+pdb2name+'.pdb'):
         if localPdb=='':
           if verbose: print(" downloading file")
           pdb.Download(pdbCode=pdb2name, verbose=True, path=tmpPath)
           pdb.ReadFile(tmpPath+'/'+pdb2name+'.pdb', skipNonAminoacid=True)
         else:
           if verbose: print(" copying file from local directory "+localPdb)
           pdb.TakeLocalFiles(pdbCode=pdb2name, sourceDir=localPdb, targetDir=tmpPath )
           pdb.ReadFile(tmpPath+'/'+pdb2name+'.pdb', skipNonAminoacid=True)
        else:
            pdb.ReadFile(tmpPath+'/'+pdb2name+'.pdb', skipNonAminoacid=True)
    
        seq = pdb.GetSequence(chain=myChain)

        self.analogs.at[i, 'sequence'] = seq
        
        self.analogs.at[i, 'organism_id'] = pdb.GetKey('ORGANISM_TAXID')
        self.analogs.at[i, 'organism_scientific'] = pdb.GetKey('ORGANISM_SCIENTIFIC')
                                 
        if len(seq)==0: sys.exit("ERROR: Empty sequence for pdb "+pdb2name+myChain)

        if verbose:
          print(' has read pdb '+pdb2name+myChain+' (length='+str(len(seq))+')')
          print(seq)

        del pdb

   def __SubstituteChar(self, string, pos, char):
                               
        if pos<0: return string
        if pos>len(string): return string 
        return string[0:pos]+char+string[pos+1:]

   def AlignToRef(self):

      self.analogs['aligned'] = '-' * self.refLength

      for i,prot in self.analogs.iterrows():
    
        pdb2name = prot['pdb']
        myChain = prot['chain']
        tmpSeq = prot['aligned']     
        pdb2Sequence = prot['sequence']

        for k,equiv in self.equivalences.iterrows():
           if equiv['pdb2']==prot['pdb'] and equiv['chain2']==prot['chain']:
              deltaj = equiv['from2']  - equiv['from1'] 

              for j in range( equiv['from1'], equiv['to1']+1 ):
                if j+deltaj>=0 and j+deltaj-1<self.refLength:
                 tmpSeq = self.__SubstituteChar( tmpSeq, j-1, pdb2Sequence[j+deltaj-1] )

        self.analogs.at[i, 'aligned'] = tmpSeq

   def PrintAlignment(self, filename=''):

      fp = self.__OpenOutput(filename,'alignment')

      print(self.refPdb + self.refChain + '   ' + self.refSequence, file=fp)
      for i,prot in self.analogs.iterrows():
         print(prot['pdb'] + prot['chain'] + '   ' + prot['aligned'], file=fp)

   


   def CalculateEntropy(self, toPrint=False, filename=''):

       
      fp = self.__OpenOutput(filename,'entropy')

      for iSite in range(self.refLength):
   
         prob = {'A':0, 'B':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'J':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0,  
                 'S':0, 'T':0, 'U':0, 'V':0, 'W':0, 'X':0, 'Y':0, 'Z':0, '-':0}

         tot=0
         for i, prot in self.analogs.iterrows():
            aa = prot['aligned'][iSite]
            prob[aa] = prob[aa] + 1
            tot = tot + 1         

         s=0
         for x in prob:
           if prob[x]>0:
              s = s - prob[x]/tot * log( prob[x]/tot )

         if toPrint == True:
             print(iSite+1,s, file=fp)

 
   def PrintAnalogs(self):

      print("%5s %5s %3s %7s %s"  % ('pdb',  'z', 'id%', 'organism', 'scientific' ))
      for i,x in self.analogs.iterrows():
        print( "%5s %5.1f %3d %7s %s" % (x['pdb']+x['chain'], x['z'], x['id'], x['organism_id'], x['organism_scientific']) )

   def __ReducedAlphabet(self, c):

      if c=='A' or c=='V' or c=='L' or c=='I' or c=='M' or c=='C': return 'A'
      if c=='F' or c=='W' or c=='Y' or c=='H': return 'R'
      if c=='S' or c=='T' or c=='N' or c=='Q': return 'P'
      if c=='D' or c=='E': return 'N'
      if c=='K' or c=='R': return 'S'
      if c=='G': return 'G'
      if c=='P': return 'P'
      if c=='-': return '-'
      if c=='X': return 'X'
      return '*'

   def ReduceAlphabet(self, verbose=False):

      if verbose: print("Alphabet reduced to 7 symbols")

      for i,x in self.analogs.iterrows():
        tmpSeq = x['aligned']
        for j in range(0, self.refLength ):
                 tmpSeq = self.__SubstituteChar( tmpSeq, j, self.__ReducedAlphabet(tmpSeq[j]) )

        self.analogs.at[i, 'aligned'] = tmpSeq


   def CallDali(self, daliPath, refPdb, verbose=False):

      if verbose: print("Calling Dali tool")


      cmd = daliPath + '/import.pl --pdbfile '+pdbFile+' --pdbid '+ refPdb+refChain + ' --dat' + daliPath + '/dat --clean'
      print(cmd)
