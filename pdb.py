# Class Pdb	
# v. 0.1
# G. Tiana, 18 May 2020
#
# atom dict = {'ia', 'atom', 'residue', 'chain', 'iResidue', 'x', 'y', 'z', 'record', 'model'} 

import sys
from atomStatistics import AtomStatistics
import textwrap
from urllib import request
from os import system

class Pdb:

   aminoAcids = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL' }
   dnaAcids = { 'A', 'C', 'G', 'T' }
   rnaAcids = { 'A', 'C', 'G', 'U' }
   ions = {'NA', 'CL', 'MN', 'MG'}

   lowerMassProtein = 3000


   def __init__(self):
   
      self.atom = []
      self.pdbName = ''
      self.nAtoms = 0
      self.chains = ''
      self.nResidues = 0
      self.nChains = 0
      self.nModels = 1
      self.verbose = False
      self.atomStatistics = AtomStatistics()
      self.header = []


   ##################################################################
   def ReadFile(self, filename, onlyFirstModel=True, skipWater=True, skipIons=True, skipNonAminoacid=False, keepChainLabel=True, pdbName='', verbose=False):
    
     iResidueOld = -1
     if pdbName=='':
       self.pdbName = filename.replace('.pdb','')
     else: self.pdbName = pdbName
     previousRecord = ""
     previousChain=""
     previousChainTmp=""
     chainTmp=""
     chain = '@'

     if verbose: 
       print("Reading file "+filename)
                            
     # open file
     try:
       fp = open(filename,"r")
     except: 
       sys.exit("ERROR: file "+filename+" not found.")

     # loop on file lines
     for i, line in enumerate(fp):

       #get keyword of pdb line
       record = line[0:6].strip()

       #get header
       if self.nAtoms == 0 and record!="ATOM":
          self.header.append(line.strip())

       #read the pdb line
       if record=="ATOM" or record=="HETATM":
         try:
          ia = int(line[6:12])
          myAtom = line[12:16].strip()
          alternate = line[16].strip()
          residue = line[17:20].strip()          
          chainTmp = line[21]
          iResidue = int(line[22:26])
          x = float(line[30:38]) 
          y = float(line[38:46]) 
          z = float(line[46:54]) 

          # discard water and ions
          if skipWater and residue == 'HOH': continue
          if skipIons and residue in self.ions and myAtom in self.ions: continue
          if skipNonAminoacid and ( residue not in self.aminoAcids ): continue

          # discard alternates other than A
          if alternate != '' and alternate !='A': continue

          # handle chain label
          if keepChainLabel:
             chain = chainTmp
          else:
             if previousRecord == 'TER' or chainTmp != previousChainTmp  :
             #if previousRecord == 'TER' or chainTmp != previousChainTmp or (record=='ATOM' and previousRecord=='HETATM') or (record=='HETATM' and previousRecord=='ATOM')  :
                chain = chr(ord(chain)+1)
            

          #fill the atom structure
          self.nAtoms = self.nAtoms + 1
          self.atom.append( {'ia':self.nAtoms, 'atom':myAtom, 'residue':residue, 'chain':chain, 'iResidue':iResidue, 
                              'x':x, 'y':y, 'z':z, 'record':record, 'model':self.nModels} )

          #count number of chains
          if self.chains.find(chain) == -1 :
            self.chains = self.chains + chain
            self.nChains = len(self.chains)   

          #count number of residue
          if iResidue != iResidueOld:
            self.nResidues = self.nResidues +1
            iResidueOld = iResidue       

         except:
           if self.verbose:
             print("WARNING: error in reading line "+str(i+1)+" of pdb. Skipped.")


       previousRecord = record
       previousChainTmp = chainTmp

       if record == "ENDMDL":
         self.nModels = self.nModels + 1 
         if onlyFirstModel:
           break
 
     if self.verbose:
        print(" Read "+str(self.nAtoms)+" atoms.")      
        print(" There are "+str(self.nChains)+" chains ("+self.chains+")")      

     if self.nModels > 1:
         print(" WARNING: There are "+str(self.nModels)+" models in the pdb.")

     fp.close()

     if self.nAtoms == 0:
       sys.exit("ERROR: Could not read any atom in the pdb.")



   ##################################
   def MolecularMass(self, chain=""):

      m = 0

      for a in self.atom:
        if chain=="" or a['chain']==chain:
          
            if a['atom'][0] == 'H': m = m + 1
            if a['atom'][0] == 'C': m = m + 13
            if a['atom'][0] == 'O': m = m + 16
            if a['atom'][0] == 'N': m = m + 14
            if a['atom'][0] == 'P': m = m + 31
            if a['atom'][0] == 'S': m = m + 32

      return m


   ######################################
   # returns: { protein, peptide, dna, rna, chemical }
   def IdentifyChainType(self, chain="", verbose=True):
  
      nProt= 0
      nDna = 0
      nRna = 0
      nOther = 0
      chainType="unknown"

      if verbose:
        print("Identify chain type:")
 
      if chain not in self.chains:
        sys.exit("ERROR: chain "+chain+" not present in pdb "+self.pdbName)

      # molecular mass of system or chain
      mass = self.MolecularMass(chain)

      # count residues types
      for a in self.atom:
       if chain=="" or chain == a['chain']:
         if a['residue'] in self.aminoAcids: nProt = nProt + 1
         elif a['residue'] in self.dnaAcids: nDna = nDna + 1
         elif a['residue'] in self.rnaAcids: nRna = nRna + 1
         else: nOther = nOther + 1

      if verbose:
         if chain != "": print(" Chain "+chain+":", end='')
         print( "  mass="+str(mass)+" nProt="+str(nProt)+" nDna="+str(nDna)+" nRna="+str(nRna)+" nOther="+str(nOther) )

      # selection rules
      if nProt>0 and nDna==0 and nRna==0 and nOther==0 and mass > self.lowerMassProtein: chainType="protein"
      elif nProt>0 and nDna==0 and nRna==0 and nOther==0 and mass < self.lowerMassProtein: chainType="pepetide"
      elif nProt==0 and nDna>0 and nDna>nRna and nOther==0 : chainType="dna"
      elif nProt==0 and nDna>0 and nDna<nRna and nOther==0 : chainType="rna"
      elif nOther>0 and nDna==0 and nRna==0 and nProt==0 and mass < self.lowerMassProtein: chainType="chemical"
      else: chainType="unknown"

      if verbose:
         print("Chain type is "+chainType)

      return chainType


   #####################################
   def ProteinChains(self, verbose=True):
                       
     chain = ''

     for c in self.chains:
       ct = self.IdentifyChainType(c, verbose=False)
       if ct == 'protein':
         chain = chain + c

     if chain == '' and verbose:
       print('WARNING: no protein chain found')
     
     return chain

   #####################################
   def Print(self, filename=""):
    

     if filename != "":
       try: 
         fp = open(filename,"w")
       except:
         sys.exit("ERROR: cannot open file to write pdb")
     else:
       fp = sys.stdout

     if self.pdbName != "":
       print("TITLE "+self.pdbName, file=fp)
 
     for a in self.atom:
        print("%-6s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f" % ("ATOM", a['ia'], a['atom'], "", a['residue'], a['chain'], a['iResidue'], "",a['x'], a['y'], a['z']) ,file=fp)

     print("END",file=fp)

     if filename != sys.stdout:
       fp.close()


   
   #####################################
   def Copy(self, chain=""):
     
      if chain not in self.chains: 
        sys.exit("Chain "+chain+" not in pdb to be copied")
 
      newPdb = Pdb()

      newPdb.pdbName = self.pdbName + chain
      newPdb.nAtoms = 0
      iResidueOld = -1
      newPdb.nResidues = 0
      newPdb.verbose = self.verbose

      if chain == "": 
        newPdb.chains = self.chains
        newPdb.nChains = self.nChains
        newPdb.nModels = self.nModels
      else: 
        newPdb.chains = chain
        newPdb.nChains = 1
        newPdb.nModels = 1

      print(len(self.atom))
      for a in self.atom:
        if chain == "" or a['chain'] == chain:
          newPdb.nAtoms = newPdb.nAtoms + 1 
          newPdb.atom.append( {'ia':newPdb.nAtoms, 'atom':a['atom'], 'residue':a['residue'], 'chain':a['chain'], 
                             'iResidue':a['iResidue'], 'x':a['x'], 'y':a['y'], 'z':a['z'], 
                              'record':a['record'], 'model':a['model']} )

          if a['iResidue'] != iResidueOld:
            newPdb.nResidues = newPdb.nResidues +1
            iResidueOld =  a['iResidue']       




      return newPdb

   #####################################
   def __AtomDistanceSquared(self, a, b):
    
      return (a['x']-b['x'])**2 + (a['y']-b['y'])**2 + (a['z']-b['z'])**2  


   #####################################
   def SelectChainByDistance(self, chainFromList, chainRefList, distance):


      if self.verbose:
        print(" Select atoms with distance in chain(s) "+chainFromList+" closer than "+str(distance)+"A from chain(s) "+chainRefList)


      newPdb = Pdb()
      iResidueOld = -1
      newPdb.nAtoms = 0
      newPdb.pdbName = self.pdbName + chainFromList + 'selected'
      newPdb.verbose = self.verbose
      newPdb.nResidues = 0

      for b in self.atom:
        if b['chain'] in chainFromList:
          for a in self.atom:
             if a['chain'] in chainRefList:
               d = self.__AtomDistanceSquared(a, b)
               if d < distance**2 :
                    newPdb.atom.append( {'ia':b['ia'], 'atom':b['atom'], 'residue':b['residue'], 'chain':b['chain'], 
                             'iResidue':b['iResidue'], 'x':b['x'], 'y':b['y'], 'z':b['z'], 
                              'record':b['record'], 'model':b['model']} )
                    newPdb.nAtoms = newPdb.nAtoms + 1 
                    if a['iResidue'] != iResidueOld:
                      newPdb.nResidues = newPdb.nResidues +1
                      iResidueOld =  a['iResidue']       

      if self.verbose:
         print(" selected "+str( newPdb.nAtoms)+" atoms")

      return newPdb


   #####################################
   def GetResidues(self):
   
     residues = []
     iRes = -1

     for a in self.atom:
       if a['residue'] != iRes:			#only one per residue
         residues.append(a['residue'])
         iRes = a['residue']

     return residues

   #####################################
   def CalculateAtomStatistics(self, verbose=False):

      for a in self.atom:
        atomType = a['atom'][0]
        try:      
            self.atomStatistics.histo[atomType] = self.atomStatistics.histo[atomType] + 1
        except:
            if verbose:
              print("WARNING: atom type "+a['atom']+" not recognized")

   
   #####################################
   def GetChainWithResidue(self, residue, verbose=True):

      chain = ''      

      for a in self.atom:
        if a['residue'] == residue:
          if chain == '':
            chain = a['chain']      
          elif chain != a['chain'] and verbose:
            print('WARNING: more chains contain residue '+residue+', namely '+chain+' and '+a['chain'])
            print('Keeping only '+chain)
            break

      if chain == '' and verbose:
        print('WARNING: no residue '+residue+' found')

      return chain 

   #####################################
   def PrintFasta(self, chain='', filename='', fileOpt='w', cut=True):

      if filename != "":
       try: 
         fp = open(filename,fileOpt)
       except:
         sys.exit("ERROR: cannot open file to write pdb")
      else:
       fp = sys.stdout

      print(">"+self.pdbName, file=fp)
      seq = self.GetSequence(chain)
      if cut == False:
        print(seq, file=fp)
      else:
        print(textwrap.fill(seq,80), file=fp)

      if filename != "": fp.close()

      
   #####################################
   def GetSequence(self, chain="", verbose=False):
 
     seq=''
     oldIResidue = -1

     for a in self.atom:
        if ( chain=='' or a['chain']==chain ):
          if ( a['iResidue'] != oldIResidue ):
            oldIResidue = a['iResidue']  
            aa = a['residue']

            c='X'
            if ( aa == 'ALA'):  c='A'
            if ( aa == 'ARG'):  c='R'
            if ( aa == 'ASN'):  c='N'
            if ( aa == 'ASP'):  c='D'
            if ( aa == 'ASX'):  c='B'
            if ( aa == 'CYS'):  c='C'
            if ( aa == 'GLU'):  c='E'
            if ( aa == 'GLN'):  c='Q'
            if ( aa == 'GLX'):  c='Z'
            if ( aa == 'GLY'):  c='G'
            if ( aa == 'HIS'):  c='H'
            if ( aa == 'ILE'):  c='I'
            if ( aa == 'LEU'):  c='L'
            if ( aa == 'LYS'):  c='K'
            if ( aa == 'MET'):  c='M'
            if ( aa == 'PHE'):  c='F'
            if ( aa == 'PRO'):  c='P'
            if ( aa == 'SER'):  c='S'
            if ( aa == 'THR'):  c='T'
            if ( aa == 'TRP'):  c='W'
            if ( aa == 'TYR'):  c='Y'
            if ( aa == 'VAL'):  c='V'

            seq = seq+c

            if verbose:
               print( str(a['iResidue'])+aa+' --> '+c )

     return seq

  
   #####################################
   def Download(self, pdbCode, path='.', verbose=False):

      url = 'http://files.rcsb.org/download/'+pdbCode+'.pdb'

      if verbose: print(url)
      

      try:  
        request.urlretrieve(url, path+'/'+pdbCode+'.pdb')
      except:
        sys.exit('Canot download pdb '+pdbCode)
      
      self.ReadFile(filename=path+'/'+pdbCode+'.pdb')


   #####################################
   def GetKey(self, key, verbose=False):

      for line in self.header:
         words = line.split()
         for i,w in enumerate(words):
               if w.replace(':','') == key:
                 if len(words)<i+2: 
                    if verbose: print("No value found for key "+key)
                    return ''
                 foundKey = ' '.join(words[i+1:]).replace(';','').replace(',','')
                 if verbose: print("Key "+key+' is '+foundKey)
                 return foundKey

      if verbose: print("No key "+key+" found")

      return ''

   #####################################
   def TakeLocalFiles(self, pdbCode, sourceDir, targetDir):

      sys.exit('IT DOES NOT WORK')
      
      subDir = pdbCode[1:3]
      sourceFile = sourceDir+'/'+subDir+'/'+'pdb'+pdbCode+'.ent.gz'
      targetFile = targetDir+'/'+pdbCode+'.pdb'
                                    
    #  try:
      cmd = 'cp '+sourceFile+' '+targetDir+'/'
      system(cmd)
      cmd = 'gunzip '+targetDir+'/'+'pdb'+pdbCode+'.ent.gz' 
      system(cmd)
      cmd = 'mv '+ targetDir+'/'+'pdb'+pdbCode+'.ent.gz' +' '+targetFile
      system(cmd)
     # except:
  #      sys.exit('Cannot get'+targetFile+' from local directory')
 
   def PrintSummary(self):

      print('Summary of pdb '+self.pdbName)
      print(self.nAtoms,'atoms')
      print(self.nChains,'chains',self.chains)
      print(self.nModels,'models')


