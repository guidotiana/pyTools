from pdb import Pdb
import sys, os
from Bio.Align.Applications import MuscleCommandline

verbose = True

if len(sys.argv) != 2:
  sys.exit('ERROR: use analyzePdbBunch pdblist.txt')

os.remove('tmp.fasta')

# open pdbList
try:
   fp = open(sys.argv[1],"r")
except:
   sys.exit("ERROR: file "+pdbList+" not found.")

# loop on pdbList
for i,line in enumerate(fp):

  # read pdbCode and chain
  pdbFile = line.strip()
  path, pdbFile = os.path.split(pdbFile)
  if pdbFile[0:3]=='pdb': pdbCode=pdbFile[3:]
  else: pdbCodeC=pdbFile
  pdbCodeC = pdbCodeC.replace('.pdb','')
  pdbCodeC = pdbCodeC.replace('.ent','')
  if len(pdbCodeC) == 4:
    pdbCode = pdbCodeC[0:4]
    chain = pdbCode[4]
  else:
    pdbCode = pdbCodeC
    chain = ''

  if verbose: print('Reading '+pdbCode+chain)


  # read pdb
  pdb = Pdb()
  pdb.ReadFile(filename=pdbFile,  keepChainLabel=True)
  pdb.PrintFasta(filename='tmp.fasta', fileOpt='a', chain=chain)
  del pdb

# sequence alignment
muscle_cline = MuscleCommandline('/Users/guido/prog/muscle-3.8.31', input='tmp.fasta', out='align.fasta')
print(muscle_cline)
muscle_cline()

