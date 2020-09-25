import sys,os
sys.path.append('~/prog/pyTools')
from pdb import Pdb



pdbName = sys.argv[1]
path, fname = os.path.split(pdbName)
pdbCode = fname.replace('.pdb','')

pdb = Pdb()
pdb.verbose = True
pdb.ReadFile(pdbName)

for c in pdb.chains:
  chainType = pdb.IdentifyChainType(c)
  if chainType == 'protein': 
    newPdb = pdb.Copy(chain=c)
    
    nFileName = pdbCode + c + '.pdb'

    newPdb.Print(filename=nFileName)

