#!/usr/bin/env python3.7

from dali import Dali
from pdb import Pdb
from os import path
import sys

tmpDir='.'
daliFile='1tenA.txt'
myRefPdb='1ten'
myRefChain='A'
myZ=2.
myTmpPath='./pdb'+myRefPdb+myRefChain
v=True

dali = Dali(refPdb=myRefPdb, refChain=myRefChain, tmpPath=myTmpPath, verbose=v)
#dali.CallDali(

dali.ReadFile(daliFile, verbose=v)
dali.RemoveDissimilar(z=myZ)
dali.GetSequenceFromPdb(tmpPath=myTmpPath, localPdb='', verbose=v, paranoid=False)
dali.RemoveInternalHomologs(byOrganism=True, dropNoOrganism=True, orderBy='lali')
dali.AlignToRef()
dali.PrintAlignment()
dali.CalculateEntropy(toPrint=True, filename='entropy_'+myRefPdb+myRefChain+'.dat')
dali.PrintAnalogs()
print('Number of sequences=',len(dali.analogs['aligned'].index))

