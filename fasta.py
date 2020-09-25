import sys
import pandas as pd
import textwrap

class Fasta:

  def __init__(self):

    self.df = pd.DataFrame(columns = [ 'seq', 'id' ])
    self.nseq = 0
    self.length = 0
  
  def _ReplaceChar(self, string, i, char):
    return string[:i-1] + char + string[i:] 

  def ReadFile(self, filename, verbose=True):

    try: 
       fp = open(filename,"r")      
    except:
      sys.exit("Cannot read alignment "+filename)

    title = '' 
    tmpSeq = ''
    for i, line in enumerate(fp):

      line = line.strip()

      if line[0:1] == '>':

        if title !='':
          self.df = self.df.append({'seq':tmpSeq, 'id':title}, ignore_index=True)
          if self.nseq > 1 and len(tmpSeq) != self.length:
            sys.exit("Fasta alinment contains sequences of different lengths")
          self.length = len(tmpSeq)
          self.nseq = self.nseq + 1
          tmpSeq=''
          title=''
   
        title = line[1:]

      elif title !='':
        tmpSeq = tmpSeq + line
        

    if title!='':
      self.df = self.df.append({'seq':tmpSeq, 'id':title}, ignore_index=True)

    if verbose:
     print('Read '+str(self.nseq)+' sequences of length '+str(self.length))



  def RemoveGapColumns(self):
    
    mask = []
    for i in range(self.length): mask.append(0)

    for i, x in self.df.iterrows():
      for j in range(self.length):
       if x['seq'][j]!='-' and x['seq'][j]!='X': mask[j] = 1
 
    
    for i, x in self.df.iterrows():
      for j in range(self.length):
        if mask[j]==0: x['seq'] = self._ReplaceChar(x['seq'],j,'*')
      x['seq'] = x['seq'].replace('*','')
   
    self.length = 


  def Print(self, filename=''):
   
    if filename != "":
      try: 
         fp = open(filename,"w")
      except:
         sys.exit("ERROR: cannot open file to write pdb")
    else:
      fp = sys.stdout
 
    
    for i, x in self.df.iterrows():
      print('>'+x['id'], file=fp)
      print(textwrap.fill(x['seq'],80), file=fp) 
   
    if (fp != sys.stdout):
      fp.close() 
