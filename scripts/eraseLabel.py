#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "21 Juillet 2011"
__version__ = "$Revision: 2.6.2 $"
__credits__ = """Script that suppress all labels for loops in a F77 style source code. 
This only works if no do while loops, or exotics loops exists"""

import shutil
import re
import pdb

files = ["opkda1.f90", "opkda2.f90", "opkdmain.f90"]
#~ 
#~ files = ["test.f90"]
#~ 
#~ shutil.copy2("opkda1.f90", "test.f90")

def changeDo(line):
  """create a new string that remove the f77 label of the loop"""
  tmp = line.split()
  del(tmp[1])
  # When splitting, the \n disappear
  tmp.append("\n")

  
  return " ".join(tmp)

def changeContinue(line):
  """create a new string that modify the end of an old f77 do-loop : Version 1, with continue"""
  
  return "enddo\n"
  
def changeEnddo2(line):
  """create a new string that modify the end of an old f77 do-loop : Version 2, the line past the label is inside the loop"""
  tmp = line.split()
  del(tmp[0])
  
  tmp = " ".join(tmp)
  tmp += "\nenddo\n"
  
  return tmp

start_do = re.compile('^do[ ]+[0-9]+', re.IGNORECASE) # If any, we want the label of the do-loop
#~ end_do = re.compile('^[0-9]+[ ]+continue', re.IGNORECASE) # If any, we want the label of the do-loop
end_do = re.compile('^[0-9]+', re.IGNORECASE) # If any, we want the label of the do-loop
for filename in files:
  
  objectFile = open(filename, 'r')
  
  lines = objectFile.readlines()
  objectFile.close()
  
  nb_lines = len(lines)
  
  #~ pdb.set_trace()
  
  for i in range(nb_lines):
    
    line = lines[i]
    
    if (start_do.match(line) != None):
      words = line.split()
      labelID = int(words[1])
      
      #~ if labelID==140:
        #~ pdb.set_trace()
      
      # We create another index, because we will start from the actual line and search patterns in the rest of the code
      current_idx = i
      ongoing = True
      while (ongoing):
        # The next step will be the last
        if (current_idx>=nb_lines-2):
          ongoing = False
        current_idx += 1
        
        current_line = lines[current_idx]
        
        if (end_do.match(current_line)):
          words = current_line.split()
          
          #~ if (current_idx==143):
            #~ pdb.set_trace()
          
          # if this is the right ID
          tmp_label = int(words[0])
          if (tmp_label == labelID):
            ongoing = False
            lines[i] = changeDo(lines[i])
            
            if (words[1].lower() == "continue"):
              tmp_line = changeContinue(lines[current_idx])
            else:
              tmp_line = changeEnddo2(lines[current_idx])
            lines[current_idx] = tmp_line
          
  
  # We erase the old file with the new one
  objectFile = open(filename, 'w')
  
  for line in lines:
    objectFile.write(line)
  
  objectFile.close()
