#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that run a mercury simulation and test if the outputs and binaries have correct behaviour. 
# Of course, everything is not tested, but it is planed to test as many things as possible

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "21 Juillet 2011"
__version__ = "$Revision: 2.6.2 $"
__credits__ = """We run a test simulation and erase all the files created after the tests. The simulations files are thought to be 
in a "simu_test" subdirectory of the directory were are the sources (and binaries) of mercury (and this script)"""

files = ["opkda1.f90", "opkda2.f90", "opkdmain.f90"]

start_do = re.compile('do[ ]+[0-9]+', re.IGNORECASE) # If any, we want the label of the do-loop
end_do = re.compile('[0-9]+[ ]+continue', re.IGNORECASE) # If any, we want the label of the do-loop
for filename in files:
  
  objectFile = open(filename, 'r')
  
  lines = objectFile.readlines()
  
  nb_lines = len(lines)
  
  for i in range(nb_lines):
    
    line = lines[i]
    
    if (start_do.match(line) != None):
      words = line.split()
      labelID = int(words[1])
      
      # We create another index, because we will start from the actual line and search patterns in the rest of the code
      current_ID = i
      while (ongoing):
        current_ID += 1
        current_line = lines[current_ID]
        
        if (end_do.match(current_line)):
          TODO à finir
        
        
