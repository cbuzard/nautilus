#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import os
import sys
import pdb
import glob

isProblem = False
problem_message = "The aim is to update an old simulation folder so that the new binary will work.\n" + 
"The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" Example : " + "\n" + \
" updateSimulationFolder.py"

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'help'):
    isProblem = True
  else:
    print("the key '%s' does not match" % key)
    isProblem = True

if isProblem:
  print(problem_message)
  exit()

# We rename files when needed
ABUNDANCES_PREFIX = "output_1D."
ABUNDANCES_FILENAMES_OLD = glob.glob("%s*" % ABUNDANCES_PREFIX)
ABUNDANCES_FILENAMES_OLD.sort()
for filename in ABUNDANCES_FILENAMES_OLD:
  new_name = "abundances.%s.out" % filename.lstrip(ABUNDANCES_PREFIX)
  os.move(filename, new_name)

RATES_PREFIX = "rates1D."
RATES_FILENAMES_OLD = glob.glob("%s*" % RATES_PREFIX)
RATES_FILENAMES_OLD.sort()
for filename in RATES_FILENAMES_OLD:
  new_name = "rates.%s.out" % filename.lstrip(RATES_PREFIX)
  os.move(filename, new_name)
  
os.move('nls_init.d', 'chemical_composition.in')
os.move('nlso_tail.d', 'chemical_composition.tmp')
os.move('nlso_spec.d', 'species.out')
os.move('nls_control.d', 'parameters.in')

