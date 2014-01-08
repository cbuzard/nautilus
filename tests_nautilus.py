#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that run a mercury simulation and test if the outputs and binaries have correct behaviour. 
# Of course, everything is not tested, but it is planed to test as many things as possible

__author__ = "Christophe Cossou <cossou@obs.u-bordeaux1.fr>"
__date__ = "21 Juillet 2011"
__version__ = "$Revision: 2.6.2 $"
__credits__ = """We run a test simulation and erase all the files created after the tests. The simulations files are thought to be 
in a "simu_test" subdirectory of the directory were are the sources (and binaries) of mercury (and this script)"""

import sys
import os
import difflib # To compare two strings
import subprocess # To launch various process, get outputs et errors, returnCode and so on.
import pdb # To debug
import glob # to get list of file through a given pattern

# Parameters
force = False # To force the compilation of every module

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * force : To force generation of outputs for the original program" + "\n" + \
" Example : " + "\n" + \
" tests_nautilus.py force"

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'force'):
    force = True
  elif (key == 'help'):
    isProblem = True
  else:
    print("the key '%s' does not match" % key)
    isProblem = True

if isProblem:
  print(problem_message)
  exit()


def run(commande):
  """lance une commande qui sera typiquement soit une liste, soit une 
  commande seule. La fonction renvoit un tuple avec la sortie, 
  l'erreur et le code de retour"""
  if (type(commande)==list):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  elif (type(commande)==str):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  else:
    raise TypeError("The command is neither a string nor a list.")
  (process_stdout, process_stderr) = process.communicate()
  returncode = process.poll()
  # there is .poll() or .wait() but I don't remember the difference. For some kind of things, one of the two was not working
  return (process_stdout, process_stderr, returncode)

def clean():
  """Delete all outputs files for a given nautilus simulation"""
  run("rm output_1D.*")
  run("rm rates1D.*")

def ASCIICompare(original, new):
  """function that compare and print differences between to strings that are compared line by line."""
  
  modified = ["+ ", "- ", "? "]
  
  # We want lists, because difflib.Differ only accept list of lines.
  original = original.split("\n")
  new = new.split("\n")
  
  d = difflib.Differ()

  result = list(d.compare(original, new))
  
  differences = []
  line_number_original = 0
  line_number_new = 0
  for line in result:
    if (line[0:2] == "  "):
      line_number_original += 1
      line_number_new += 1
    elif (line[0:2] == "- "):
      line_number_original += 1
      differences.append("[ori] l"+str(line_number_original)+" :"+line[2:])
    elif (line[0:2] == "+ "):
      line_number_new += 1
      differences.append("[new] l"+str(line_number_new)+" :"+line[2:])
    elif (line[0:2] == "? "):
      differences.append("      l"+str(max(line_number_new, line_number_original))+" :"+line[2:])

  # We print separately because it is more convenient if we want to store in a file instead.
  if (differences != []):
    return "\n".join(differences)
  else:
    return None

def compare2files(ori_files,new_files):
  """Function that will use compare to see differences between 'original' 
  that is thought to be a variable and 'new_file' that is the name of a 
  file to read then use as input
  """
  no_diff = []
  diff = []
  
  for (original, new) in zip(ori_files, new_files):
    f_old = open(original, 'r')
    old_lines = f_old.readlines()
    f_old.close()
    
    f_new = open(new, 'r')
    new_lines = f_new.readlines()
    f_new.close()
    
    
    difference = compare(''.join(old_lines), ''.join(new_lines))
    if (difference == None):
      no_diff.append(new)
    else:
      diff.append([new, difference])
  
  # Now we output results
  if (diff != []):
    for (file, comp) in diff:
      print("\nFor "+file)
      print(comp)
      
    if (no_diff != []):
      print "No differences seen on :",', '.join(no_diff)
  else:
    print("Everything OK")  
  
  return 0

def compare2Binaries(ori_files, new_files):
  """Function that will use compare to see differences between 'original' 
  that is thought to be a variable and 'new_file' that is the name of a 
  file to read then use as input
  """
  no_diff = []
  diff = []
  
  for (original, new) in zip(ori_files, new_files):
    (stdout, stderr, returnCode) = run("md5sum %s" % original)
    md5_ori = stdout.split()[0]
    
    (stdout, stderr, returnCode) = run("md5sum %s" % new)
    md5_new = stdout.split()[0]
    
    if (md5_new != md5_ori):
      diff.append(original)
    else:
      no_diff.append(original)
  
  # Now we output results
  if (diff != []):
    for filename in diff:
      print("\ndifferences with binary  %s" % filename)
      
    if (no_diff != []):
      print("No differences seen on :%s" % ', '.join(no_diff))
  else:
    print("Everything OK")
  
  return 0

ABUNDANCES_FILENAMES = glob.glob("output_1D.*")
RATES_FILENAMES = glob.glob("rates1D.*")

EXTENTION_ORIGINAL = ".ori"

ABUNDANCES_FILENAMES_OLD = ["%s%s" % (name, EXTENTION_ORIGINAL) for name in ABUNDANCES_FILENAMES]
RATES_FILENAMES_OLD = ["%s%s" % (name, EXTENTION_ORIGINAL) for name in RATES_FILENAMES]

NEW_TEST = "example_simulation"
ORIGINAL_TEST = "nautilus_original"

##################
# Outputs of various binaries and tests to compare with the actual ones. 
# Theses outputs are those of the original version of mercury, that is, mercury6_2.for
##################


os.chdir(NEW_TEST)
# We clean undesirable files. 
clean()

print("##########################################")
sys.stdout.write("Running new binaries ...\r")
sys.stdout.flush()

(naut_new_stdout, naut_new_stderr, returnCode) = run("../nautilus")

ABUNDANCES_FILENAMES = glob.glob("output_1D.*")
RATES_FILENAMES = glob.glob("rates1D.*")

os.chdir("..")
print("Running new binaries ...ok")
print("##########################################")
os.chdir(ORIGINAL_TEST)

if force:
  # We clean the output files
  clean()

if not(os.path.isfile("output_1D.000001")):
  sys.stdout.write("Running original binaries ...\r")
  sys.stdout.flush()
  (naut_or_stdout, naut_or_stderr, returnCode) = run("./nautilus")
  print("Running original binaries ...ok")
else:
  print("Skipping running original Nautilus, output already exists")
  # To prevent finding differences in the standard output and error
  naut_or_stdout = naut_new_stdout
  naut_or_stderr = naut_new_stderr

ABUNDANCES_FILENAMES_OLD = glob.glob("output_1D.*")
RATES_FILENAMES_OLD = glob.glob("rates1D.*")

os.chdir("..")

print("##########################################")

# We make the comparison

if (len(ABUNDANCES_FILENAMES_OLD) != len(ABUNDANCES_FILENAMES)):
  print("Error: number of abundances files is different")

if (len(RATES_FILENAMES_OLD) != len(RATES_FILENAMES)):
  print("Error: number of abundances files is different")

diff = ASCIICompare(naut_or_stdout, naut_new_stdout)
if (diff != None):
  print("\nTest of nautilus")
  print("\tFor the Output of nautilus")
  print diff

# We create names including the folder in which they are
RATES_FILENAMES = [os.path.join(NEW_TEST, rate) for rate in RATES_FILENAMES]
ABUNDANCES_FILENAMES = [os.path.join(NEW_TEST, abundance) for abundance in ABUNDANCES_FILENAMES]

# We create names including the folder in which they are
RATES_FILENAMES_OLD = [os.path.join(ORIGINAL_TEST, rate) for rate in RATES_FILENAMES_OLD]
ABUNDANCES_FILENAMES_OLD = [os.path.join(ORIGINAL_TEST, abundance) for abundance in ABUNDANCES_FILENAMES_OLD]

compare2Binaries(ABUNDANCES_FILENAMES_OLD, ABUNDANCES_FILENAMES)

compare2Binaries(RATES_FILENAMES_OLD, RATES_FILENAMES)
