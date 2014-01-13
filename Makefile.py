#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that is a Makefile-like, but in python, and for the original mercury files

import subprocess
import os
import sys
import pdb

COMPILATOR = "gfortran"
DEBUG_OPTIONS = "-pedantic-errors -Wall -Wconversion -Wunderflow -Wextra -Wunreachable-code -fbacktrace" + \
  " -ffpe-trap=invalid,zero,overflow,underflow -g3 -fbounds-check -O0" + \
  " -fstack-protector-all -fno-automatic -Wuninitialized -ftrapv -fno-automatic"
OPTIMIZATIONS = "-O0 -march=native -ffast-math -pipe -finit-real=nan"
GDB_OPTIONS = "-g3"
PROFILING_OPTIONS = "-pg"

LOG_NAME = "compilation.log"

# Parameters
debug = False
gdb = False
profiling = False
force = False # To force the compilation of every module

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * force : To force the compilation of every module even those not modified" + "\n" + \
" * debug : [%s] activate debug options" % debug + "\n" + \
" * gdb : [%s] activate options for gdb" % gdb + "\n" + \
" * profiling : [%s] activate options for profiling" % profiling + "\n" + \
" Example : " + "\n" + \
" Makefile.py gdb"

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
  if (key == 'debug'):
    debug = True
  elif (key == 'force'):
    force = True
  elif (key == 'gdb'):
    gdb = True
  elif (key == 'profiling'):
    profiling = True
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
  returnCode = process.poll()
  
  # If returnCode is not 0, then there was a problem
  if (returnCode != 0):
    # We write compilation errors in the following file.
    f = open(LOG_NAME,'w')
    f.write(process_stderr)
    f.close()
    
    print("Compilation error, see '%s'" % LOG_NAME)
    sys.exit(1)
  else:
    if (len(process_stderr) != 0):
      # We write compilation errors in the following file.
      f = open(LOG_NAME,'a')
      f.write(process_stderr)
      f.close()
    
    print("Warnings: see '%s'" % LOG_NAME)
  
  # there is .poll() or .wait() but I don't remember the difference. For some kind of things, one of the two was not working
  return (process_stdout, process_stderr, returnCode)

if debug:
  OPTIONS = DEBUG_OPTIONS
else:
  OPTIONS = OPTIMIZATIONS

if gdb:
  OPTIONS = GDB_OPTIONS

if profiling:
  OPTIONS = PROFILING_OPTIONS

# Before compiling, we delete the previous compilation log. Indeed, we need to append the several warnings in the same file
# But we do not want to have infos of the previous compilation in it.
if os.path.isfile(LOG_NAME):
  os.remove(LOG_NAME)

command = "%s %s -c nls_header_mod.f90 ode_solve.f90" % (COMPILATOR, OPTIONS)
print(command)
(process_stdout, process_stderr, returncode) = run(command)

#~ pdb.set_trace()

#~ command = "%s %s -o nautilus opk*.f90 nautilus.f90 nls*.f90" % (COMPILATOR, OPTIONS)
command = "%s %s -o nautilus ode_solve.o nautilus.f90 nls*.f90" % (COMPILATOR, OPTIONS)
print(command)
(process_stdout, process_stderr, returncode) = run(command)

#~ pdb.set_trace()
