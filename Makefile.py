#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that is a Makefile-like, but in python, and for the original mercury files

import subprocess
import os
import sys
import pdb

COMPILATOR = "gfortran"
DEBUG_OPTIONS = "-pedantic-errors -Wall -Wconversion -Wunderflow -Wextra -Wunreachable-code -fbacktrace" + \
  " -ffpe-trap=invalid,zero,overflow -g3 -fbounds-check -O0" + \
  " -fstack-protector-all -fno-automatic -Wuninitialized -ftrapv -fno-automatic"
OPTIMIZATIONS = "-O3 -ffast-math -pipe -finit-real=nan"
GDB_OPTIONS = "-g3"
PROFILING_OPTIONS = "-g -pg"

LOG_NAME = "compilation.log"

nautilus_order = ["-c numerical_types.f90 iso_fortran_env.f90 utilities.f90 git_infos.f90", "-c global_variables.f90", 
"-c shielding.f90 diffusion.f90 input_output.f90 model_1D.f90 ode_solver.f90", 
"-o nautilus nautilus.f90 opk*.f90 *.o"]

output_order = ["-c numerical_types.f90 iso_fortran_env.f90 utilities.f90", "-o nautilus_outputs nautilus_outputs.f90 *.o"]

# Parameters
debug = False
gdb = False
profiling = False
force = False # To force the compilation of every module
compilation_order = nautilus_order

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * force : To force the compilation of every module even those not modified" + "\n" + \
" * output : To compile binary abundances instead of Nautilus" + "\n" + \
" * debug : [%s] activate debug options" % debug + "\n" + \
" * gdb : [%s] activate options for gdb" % gdb + "\n" + \
" * profiling : [%s] activate options for profiling" % profiling + "\n" + \
" Example : " + "\n" + \
" Makefile.py gdb"

def run_command(commande):
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

def get_current_branch():
  """function that return as a string the current branch of the git repository"""
  (stdout, stderr, returnCode) = run_command("git branch")
  
  if (returnCode != 0):
    return None
  
  lines = stdout.split("\n")
  for line in lines:
    if (line[0] == '*'):
      return line[2:]

def get_current_revision():
  """function that return as a string the current revision of the git repository"""
  (stdout, stderr, returnCode) = run_command("git log|head -1")
  
  if (returnCode != 0):
    return None
  
  commit = stdout.split()[1]
  return commit

def is_non_committed_modifs():
  """function that return as a boolean if theere is non committed modifications in the repository"""
  (stdout, stderr, returnCode) = run_command("git diff|wc -l")
  
  if (returnCode != 0):
    return None
  
  nbLines = int(stdout)
  
  return (nbLines != 0)
  
def list_tag(commit):
  """list the tags that exists linking towar the considered commit
  
  Return :
  The list of tags corresponding to 'commit'. If none, an empty list is returned.
  """
  (stdout, stderr, returnCode) = run_command("git tag -l --contains %s" % commit)
  
  tags = stdout.split("\n")[0:-1] # We do not include the extra "" in the end.
  
  return tags 

def write_infos_in_f90_file(main_branch='master'):
  """This function will create a fortran file that will store, as variable, some infos about a git repository"""
  
  F90_BEGIN = """!******************************************************************************
! MODULE: git_infos
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Automatically generated module (with Makefile.py) that contain
!! information on code version and branch name. 
!
!> @warning Do not modify this file manually !
!
!******************************************************************************
module git_infos
implicit none

"""

  F90_END = "\nend module git_infos"

  
  branch = get_current_branch()
  commit = get_current_revision()
  isModifs = is_non_committed_modifs()
  tags = list_tag(commit)
  
  if (branch != main_branch):
    print("Warning: The current branch is %s" % branch)
  
  f90source = open("git_infos.f90", 'w')
  f90source.write(F90_BEGIN)
  f90source.write("character(len=40), parameter :: commit = '%s' !< commit ID when binary was compiled \n" % commit)
  f90source.write("character(len=%d), parameter :: branch = '%s' !< branch name when compilation occured\n" % (len(branch), branch))
  if (isModifs):
    f90source.write("character(len=80), parameter :: modifs = '/!\ There is non committed modifications'\n")
  else:
    f90source.write("character(len=80), parameter :: modifs = 'This is a pure version (without any local modifs)'\n")
    
  
  if (len(tags)==0):
    tag_text = "There is no tag"
  else:
    tag_text = " ; ".join(tags)
  
  f90source.write("character(len=%d), parameter :: tags = '%s'\n" % (len(tag_text), tag_text))
  f90source.write(F90_END)
  f90source.close()
  
def run_compilation(commande):
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
  
  # if returnCode is not 0, then there was a problem
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
  elif (key == 'output'):
    compilation_order = output_order
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

write_infos_in_f90_file(main_branch='master')

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



for order in compilation_order:
  command = "%s %s %s" % (COMPILATOR, OPTIONS, order)
  print(command)
  (process_stdout, process_stderr, returncode) = run_compilation(command)

#~ pdb.set_trace()
