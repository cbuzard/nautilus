#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that is a Makefile-like, but in python, and for the original mercury files

import subprocess
import os

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

FOLDER = "nautilus_original"

os.chdir(FOLDER)

command = "gfortran -O3  -ffast-math -c nls_header_mod.f90"
run(command)

command = "gfortran -O3  -ffast-math -o nautilus opk*.f nautilus.f90 nls*.f90 "
run(command)