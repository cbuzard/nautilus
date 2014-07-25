#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that is a Makefile-like, but in python, and for fortran90

import subprocess
import git_infos
import os
import sys
import glob
import string

LOG_NAME = 'compilation.log'
SIMULATION_FOLDER = "tests"
GNUPLOT_EXTENSION = "gnuplot"

VISUALISEUR = {"svg":"gthumb", "png":"gthumb", "pdf":"gv", "jpg":"gthumb"}

# Parameters
debug = False
gdb = False
profiling = False
force = False # To force the compilation of every module
indice_script = None # to allow default actions if nothing is given in parameters
generate_all = False

isProblem = False
problem_message = "The script can take various arguments :" + "\n" + \
"(no spaces between the key and the values, only separated by '=')" + "\n" + \
" * help : display a little help message on HOW to use various options" + "\n" + \
" * force : To force the compilation of every module even those not modified" + "\n" + \
" * debug : [%s] activate debug options" % debug + "\n" + \
" * script=1 to avoid be asked what test we want. (index refer to a test, 'all' if we want to execute them all)" + "\n" + \
" * gdb : [%s] activate options for gdb" % gdb + "\n" + \
" * profiling : [%s] activate options for profiling" % profiling + "\n" + \
" Example : " + "\n" + \
" unitary_tests.py gdb"

class sourceFile(object):
  """Define an object linked to a fortran 90 source code that will
  store dependencies, including modules included, used or included
  
  Parameters : 
  filename : the name of an existing source code filename (for example "main.f90")
  name='default' : the name we want for a program (if it is a program, if not, there will be no use)
  isProgram=False : is this source file a main source for a binary we want, or just a module or sub-program?
  
  Attribute : 
  self.filename : the filename of the source code file.
  self.defined : a list of procedures defined in the fortran source code
  self.used : a list of modules that are used by the code
  self.included : a list of things included in the code
  self.isCompiled : a boolean to say if the source file has already been compiled.
  self.dependencies : a list of object (*.o) filenames we need to compile
  self.isProgram : a boolean to say if we want to have a binary, or if it's just a module or a subprogram
  self.name : the name we want for the binary file if it is a program. This name is by default the filename without the extension
  
  Methods :
  .compile() : Compile the current program and all the required dependencies
  .writeArchitecture(filename,excluded=[], direction="leftright") : To write a .svg file showing the tree of dependencies of the current sourceFile.
    filename : The filename of the .svg file, without the extension
    excluded : a list of names of the modules you don't want to print in the diagram
    direction="leftright" : How the diagram will be constructed
      leftright -> The main program will be on the left, and the dependencies will appear on the right, row by row for each 
        'generation' of children dependencies
      topbottom -> The main program will be on the top, and the dependencies will appear on the bottom, line by line for each 
        'generation' of children dependencies

  """
  
  ## We define a dictionnary where we store, for each name of a module,
  # the link toward the object file where it is defined
  findModule = {}
  
  # We define a dictionary to make the correspondance between a source filename and the object associated
  findSource = {}
  #~ COMPILATOR = "ifort"
  #~ OPTIONS = "-vec-report0 -i-dynamic -mcmodel=medium -shared-intel -L/usr/lib64/atlas -llapack"
  
  COMPILATOR = "gfortran"
  OPTIONS = "-O3 -march=native"
  
  # Even with theses compilation options, it misses some warning. I had a weird behavior where some warnings only showed when I 
  # had more serious issues.
  DEBUG = "-pedantic-errors -Wall -Wconversion -Wunderflow -Wextra -Wunreachable-code -fbacktrace" + \
  " -g3 -fbounds-check -O0" + \
  " -fstack-protector-all -fno-automatic -Wuninitialized -ftrapv -fno-automatic"
  GDB = "-g3"
  
  
  
  #-Wextra : batterie supplémentaire de vérifications
  #-ffast-math : je l'ai enlevé car les résultats ne sont pas identiques, les derniers chiffres significatifs sont différents.
  #-march=native : permet d'optimiser pour le processeur de la machine sur laquelle on compile. 'native' permet d'aller chercher 
  #   cette information sur la machine au lieu de la spécifier à la main. Si l'option ne fonctionne pas, typiquement si 'native' ou 
  #   le type de processeur spécifié n'existe pas, une erreur est retournée.
  #-fimplicit-none : empêche les déclarations implicites à moins que le mot clé "implicit" ne soit explicitement utilisé.
  #-finit-real=zero : initialise tous les réels à 0
  # farfadet spatial m'a conseillé -O2 au lieu de -O3 mais je ne comprends pas encore pourquoi.

  # Boolean that say if we want to activate debug or not
  isDebug = False
  isGDB = False
  isProfiling = False
  
  def __init__(self, filename, name='default', isProgram=False):
    """Will check everything that is included in the source code
    and initialize the object"""
    self.filename = filename
    
    if (name == 'default'):
      self.name = self.filename.rstrip(".f90")
    else:
      self.name = str(name)
    
    
    
    # By default, nothing is a program
    self.isProgram = isProgram
    
    
    self.isCompiled = False
    
    if (not(self.isProgram)):
      # If the source file is newer than the object file, we need to compile it
      object_file = "%s.o" % self.name
      if (os.path.isfile(object_file) and os.path.isfile("%s.mod" % self.name)):
        self.toBeCompiled = os.path.getmtime(self.filename) > os.path.getmtime(object_file)
      else:
        # If the object file do not exist
        self.toBeCompiled = True
    else:
      self.toBeCompiled = True # We always compile the programs
    
    (self.defined, self.used, self.included) = self.__getModules()
    
    for module in self.defined:
      sourceFile.findModule[module] = self
    
    sourceFile.findSource[self.filename] = self
  
  @classmethod
  def setDebug(cls, isDebug):
    """method that define cls.isDebug parameter to True or False.
    
    Parameter: isDebug (boolean)
    """
    
    cls.isDebug = isDebug

  @classmethod
  def setGDB(cls, isGDB):
    """method that define cls.isGDB parameter to True or False.
    
    Parameter: isGDB (boolean)
    """
    
    cls.isGDB = isGDB
  
  @classmethod
  def setProfiling(cls, isProfiling):
    """method that define cls.isProfiling parameter to True or False.
    
    Parameter: isProfiling (boolean)
    """
    
    cls.isProfiling = isProfiling
  
  @classmethod
  def setModColors(cls):
    """class function that gives a dictionnary of colors for each sourceFile defined so far
    """
    try:
      cls.mod_colors = dict(zip(cls.findSource.keys(), autiwa.colorList(len(cls.findSource))))
    except:
      pdb.set_trace()
    return 0
  
  @classmethod
  def setCompilingOptions(cls, options):
    """method that set the 'OPTIONS' value. 
    
    Parameter:
    options : a string with the compilation options for the binary construction
    """
    
    cls.OPTIONS = options
  
  @classmethod
  def setCompilator(cls, compilator):
    """method that set the 'COMPILATOR' value. 
    
    Parameter:
    options : a string with the compilator you want to use
    """
    
    cls.COMPILATOR = compilator
  
  def __getModules(self):
    """returns a tuple containing the list of defined modules and 
    the list of used modules of a fortran source file
    
    Return 
    defined : a list of procedures defined in the fortran source code
    used : a list of modules that are used by the code
    included : a list of things included in the code
    """

    f=open(self.filename,'r')
    lines = f.readlines()
    f.close()

    defined=[]
    used=[]
    included=[]

    for lsave in lines:
      l=string.expandtabs(string.lower(lsave)[:-1],1)
      words=string.split(string.lstrip(l))
      if len(words) > 0:
        if words[0] == 'use':
          used.append(string.split(words[1],',')[0])
        if words[0] == 'module':
          if len (words) == 2 or words[1] != "procedure":
            defined.append(words[1])
        if words[0] == 'include':
          newstring = string.replace(words[1],'\'','')
          newstring = string.replace(newstring,'\"','')
          included.append(newstring)
      l=string.expandtabs(lsave[:-1],1)
      words=string.split(string.lstrip(l))
      if len(words) > 0:
        if words[0] == '#include':
          newstring = string.replace(words[1],'\"','')
          included.append(newstring)

# We delete all dependencies that are present several number of times.
    used = list(set(used))

    return defined,used,included
  
  def __getFirstOrderDependence(self):
    """return a list of *.o files we need to compiled. They are 
    extracted from direct used files, that why we call them 
    "first order" dependance."""
    
    dependances = []
    for mod in self.used:
      try:
        source = sourceFile.findModule[mod]
      except:
        print("Error: Unable to locate the module '"+mod+"'")
      obj = string.replace(source.filename,'.f90','.o')
      dependances.append(obj)
    
    return dependances
  
  def setProgram(self, boolean):
    """method to defined the current object as a program, that is, if
    we want to have a binary from this source"""
    
    self.isProgram = boolean
  
  def writeArchitecture(self,filename,excluded=[], direction="leftright"):
    """write the architecture in a .svg file"""
    
    if (direction == "topbottom"):
      constructArch = self.__getArchitectureTopBottom
    elif (direction == "leftright"):
        constructArch = self.__getArchitectureLeftRight
    else:
        raise ValueError("the 'direction' you typed do not exist. You must choose between: \n_ 'topbottom'\n_ 'leftright'")
    
    # We define the dictionary of colors for each source file defined.
    self.setModColors()
    
    architecture = []
    (architecture, x, y, total_width) = constructArch(architecture=architecture,excluded=excluded)
    
    createSVG(filename, *architecture)

  def __getArchitectureTopBottom(self, architecture, excluded=[], x=0, y=0):
    """Retrieve the architecture of the program
    
    Parameter:
    x : The x position of the first box
    y : The y position of the parent box"""
    
    width = 150
    height = 30
    
    width_space = 10
    height_space = 50
    
    # we search for the dependencies of the current source file.
    dependencies = []
    used = list(self.used) # list() is here to avoid virtual link

    # We do not print excluded modules
    for mod in excluded:
      try:
        used.remove(mod)
      except:
        pass
    
    #~ for source in sourceFile.findSource.keys():
      # TODO là je mettais la ligne pour générer le dictionnaire des couleurs, mais il faut que ce soit un paramètre de classe
      # sinon il n'est pas connu des autres instances

    # For each modules, we search the defined object in the base of sourcefiles
    for mod in used:
      dependencies.append(sourceFile.findModule[mod])

    nb_dep = len(dependencies)
    
    # We search for children and print lines to join them to the present
    # box if there are any dependencies (except thoses excluded)
    if (nb_dep != 0):
      
      tree_width = 0
      x_cursor = x
      y_cursor = y + height + height_space# * len(dependencies) # may cause problems for close boxes that might overlap
      child_coord = []
      for dep in dependencies:
        (architecture, x_cursor, y_cursor, childchild_width) = dep.__getArchitectureTopBottom(architecture, excluded=excluded, x=x_cursor, y=y_cursor)
        child_coord.append((x_cursor, y_cursor))
        tree_width += childchild_width + width_space
        x_cursor += childchild_width/2. + width/2. + width_space
      
      # For the last element, we don't want the space between boxes
      tree_width = tree_width - width_space
      parent_y = y

      # We calculate the coords for the parent box
      parent_x = x + tree_width / 2. - width /2.
      
      # We write lines between the parent and each child
      for (child_x, child_y) in child_coord:
        temp = Path((parent_x, parent_y+height/2.), (child_x, child_y-height/2.))
        architecture.append(temp)
      
    else:
      tree_width = width
      
      parent_x = x
      parent_y = y
    
        
    architecture.append(TextBox(self.name, parent_x-width/2., parent_y-height/2., color=sourceFile.mod_colors[self.filename]))

    return (architecture, parent_x, parent_y, tree_width)

  def __getArchitectureLeftRight(self, architecture, excluded=[], x=0, y=0):
    """Retrieve the architecture of the program
    
    Parameter:
    x : The x position of the first box
    y : The y position of the parent box"""
    
    width = 150
    height = 30
    
    width_space = 50
    height_space = 10
    
    # we search for the dependencies of the current source file.
    dependencies = []
    used = list(self.used) # list() is here to avoid virtual link

    # We do not print excluded modules
    for mod in excluded:
      try:
        used.remove(mod)
      except:
        pass

    # For each modules, we search the defined object in the base of sourcefiles
    for mod in used:
      dependencies.append(sourceFile.findModule[mod])

    nb_dep = len(dependencies)
    
    # We search for children and print lines to join them to the present
    # box if there are any dependencies (except thoses excluded)
    if (nb_dep != 0):
      
      tree_width = 0
      x_cursor = x + width + width_space
      y_cursor = y
      child_coord = []
      for dep in dependencies:
        (architecture, x_cursor, y_cursor, childchild_width) = dep.__getArchitectureLeftRight(architecture, excluded=excluded, x=x_cursor, y=y_cursor)
        child_coord.append((x_cursor, y_cursor))
        tree_width += childchild_width + height_space
        y_cursor += childchild_width/2. + height/2. + height_space
      
      # For the last element, we don't want the space between boxes
      tree_width = tree_width - height_space
      parent_x = x

      # We calculate the coords for the parent box
      parent_y = y + tree_width / 2. - height /2.
      
      # We write lines between the parent and each child
      for (child_x, child_y) in child_coord:
        temp = Path((parent_x + width/2., parent_y), (child_x-width/2., child_y))
        architecture.append(temp)
      
    else:
      tree_width = height
      
      parent_x = x
      parent_y = y
    
        
    architecture.append(TextBox(self.name, parent_x-width/2., parent_y-height/2., color=sourceFile.mod_colors[self.filename]))

    return (architecture, parent_x, parent_y, tree_width)

  def compile(self, parent_dependencies=[]):
    """method that check dependencies and try to compile
    
    Parameter : 
    parent_dependencies=[] : list that store all the parent dependencies of the current file, namely, all the module 
    that MUST NOT be used inside the current module, else we will have an infinite loop.
    
    Beware, there MUST NOT be inter dependant modules.
    """
    
    parent_dependencies.append(self.name)
    
    if not(self.isCompiled):
      # We only store, for the moment, the first order dependencies.
      self.dependencies = self.__getFirstOrderDependence()
      
      # We store links towards all the object sourceFile that defined the modules we are interested in.
      module_sources = []
      for module in self.used:
        module_sources.append(sourceFile.findModule[module])
        
      # For each object, we check if there is loop call of modules. If that's the case, return an error
      for name in self.used:
        if name in parent_dependencies:
          error_message = "The module '"+self.name+"' try to use the module '"+name+"' that already "+ \
                          "use the module '"+self.name+"'. So there is an infinite loop that is not correct."
          raise NameError(error_message)
      
      # For each object, we compile it if it's not already the case.
      for source in module_sources:
        # We must launch source.compile() even if the file will not be compilated, just to get the right dependencies.
        if not(source.isCompiled):
          # the list() is here to ensure not to have a pointer and share the list. If not, the list of parent_dependencies will 
          # not be correct and contains all the previous parent dependencies. 
          source.compile(list(parent_dependencies)) 
          
          # We compile the current file if any of the dependencies has been effectively compiled
          if source.toBeCompiled:
            self.toBeCompiled = True
          
      # We complete the dependencies list now that all used modules
      # have been compiled, they must have a complete list of 
      # their own dependencies.
      for source in module_sources:
        self.dependencies.extend(source.dependencies)
        
      # We delete all dependencies that are present several number of times.
      self.dependencies = list(set(self.dependencies))
      
      if self.toBeCompiled:
        # Now that all the dependencies have been compiled, we compile 
        # the current source file.
        options = sourceFile.OPTIONS
        if (sourceFile.isDebug):
          options += " "+sourceFile.DEBUG
        
        if (sourceFile.isGDB or sourceFile.isProfiling):
          options += " "+sourceFile.GDB
        
        if (sourceFile.isProfiling):
          # We deactivate all other options except GDB => not True anymore
          options += " "+" -pg"
          
        if not(self.isProgram):
          commande = sourceFile.COMPILATOR+" "+options+" -c "+self.filename
        else:
          commande = sourceFile.COMPILATOR+" "+options+" -o "+self.name+" "+self.filename+" "+" ".join(self.dependencies)
          print(commande)
        
        process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

        (process_stdout, process_stderr) = process.communicate()
        
        print("Compiling "+self.filename+"...")
        returnCode = process.poll()
      
        self.isCompiled = True
      
        # If returnCode is not 0, then there was a problem
        if (returnCode==0):
          if (len(process_stderr) != 0):
            LOG_NAME = "compilation.log"
            # We write compilation errors in the following file.
            f = open(LOG_NAME,'a')
            f.write(process_stderr)
            f.close()
      
            print("Warnings: see '%s'" % LOG_NAME)
          
          return process_stderr
        else:        
          logname = "compiling_"+self.filename+".log"

          # We write compilation errors in the following file.
          f = open(logname,'w')
          f.write(process_stderr)
          f.close()
          
          print("Compilation error, see '"+logname+"'")
          sys.exit(1)
        
      self.isCompiled = True
      
  def __str__(self):
    """overload the str method. As a consequence, you can print the object via print name_instance

    return : a string that represent the properties of the object
    """

    texte = "filename: "+self.filename+"\n"
    
    texte += "defined: "+str(self.defined)+"\n"
    texte += "used: "+str(self.used)+"\n"
    texte += "included: "+str(self.included)+"\n"

    return texte

def make_binaries(sources_filename, mains, debug=False, gdb=False, profiling=False):
  """function that will compile every needed sourceFile and get a 
  binary for the defined source files
  
  Parameters : 
  sources_filename : a list of source filenames that must be present in the current working directory for which we will define an object 'sourceFile'
  mains : either a list of filename or a dictionnary for which keys are filename and values are the name of the binary file we want.
   i.e we define a list of filename if we want that the binary has the same name (without extension) as the source file, 
   or a dictionnary to make the correspondance between the two.
  debug=False : (boolean) Whether we want or not debug options for compilation. There is no debug by default
  gdb=False : If set to True, will add a compilation option to run the program under the GNU debugger gdb. 
  profiling=False : If set to True, will change compilation options to profile the binary files (using gprof). 
                    As a consequence, this will deactivate all optimization options. 
  
  Examples : 
  make_binaries(sources_filename, {"mercury6_2.for":"mercury", "element6.for":"element", "close6.for":"close"})
  or
  make_binaries(sources_filename, ["mercury.f90", "element.f90", "close.f90"])
  where 'sources_filename' is a list of all the sources file with the extension '*.for' and '*.f90' respectively.
  
  """
  
  # if the mains are defined as a list, that means the default name will be the default one, controlled by the class. 
  if (type(mains) == list):
    main_list = list(mains)
    mains = {}
    for main in main_list:
      mains[main] = 'default'
  elif (type(mains) == dict):
    main_list = mains.keys()
  else:
    raise TypeError("'mains' must be a dict or a list")
  
  # We define the objects for each source file.
  sources = []
  main_source = []
  for filename in sources_filename:
    if (filename in main_list):
      source = sourceFile(filename, name=mains[filename], isProgram=True)
    else:
      source = sourceFile(filename)
    sources.append(source)
  
  sourceFile.setDebug(debug)
  sourceFile.setGDB(gdb)
  sourceFile.setProfiling(profiling)
  
  # We compile the programs (dependencies are automatically compiled if needed.
  for source in sources:
    if (not(source.isCompiled) and source.isProgram):
      print source.compile()

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
  
def enumeration(liste):
  """liste, un élément par ligne, les éléments de l'argument
  avec les indices associés"""
  output_list = []
  for (indice,element) in enumerate(liste):
    txt = "%2i : %-34s" % (indice, os.path.splitext(element)[0])
    output_list.append(txt)
  
  print(" ".join(output_list))

def getOutput(scriptname):
  """function that return the name of the output file of the given gnuplot script, if there is any"""
  
  f=open(scriptname,'r')
  lines = f.readlines()
  f.close()
  
  output = None
  
  for line in reversed(lines):
    words = line.split()
    if (len(words) > 1):
      if (words[0] == "set") and (words[1] == "output"):
        output = words[2][1:-1]
        break
        
  # in case of output file printed not in the directory of the .dat file
  temp = output.split("/")
  output = temp[-1]
  
  return output
  
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
  elif (key == 'script'):
    if (value == 'all'):
      generate_all = True
      indice_script = 0
    else:
      generate_all = False
      indice_script = int(value)
  elif (key == 'profiling'):
    profiling = True
  elif (key == 'help'):
    isProblem = True
  else:
    print("the key '"+key+"' does not match")
    isProblem = True



if isProblem:
  print(problem_message)
  exit()

git_infos.write_infos_in_f90_file(main_branch='master')

isModifs = git_infos.is_non_committed_modifs()

# We clean undesirable files. Indeed, we will compile everything everytime.
if force:
  clean(["o", "mod"])

sourceFile.setCompilator("gfortran")
sourceFile.setCompilingOptions("-O3 -march=native -pipe -finit-real=nan")
#~ sourceFile.setCompilingOptions("")

# pour tester les bornes des tableaux : -fbounds-check (il faut ensuite faire tourner le programme, des tests sont effectués au cours de l'exécution)

sources_filename = glob.glob("*.f90")

# Before compiling, we delete the previous compilation log. Indeed, we need to append the several warnings in the same file
# But we do not want to have infos of the previous compilation in it.
if os.path.isfile(LOG_NAME):
  os.remove(LOG_NAME)

# We create the binaries
make_binaries(sources_filename, ["unitary_tests.f90"], debug=debug, gdb=gdb, profiling=profiling)

if (isModifs):
  print("Warning: There is non committed modifs!")
  

# We run the unitary tests
os.chdir(SIMULATION_FOLDER)

# We delete the old gnuplot scripts
run("rm *.gnuplot")

(process_stdout, process_stderr, returncode) = run("../unitary_tests")
print(process_stdout)

if (returncode != 0):
  print(process_stderr)

# We run the gnuplot scripts we wants
scripts = glob.glob("*."+GNUPLOT_EXTENSION)
scripts.sort()
nb_scripts = len(scripts)

# If the action is not given in parameters, we ask the user through standard input
if (indice_script == None):
  enumeration(scripts)

  while not(0<=indice_script<nb_scripts):
    try:
      txt_input = raw_input("What test do you want to display? (0-"+str(nb_scripts-1)+" ; 'all' treat them all ; 'l' display list again)\n")
      
      indice_script = int(txt_input)
    except ValueError:
      if (txt_input == 'l'):
        enumeration(scripts)
      elif (txt_input == 'all'):
        generate_all = True
        indice_script = 0
      else:
        print("The parameter must be between 0 and %i" % (nb_scripts-1))

if not(generate_all):
  indexes = [indice_script]
else:
  indexes = range(nb_scripts)

# We run gnuplot for every script needed. 
for index in indexes:
  script = scripts[index]
  sys.stdout.write("Running %s...\r" % script)
  sys.stdout.flush()
  (stdout, stderr, returncode) = run("gnuplot "+script)
  sys.stdout.write("Running %s...OK\n" % script)
  sys.stdout.flush()

# We only display the output file for the last script (in order to display only one graphic if we choose to run them all
output_file = getOutput(script)
if (output_file != None):
  output_extension = os.path.splitext(output_file)[1][1:] # [1:] is here to get rid of the prefixed dot  "." of the extension

  run(VISUALISEUR[output_extension]+" "+output_file)
else:
  print("Warning: The current gnuplot file has no output file")