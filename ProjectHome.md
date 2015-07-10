A fast 1D gas-grain chemical model by FH (2008). Based upon the OSU gas-grain chemical model. Updated from gg\_osu\_2006v1d by RTG/VW.
Rate equations from Hasegawa & Herbst (1992). Modified rates following Caselli et al. (1998)\n\n
Stiff solver for sparse Jacobians: LSODES/ODEPACK (Hindmarsh 1983)\n
Turbulent mixing implemented through first order operator splitting\n

# Documentation #
For user documentation, please refer to [PDF documentation](https://nautilus.googlecode.com/git/DOC/nautilus_documentation.pdf)

For developper information, check the [Doxygen documentation](http://nautilus.googlecode.com/git/html/index.html).

# How To Use #
A folder containing an example simulation is available in example\_simulation/

To use the code, my advice would be to have one folder for the code, in a specific location. Then, in another location, you can
create as many folders as you want, with one specific simulation in each folder.

Be sure to run the installation script (to add useful shortcuts to your .bash\_profile):
```
configure_nautilus.py
```

You can jump to your program directory via :
```
cd $nautilus
```

To launch Nautilus wherever you want in your server via :
```
nautilus_code
```
**Be careful, 'nautilus' will not launch the code, but rather the file browser in a GNU/Linux environement**

Keep in mind that typing
```
> cd -
```
in your Terminal, you will go back to the previous directory (not the parent one, the previous one). Typing this command several
times allow you to switch easily between two folders.

# Basic commands #

Compile Nautilus
```
Makefile.py
```

See all compilation options
```
Makefile.py help
```

Generate Doxygen documentation (all errors will be stored in a <b>doxygen.log</b> file)
```
doxygen doxygen.conf 2>doxygen.log
```

# Warnings #
A lot of warnings are launched by the ODEPACK files. Some of them are more problematic than others that are just deprecated feature
(even though this is a problem for readability).

The most problematic is a call of DPREP in the subroutine DIPREP, in opkda1.

The error returned is the following:
```
LYH), RWORK(LSAVF), RWORK(LEWT), RWORK(LACOR), IA, JA, IC, JC, RWORK(LWM), RWOR
                                                                           1                              
Error: Type mismatch in argument 'iwk' at (1); passed REAL(8) to INTEGER(4)
opkda1.f90:1214.93:

, RWORK(LYH), RWORK(LSAVF), RWORK(LEWT), RWORK(LACOR), IA, JA, RWORK(LWM), RWOR
                                                                           1                  
Error: Type mismatch in argument 'iwk' at (1); passed REAL(8) to INTEGER(4)
```