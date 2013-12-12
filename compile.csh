gfortran -fbounds-check -ffast-math nls_header_mod.f90
gfortran -fbounds-check -ffast-math -o nautilus opk*.f nautilus.f90 nls*.f90 
