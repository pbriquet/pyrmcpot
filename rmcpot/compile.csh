#!/bin/tcsh -xvf

gfortran -g -fcheck=all -O2 -o rmcpot globals.f90 buckingham.f90 eamalloy.f90 meam.f90 tersoff.f90 watmodel.f90 potfunc.f90 datafunc.f90 localmin.f90 rmcpot.f90 >& compile.out

rm *.mod
