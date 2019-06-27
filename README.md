# PIMAIM

## Compilation

There is a helper script `compile_pimaim` for compiling the code. This uses `makefile_gww`
This sorts out most of the conditional compilation.
e.g. to compile a dipole version:
```
./compile_pimaim dipole
```
To clean (remove) the dipole version
```
./compile_pimaim dipole clean
```
Object code is stored in a subdirectory so you can have a number of compilations present in the directory at the same time.



Changes in this branch

Modified output of paramters for potential fitting - triggered by setting verbose true in runtime.inpt 

outputfiles 

forces.out, dipoles.out quads.out stresses.out 

