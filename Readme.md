### The Tenstream solver approximately solves the radiative transfer equation in 3D and computes irradiance and heating rates from optical properties.

  At the moment the conversion from physical variables like liquid water
  to optical properties is not included.

#### The solver is coupled to the software packages:
  * LibRadtran  -- all purpose radiative transfer Library(partially free)
  * UCLA-LES    -- Large Eddy Simulation code for cloud resolving simulations

##### Prerequisites for the Tenstream solver are
  * cmake
  * PETSc
  * MPI
 Instructions on how to install these is beyond the scope of this text.

#### The solver is currently tested and running at:
* Linux Cluster University Munich
* Linux Cluster LRZ
* Linux Cluster MaxPlanckInstitute Hamburg
* IBM Power6  Machine "blizzard" at DKRZ
 
You may find a hint for a suitable config file in the config/ directory:

#### Currently the solver is tested for Compilers:
* GFortran
* Intel Fortran Compiler
* XLF IBM Compilers

 If you found yourself go nuts because of compile errors,
 please consider providing installation instructions
 for your particular enviroment.

#### Optional steps before you start:
```
 You may want to override the lut_basename, the path where the tenstream solver will put the lookuptables
 the default is otherwise local job directory. -- see src/optprop_parameters.f90
 lut_basename should be a path, that is reachable for mpi_rannk==0.
 and if possible globally reachable so that we dont have to compute the lookuptables over and over again.
 The lut_basename may also be provided as a commandline option: -lut_basename /global_reachable_path/LUT
 or simply put in a ~/.petscrc or ./.petscrc
```
>echo "-lut_basename /global_reachable_path/LUT" >> ~/.petscrc

 Next step is to set the paths for additional libraries and thus provide a config file
 Check the appropriate files in config/*

##### Then build tenstream
>cd $TENSTREAM_DIR
>mkdir build && cd build
>cmake .. && make -j

Check the tenstream runtime options with
>./bin/petsc_solver -show_options

For a first example you could run
>mpirun -np 2 ./bin/createLUT_8_10 -ident run_test -dx 70 -sza 60

Note that creating the Lookuptables takes a loong time.... dont be discouraged if you run it for the first time.

If that does not throw a lot of errors ---
well, then... Congratulations !!, you successfully installed the tenstream solver.

You may head back to the installation of the host model.

