### The Tenstream solver approximately solves the radiative transfer equation in 3D and computes irradiance and heating rates from optical properties.

  Overview Documentation is available in the corresponding paper:
  * A three-dimensional parallel radiative transfer model for atmospheric heating rates for use in cloud resolving models—The TenStream solver
  <http://dx.doi.org/10.1016/j.jqsrt.2015.05.003>

#### The solver is coupled to the software packages:
  * [LibRadtran](http://www.libradtran.org)  -- all purpose radiative transfer Library(partially free)
  * [UCLA-LES](http://www.github.com/uclales/uclales)    -- Large Eddy Simulation code for cloud resolving simulations
  * [COSMO](http://www.cosmo-model.org)       -- Numerical Weather Prediction model 

##### Further installation instructions and quick tips are available at the [TenStream-Wiki](https://github.com/tenstream/tenstream/wiki)

##### Prerequisites for the Tenstream solver are
  * cmake
  * PETSc
  * MPI

 Instructions on how to install these is beyond the scope of this text.

#### The solver is currently tested and running at:
* Linux Cluster University Munich
* Linux Cluster LRZ
* Linux Cluster MaxPlanckInstitute Hamburg
* IBM Power6  Machine "blizzard" at DKRZ, Hamburg
* Mistral supercomputer, Intel Haswell, at DKRZ, Hamburg
 
You may find a hint for a suitable config file in the config/ directory:

#### Currently the solver is tested for Compilers:

* GFortran
* Intel Compiler
* XLF IBM Compilers

If you found yourself go nuts because of compile errors,
please consider providing installation instructions
for your particular enviroment.



