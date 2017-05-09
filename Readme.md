### The Tenstream solver approximately solves the radiative transfer equation in 3D and computes irradiance and heating rates from optical properties.

[![Documentation](https://codedocs.xyz/tenstream/tenstream.svg)](https://codedocs.xyz/tenstream/tenstream/)

  Overview Documentation is available in the corresponding paper:
  * A three-dimensional parallel radiative transfer model for atmospheric heating rates for use in cloud resolving models—The TenStream solver
  <http://dx.doi.org/10.1016/j.jqsrt.2015.05.003>
  * 3-D Radiative Transfer in Large-Eddy Simulations – Experiences Coupling the TenStream Solver to the UCLA–LES
  <http://dx.doi.org/10.5194/gmd-9-1413-2016>

#### The solver is coupled to the software packages:
  * [LibRadtran](http://www.libradtran.org)  -- all purpose radiative transfer Library(partially free)
  * [UCLA-LES](http://www.github.com/uclales/uclales)    -- Large Eddy Simulation code for cloud resolving simulations
  * [COSMO](http://www.cosmo-model.org)       -- Numerical Weather Prediction model
  * [DALES](https://github.com/dalesteam/dales) -- Dutch Atmospheric Large-Eddy Simulation model 

#### Note concerning the usage
The code is distributed under the GPL, and you are therefore free to use, change and redistribute it.  
I do however highly encourage you to participate in the development of the codebase.
If you are using the code in your work, please consider sharing bugfixes and experiences.
Given the experimental status of the solver, I kindly ask that you get in touch before publishing any results concerning the TenStream solver to ensure correctness of the results.  
It would also be appreciated to discuss co-authorship for research publications conducted with the TenStream solver.
#### Contact
Don't hesitate to ask, fabian@jakub.com (<a href="http://jakub.com/" target="_blank">www.jakub.com</a>
)

---

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



