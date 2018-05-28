module m_pprts_base
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, mpiint, &
    zero, one, i0, i1, i2, i3, i4, i5, i6, i7, i8, i10, pi, &
    default_str_len

  use m_optprop, only: t_optprop

  public :: t_solver, t_solver_1_2, t_solver_3_6, t_solver_8_10, t_solver_3_10, &
    t_state_container, t_coord, t_opticalprops, t_sunangles, t_suninfo, &
    t_dof, t_solver_log_events, E_up, E_dn

  type t_coord
    integer(iintegers)      :: xs,xe                   ! local domain start and end indices
    integer(iintegers)      :: ys,ye                   ! local domain start and end indices
    integer(iintegers)      :: zs,ze                   ! local domain start and end indices
    integer(iintegers)      :: xm,ym,zm                ! size of local domain
    integer(iintegers)      :: gxs,gys,gzs             ! domain indices including ghost points
    integer(iintegers)      :: gxe,gye,gze             !
    integer(iintegers)      :: gxm,gym,gzm             ! size of local domain including ghosts

    integer(iintegers)      :: glob_xm,glob_ym,glob_zm ! global domain size

    integer(iintegers)      :: dof,dim                 ! degrees of freedom of Petsc Domain, dimension of dmda
    type(tDM)               :: da                      ! The Domain Decomposition Object
    PetscMPIInt,allocatable :: neighbors(:)            ! all 3d neighbours((x=-1,y=-1,z=-1), (x=0,y=-1,z=-1) ...), i.e. 14 is one self.
    integer(mpiint)         :: comm                    ! mpi communicatior for this DMDA
  end type

  type t_opticalprops
    real(ireals) :: kabs,ksca,g
  end type

  type t_atmosphere
    type(t_opticalprops) , allocatable , dimension(:,:,:) :: op
    real(ireals)         , allocatable , dimension(:,:,:) :: planck
    real(ireals)         , allocatable , dimension(:,:,:) :: a11, a12, a21, a22, a13, a23, a33
    real(ireals)         , allocatable , dimension(:,:,:) :: g1,g2
    real(ireals)         , allocatable , dimension(:,:,:) :: dz
    logical              , allocatable , dimension(:,:,:) :: l1d
    real(ireals)         , allocatable , dimension(:,:)   :: albedo
    real(ireals)                                          :: dx,dy
    integer(iintegers)                                    :: icollapse=1
    logical                                               :: lcollapse = .False.
  end type

  type t_sunangles
    real(ireals)        :: symmetry_phi
    integer(iintegers)  :: yinc,xinc
    real(ireals)        :: theta, phi, costheta, sintheta
  end type

  type t_suninfo
    type(t_sunangles),allocatable :: angles(:,:,:) ! defined on DMDA grid
    logical                       :: luse_topography=.False.
  end type

  type t_state_container
    integer(iintegers)  :: uid ! dirty hack to give the solution a unique hash for example to write it out to disk -- this should be the same as the index in global solutions array
    Vec                 :: edir,ediff,abso

    logical             :: lset        = .False. ! initialized?
    logical             :: lsolar_rad  = .False. ! direct radiation calculated?
    logical             :: lchanged    = .True.  ! did the flux change recently? -- call restore_solution to bring it in a coherent state

    ! save state of solution vectors... they are either in [W](true) or [W/m**2](false)
    logical             :: lintegrated_dir=.True. , lintegrated_diff=.True.

    !save error statistics
    real(ireals)        :: time   (30) = -one
    real(ireals)        :: maxnorm(30) = zero
    real(ireals)        :: twonorm(30) = zero
    real(ireals),allocatable :: ksp_residual_history(:)
  end type

  type t_dof
    integer(iintegers) :: dof
    logical, allocatable :: is_inward(:)
  end type

  type t_solver_log_events
    PetscLogStage :: stage_solve_pprts
    PetscLogEvent :: set_optprop
    PetscLogEvent :: compute_edir
    PetscLogEvent :: solve_Mdir
    PetscLogEvent :: setup_Mdir
    PetscLogEvent :: compute_ediff
    PetscLogEvent :: solve_Mdiff
    PetscLogEvent :: setup_Mdiff
    PetscLogEvent :: solve_twostream
    PetscLogEvent :: solve_mcrts
    PetscLogEvent :: get_coeff_dir2dir
    PetscLogEvent :: get_coeff_dir2diff
    PetscLogEvent :: get_coeff_diff2diff
    PetscLogEvent :: scatter_to_Zero
    PetscLogEvent :: scale_flx
  end type

  type, abstract :: t_solver
    integer(mpiint)                 :: comm, myid, numnodes     ! mpi communicator, my rank and number of ranks in comm
    type(t_coord), allocatable      :: C_dir, C_diff, C_one, C_one1, C_one_atm, C_one_atm1
    type(t_atmosphere),allocatable  :: atm
    type(t_suninfo)                 :: sun
    type(tMat),allocatable          :: Mdir,Mdiff
    class(t_optprop), allocatable   :: OPP

    type(t_dof)                     :: difftop, diffside, dirtop, dirside

    logical                         :: lenable_solutions_err_estimates=.True.  ! if enabled, we can save and load solutions.... just pass an unique identifer to solve()... beware, this may use lots of memory
    type(tVec),allocatable          :: incSolar,b

    logical                         :: linitialized=.False.
    type(t_state_container)         :: solutions(-1000:1000)
    type(t_solver_log_events)       :: logs
  end type

  type, extends(t_solver) :: t_solver_1_2
  end type
  type, extends(t_solver) :: t_solver_8_10
  end type
  type, extends(t_solver) :: t_solver_3_6
  end type
  type, extends(t_solver) :: t_solver_3_10
  end type


  integer(iintegers), parameter :: E_up=0, E_dn=1 ! for 1D Solvers

end module
