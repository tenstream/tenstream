module m_tenstr_rrtmg_base
#ifdef HAVE_PETSC
#include "petsc/finclude/petsc.h"
  use petsc
#else
#define PetscLogStage integer
#define PetscLogEvent integer
#endif
  use m_data_parameters, only: iintegers, ireals, mpiint, default_str_len, i0
  use m_helper_functions, only: CHKERR, get_arg
  implicit none

  type t_rrtmg_log_events
    PetscLogStage :: stage_rrtmg_solar = 0
    PetscLogStage :: stage_rrtmg_thermal = 0
    PetscLogEvent :: setup_tenstr_atm = 0
    PetscLogEvent :: rrtmg_optprop_lw = 0
    PetscLogEvent :: rrtmg_optprop_sw = 0
    PetscLogEvent :: smooth_surface_fluxes = 0
  end type

  private
  public :: t_rrtmg_log_events, setup_log_events

contains
  subroutine setup_log_events(logs, solvername)
    type(t_rrtmg_log_events), intent(inout) :: logs
    character(len=*), optional :: solvername
    character(len=default_str_len) :: s
#ifdef HAVE_PETSC
    PetscClassId :: cid
#endif
    integer(mpiint) :: ierr

    s = get_arg('tenstr_rrtmg.', solvername)
#ifdef HAVE_PETSC
    call PetscClassIdRegister(trim(s), cid, ierr); call CHKERR(ierr)

    call setup_stage(trim(s)//'rrtmg_solar', logs%stage_rrtmg_solar)
    call setup_stage(trim(s)//'rrtmg_thermal', logs%stage_rrtmg_thermal)

    call PetscLogEventRegister(trim(s)//'setup_tenstr_atm', cid, logs%setup_tenstr_atm, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'rrtmg_optprop_lw', cid, logs%rrtmg_optprop_lw, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'rrtmg_optprop_sw', cid, logs%rrtmg_optprop_sw, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'smooth_surface_fluxes', cid, logs%smooth_surface_fluxes, ierr); call CHKERR(ierr)
#endif

  contains
    subroutine setup_stage(stagename, logstage)
      character(len=*), intent(in) :: stagename
      PetscLogStage, intent(inout) :: logstage
#ifdef HAVE_PETSC
      call PetscLogStageGetId(stagename, logstage, ierr); call CHKERR(ierr)
      if (logstage .lt. i0) then
        call PetscLogStageRegister(stagename, logstage, ierr); call CHKERR(ierr)
      end if
#endif
    end subroutine
  end subroutine

end module
