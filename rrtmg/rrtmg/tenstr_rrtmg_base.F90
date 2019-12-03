module m_tenstr_rrtmg_base
#include "petsc/finclude/petsc.h"
  use petsc
  use m_data_parameters, only : iintegers, mpiint, default_str_len, i0
  use m_helper_functions, only : CHKERR, get_arg
  implicit none

  type t_rrtmg_log_events
    PetscLogStage :: stage_rrtmg_solar
    PetscLogStage :: stage_rrtmg_thermal
    PetscLogEvent :: setup_tenstr_atm
    PetscLogEvent :: rrtmg_optprop_lw
    PetscLogEvent :: rrtmg_optprop_sw
    PetscLogEvent :: smooth_surface_fluxes
  end type

  private
  public :: t_rrtmg_log_events, setup_log_events

  contains
    subroutine setup_log_events(logs, solvername)
      type(t_rrtmg_log_events), intent(inout) :: logs
      character(len=*), optional :: solvername
      character(len=default_str_len) :: s
      PetscClassId, parameter :: cid=0
      integer(mpiint) :: ierr

      s = get_arg('pprts.', solvername)

      call setup_stage(trim(s)//'rrtmg_solar'  , logs%stage_rrtmg_solar  )
      call setup_stage(trim(s)//'rrtmg_thermal', logs%stage_rrtmg_thermal)

      call PetscLogEventRegister(trim(s)//'setup_tenstr_atm', cid, logs%setup_tenstr_atm, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'rrtmg_optprop_lw', cid, logs%rrtmg_optprop_lw, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'rrtmg_optprop_sw', cid, logs%rrtmg_optprop_sw, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(trim(s)//'smooth_surface_fluxes', cid, logs%smooth_surface_fluxes, ierr); call CHKERR(ierr)

      contains
        subroutine setup_stage(stagename, logstage)
          character(len=*), intent(in) :: stagename
          PetscLogStage, intent(inout) :: logstage
          call PetscLogStageGetId(stagename, logstage, ierr); call CHKERR(ierr)
          if(logstage.lt.i0) then
            call PetscLogStageRegister(stagename, logstage, ierr); call CHKERR(ierr)
          endif
        end subroutine
    end subroutine
end module
