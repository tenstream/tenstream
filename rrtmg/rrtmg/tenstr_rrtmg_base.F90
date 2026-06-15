module m_tenstr_rrtmg_base
#ifdef HAVE_PETSC
#include "petsc/finclude/petsc.h"
  use petsc
#endif

  use m_tenstream_log, only: &
    & t_ts_log_event, &
    & t_ts_log_stage, &
    & ts_log_event_register, &
    & ts_log_stage_register

  use m_data_parameters, only: mpiint, default_str_len, i0
  use m_helper_functions, only: CHKERR, get_arg
  implicit none

  type t_rrtmg_log_events
    type(t_ts_log_stage) :: stage_rrtmg_solar
    type(t_ts_log_stage) :: stage_rrtmg_thermal
    type(t_ts_log_event) :: setup_tenstr_atm
    type(t_ts_log_event) :: rrtmg_optprop_lw
    type(t_ts_log_event) :: rrtmg_optprop_sw
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
#endif

    call setup_stage(trim(s)//'rrtmg_solar', logs%stage_rrtmg_solar)
    call setup_stage(trim(s)//'rrtmg_thermal', logs%stage_rrtmg_thermal)

    call reg(trim(s)//'setup_tenstr_atm', logs%setup_tenstr_atm)
    call reg(trim(s)//'rrtmg_optprop_lw', logs%rrtmg_optprop_lw)
    call reg(trim(s)//'rrtmg_optprop_sw', logs%rrtmg_optprop_sw)
  contains
    subroutine reg(name, event)
      character(len=*), intent(in) :: name
      type(t_ts_log_event), intent(inout) :: event
      call ts_log_event_register(name, event, ierr); call CHKERR(ierr)
#ifdef HAVE_PETSC
      call PetscLogEventRegister(name, cid, event%petsc_id, ierr); call CHKERR(ierr)
#endif
    end subroutine

    subroutine setup_stage(stagename, logstage)
      character(len=*), intent(in) :: stagename
      type(t_ts_log_stage), intent(inout) :: logstage
      call ts_log_stage_register(stagename, logstage, ierr); call CHKERR(ierr)
#ifdef HAVE_PETSC
      call PetscLogStageGetId(stagename, logstage%petsc_id, ierr); call CHKERR(ierr)
      if (logstage%petsc_id .lt. i0) then
        call PetscLogStageRegister(stagename, logstage%petsc_id, ierr); call CHKERR(ierr)
      end if
#endif
    end subroutine
  end subroutine

end module
