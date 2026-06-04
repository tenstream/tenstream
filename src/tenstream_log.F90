module m_tenstream_log
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_data_parameters, only: mpiint, default_str_len

  implicit none
  private

  public :: &
    & t_ts_log_event, &
    & t_ts_log_stage, &
    & ts_log_event_register, &
    & ts_log_stage_register, &
    & ts_log_begin, &
    & ts_log_end, &
    & ts_log_stage_push, &
    & ts_log_stage_pop, &
    & ts_log_view

  integer, parameter :: MAX_EVENTS = 512
  ! Stage index 1 is always "Main Stage"; user-registered stages occupy 2..MAX_STAGES.
  integer, parameter :: MAX_STAGES = 64

  type t_ts_log_event
    PetscLogEvent :: petsc_id = -1
    integer :: ts_id = -1
  end type

  type t_ts_log_stage
    PetscLogStage :: petsc_id = -1
    integer :: ts_id = -1
  end type

  integer, save :: n_events = 0
  integer, save :: n_stages = 1   ! 1 = Main Stage always present
  integer, save :: current_stage = 1

  character(len=default_str_len), save :: event_names(MAX_EVENTS)
  character(len=default_str_len), save :: stage_names(MAX_STAGES)

  real(8), save :: t_log_start = 0.0d0

  real(8), save :: event_accum(MAX_EVENTS) = 0.0d0
  real(8), save :: event_tstart(MAX_EVENTS) = 0.0d0
  integer, save :: event_count(MAX_EVENTS) = 0

  ! Stage at which each event's ts_log_begin was called; used by ts_log_end for attribution.
  integer, save :: event_stage_at_begin(MAX_EVENTS) = 1

  ! Per-stage accumulation: first index is stage (1=Main), second is event id.
  real(8), save :: event_stage_accum(MAX_STAGES, MAX_EVENTS) = 0.0d0
  integer, save :: event_stage_count(MAX_STAGES, MAX_EVENTS) = 0

  ! Actual wallclock time spent inside each stage (push→pop pairs).
  real(8), save :: stage_accum(MAX_STAGES) = 0.0d0
  real(8), save :: stage_tstart(MAX_STAGES) = 0.0d0

contains

  subroutine ts_log_event_register(name, event, ierr)
    character(len=*), intent(in) :: name
    type(t_ts_log_event), intent(inout) :: event
    integer(mpiint), intent(out) :: ierr
    integer :: i
    ierr = 0
    if (n_events == 0) then
      stage_names(1) = 'Main Stage'
      t_log_start = MPI_Wtime()
    end if
    do i = 1, n_events
      if (trim(event_names(i)) == trim(name)) then
        event%ts_id = i
        return
      end if
    end do
    n_events = n_events + 1
    if (n_events > MAX_EVENTS) then
      n_events = MAX_EVENTS
      ierr = 1_mpiint
      return
    end if
    event%ts_id = n_events
    event_names(n_events) = trim(name)
    event_accum(n_events) = 0.0d0
    event_tstart(n_events) = 0.0d0
    event_count(n_events) = 0
  end subroutine

  subroutine ts_log_stage_register(name, stage, ierr)
    character(len=*), intent(in) :: name
    type(t_ts_log_stage), intent(inout) :: stage
    integer(mpiint), intent(out) :: ierr
    integer :: i
    ierr = 0
    do i = 2, n_stages
      if (trim(stage_names(i)) == trim(name)) then
        stage%ts_id = i
        return
      end if
    end do
    n_stages = n_stages + 1
    if (n_stages > MAX_STAGES) then
      n_stages = MAX_STAGES
      ierr = 1_mpiint
      return
    end if
    stage%ts_id = n_stages
    stage_names(n_stages) = trim(name)
  end subroutine

  subroutine ts_log_begin(event, ierr)
    type(t_ts_log_event), intent(in) :: event
    integer(mpiint), intent(out) :: ierr
    integer :: id
    ierr = 0
    if (event%petsc_id >= 0) then
      call PetscLogEventBegin(event%petsc_id, ierr)
    end if
    id = event%ts_id
    if (id < 1 .or. id > n_events) return
    event_tstart(id) = MPI_Wtime()
    event_stage_at_begin(id) = current_stage
  end subroutine

  subroutine ts_log_end(event, ierr)
    type(t_ts_log_event), intent(in) :: event
    integer(mpiint), intent(out) :: ierr
    integer :: id, s
    real(8) :: dt
    ierr = 0
    if (event%petsc_id >= 0) then
      call PetscLogEventEnd(event%petsc_id, ierr)
    end if
    id = event%ts_id
    if (id < 1 .or. id > n_events) return
    dt = MPI_Wtime() - event_tstart(id)
    event_accum(id) = event_accum(id) + dt
    event_count(id) = event_count(id) + 1
    ! Use stage recorded at begin-time so events that span a stage boundary
    ! are attributed to the stage where they started, not where they ended.
    s = event_stage_at_begin(id)
    if (s >= 1 .and. s <= MAX_STAGES) then
      event_stage_accum(s, id) = event_stage_accum(s, id) + dt
      event_stage_count(s, id) = event_stage_count(s, id) + 1
    end if
  end subroutine

  subroutine ts_log_stage_push(stage, ierr)
    type(t_ts_log_stage), intent(in) :: stage
    integer(mpiint), intent(out) :: ierr
    ierr = 0
    call PetscLogStagePush(stage%petsc_id, ierr)
    if (stage%ts_id >= 2 .and. stage%ts_id <= n_stages) then
      current_stage = stage%ts_id
      stage_tstart(current_stage) = MPI_Wtime()
    end if
  end subroutine

  subroutine ts_log_stage_pop(ierr)
    integer(mpiint), intent(out) :: ierr
    real(8) :: t
    ierr = 0
    call PetscLogStagePop(ierr)
    if (current_stage >= 2 .and. current_stage <= n_stages) then
      t = MPI_Wtime()
      stage_accum(current_stage) = stage_accum(current_stage) + t - stage_tstart(current_stage)
    end if
    current_stage = 1
  end subroutine

  subroutine ts_log_view(comm, ierr)
    integer(mpiint), intent(in) :: comm
    integer(mpiint), intent(out) :: ierr

    integer :: ne, ns, myid, numnodes, i, s
    integer(mpiint) :: mpierr

    real(8), allocatable :: tlocal(:), tmin(:), tmean(:), tmax(:)
    integer, allocatable :: cnt_local(:), cnt_global(:)
    ! Flattened stage×event arrays for reduction (ns*ne elements).
    real(8), allocatable :: st_local(:), st_min(:), st_mean(:), st_max(:)
    integer, allocatable :: sc_local(:), sc_global(:)

    ! Stage wallclock reductions (actual push→pop time, not sum of events).
    real(8) :: sw_local(MAX_STAGES), sw_max(MAX_STAGES)
    real(8) :: stage_total_max(MAX_STAGES), grand_max
    real(8) :: t_elapsed_local, t_elapsed_max
    character(len=48) :: sfmt

    ierr = 0
    if (n_events == 0) return

    ne = n_events
    ns = n_stages

    allocate (tlocal(ne), tmin(ne), tmean(ne), tmax(ne))
    allocate (cnt_local(ne), cnt_global(ne))

    tlocal = event_accum(1:ne)
    cnt_local = event_count(1:ne)

    call MPI_Comm_rank(comm, myid, mpierr)
    call MPI_Comm_size(comm, numnodes, mpierr)

    call MPI_Reduce(tlocal, tmin, ne, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, mpierr)
    call MPI_Reduce(tlocal, tmean, ne, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, mpierr)
    call MPI_Reduce(tlocal, tmax, ne, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, mpierr)
    call MPI_Reduce(cnt_local, cnt_global, ne, MPI_INTEGER, MPI_MAX, 0, comm, mpierr)

    ! Per-stage reductions (flattened to 1D).
    ! Layout: row-major (stage outer, event inner) to match (s-1)*ne+i access below.
    ! A plain reshape() would give column-major order, which is wrong here.
    allocate (st_local(ns * ne), st_min(ns * ne), st_mean(ns * ne), st_max(ns * ne))
    allocate (sc_local(ns * ne), sc_global(ns * ne))

    do s = 1, ns
      st_local((s - 1) * ne + 1:s * ne) = event_stage_accum(s, 1:ne)
      sc_local((s - 1) * ne + 1:s * ne) = event_stage_count(s, 1:ne)
    end do

    call MPI_Reduce(st_local, st_min, ns * ne, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, mpierr)
    call MPI_Reduce(st_local, st_mean, ns * ne, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, mpierr)
    call MPI_Reduce(st_local, st_max, ns * ne, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, mpierr)
    call MPI_Reduce(sc_local, sc_global, ns * ne, MPI_INTEGER, MPI_MAX, 0, comm, mpierr)

    ! Stage wallclock totals: user stages use push→pop time; main stage sums its events.
    sw_local = stage_accum(1:MAX_STAGES)
    sw_local(1) = sum(event_stage_accum(1, 1:ne))   ! main stage has no push/pop
    call MPI_Reduce(sw_local, sw_max, MAX_STAGES, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, mpierr)

    t_elapsed_local = MPI_Wtime() - t_log_start
    call MPI_Reduce(t_elapsed_local, t_elapsed_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, mpierr)

    if (myid /= 0) return

    tmean = tmean / real(numnodes, 8)
    st_mean = st_mean / real(numnodes, 8)

    do s = 1, ns
      stage_total_max(s) = sw_max(s)
    end do
    grand_max = t_elapsed_max
    if (grand_max <= 0.0d0) grand_max = 1.0d0

    write (*, '(a)') ''
    write (*, '(a)') 'TenStream Timing Summary (wallclock [s] per MPI rank)'
    write (*, '(a,i0)') '  Ranks: ', numnodes

    ! Stage summary header (only if user stages exist).
    if (ns > 1) then
      write (*, '(a)') ''
      write (*, '(a)') 'Summary of Stages:             ----- Time (s) -----   %Total'
      write (*, '(a)') repeat('-', 65)
      do s = 1, ns
        write (*, '(i2,a,a32,3x,f10.4,4x,f5.1,a)') &
          s - 1, ': ', trim(stage_names(s)), &
          stage_total_max(s), &
          100.0d0 * stage_total_max(s) / grand_max, '%'
      end do
      write (*, '(a)') repeat('-', 65)
      write (*, '(a37,3x,f10.4,4x,f5.1,a)') 'Total', grand_max, 100.0d0, '%'
    end if

    ! Per-stage event tables.
    do s = 1, ns
      write (*, '(a)') ''
      if (ns > 1) then
        write (sfmt, '(a,i0,a)') '--- Stage ', s - 1, ': '//trim(stage_names(s))//' ---'
        write (*, '(a)') trim(sfmt)
      end if
      write (*, '(a)') repeat('-', 80)
      write (*, '(a40,a7,3a12)') 'Event', 'Count', 'Min', 'Mean', 'Max'
      write (*, '(a)') repeat('-', 80)
      do i = 1, ne
        if (sc_global((s - 1) * ne + i) == 0) cycle
        write (*, '(a40,i7,3f12.4)') &
          trim(event_names(i)), sc_global((s - 1) * ne + i), &
          st_min((s - 1) * ne + i), st_mean((s - 1) * ne + i), &
          st_max((s - 1) * ne + i)
      end do
      write (*, '(a)') repeat('-', 80)
    end do

  end subroutine

end module m_tenstream_log
