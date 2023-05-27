!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use encap
  use pf_my_sweeper
  use probin
  use global_state
  implicit none

  interface
     function HypreMaxErr(x, t, init_cond) result(max_err) bind(c, name="HypreMaxErr")
        use iso_c_binding
        type(c_ptr), value :: x
        real(c_double), value :: t
        real(c_double), value :: init_cond
        real(c_double) :: max_err
     end function
  end interface
contains


  !>  Resize libpfasst randomly
  subroutine resize_decider(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    !integer :: max_timesteps = 8
    integer :: max_timesteps = 6
    integer :: cur_timesteps
    integer :: new_timesteps
    real :: u

    ! we only can set resize_delta at the process that calls the psetop
    if (pf%rank == 0 .and. ((.not. pf%dynprocs%global_used) .or. pf%dynprocs%horizontal_rank == 0)) then
        cur_timesteps = pf%comm%nproc
        ! get random number between 1 and max_timesteps
        ! and subtract cur_timesteps from it
        call random_number(u)
        new_timesteps = 1 + floor(u * (max_timesteps +1 - 1))
        print *, "Trying to resize to ", new_timesteps, " parallel timesteps"
        pf%dynprocs%resize_delta = new_timesteps - cur_timesteps
        print *, "Set resize_delta to ", pf%dynprocs%resize_delta
    end if
  end subroutine resize_decider

  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    real(pfdp) :: yexact
    real(pfdp) :: maxerr, error 
    real(pfdp) :: residual
    class(hypre_vector_encap), pointer :: y_end
    integer :: ierr

    !> Get the solution at the end of this step
    y_end => cast_as_hypre_vector(pf%levels(level_index)%qend)

    !>  compute error
    error = HypreMaxErr(y_end%c_hypre_vector_ptr, Tfin, init_cond)
    residual = pf%levels(level_index)%residual

    pf%results%residuals(level_index,pf%state%pfblock,pf%state%iter+1,pf%state%sweep) = residual
    pf%results%errors(level_index,pf%state%pfblock,pf%state%iter+1,pf%state%sweep) = error


    if ((pf%rank == pf%comm%nproc-1) .and. (level_index == pf%nlevels) .and. (pf%state%iter .gt. 0) &
        .and. ((pf%state%step .eq. pf%state%nsteps-1))) then
       print '("error: rank: ", i4.4," step: ",i4.4," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es18.10e4)', &
            pf%rank,pf%state%step+1, pf%state%iter,level_index, error, residual
       call flush(6)
    end if
  end subroutine echo_error


  subroutine dump_hook(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    type(hypre_vector_encap) :: values
    integer            :: space_color, time_color
    character(len=1000) :: dump_dir
    character(len=1000) :: fname

    call get_global_str("dump_dir", dump_dir)
    call get_global_int("time_color", time_color)
    call get_global_int("space_color", space_color)

    write(fname, "(A,i5.5,A,i4.4,A,i4.4,A,i4.4,A)") &
         trim(adjustl(dump_dir)) // trim("/dump_step"), pf%state%step+1, &
         "_time", time_color, "_space", space_color, "_level", level_index, ".csv"

    values = cast_as_hypre_vector(pf%levels(level_index)%qend)
    print *, "Dumping to ", trim(fname)
    call values%dump(fname)
  end subroutine dump_hook


  subroutine print_sol(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    real(pfdp) :: yexact
    real(pfdp) :: maxerr, error
    real(pfdp) :: residual
    class(hypre_vector_encap), pointer :: q0, qend
    integer :: nproc, rank, ierr

    !> Get the solution at the end of this step
    q0 => cast_as_hypre_vector(pf%levels(level_index)%q0)
    qend => cast_as_hypre_vector(pf%levels(level_index)%qend)

    print *, "Q0"
    call q0%eprint()

    print *, "QEnd"
    call qend%eprint()
  end subroutine print_sol

  subroutine print_sol_encap(q)
    class(pf_encap_t), intent(in) :: q

    class(hypre_vector_encap), pointer :: qhr

    !> Get the solution at the end of this step
    qhr => cast_as_hypre_vector(q)
    call qhr%eprint()
  end subroutine print_sol_encap

end module hooks
