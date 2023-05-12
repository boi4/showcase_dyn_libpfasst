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

  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    real(pfdp) :: yexact
    real(pfdp) :: maxerr, error 
    real(pfdp) :: residual
    class(hypre_vector_encap), pointer :: y_end
    integer :: nproc, rank, ierr

    !> Get the solution at the end of this step
    y_end => cast_as_hypre_vector(pf%levels(level_index)%qend)

    !>  compute error
    error = HypreMaxErr(y_end%c_hypre_vector_ptr, Tfin, init_cond)
    residual = pf%levels(level_index)%residual

    pf%results%residuals(level_index,pf%state%pfblock,pf%state%iter+1,pf%state%sweep) = residual
    pf%results%errors(level_index,pf%state%pfblock,pf%state%iter+1,pf%state%sweep) = error

    !call mpi_comm_rank(pf%comm%comm, rank, ierr)
    !call mpi_comm_size(pf%comm%comm, nproc, ierr)
    

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
    print *, "Dumping, space color=", space_color, " time color=", time_color

    write(fname, "(A,A,i5.5,A,i4.4,A,i4.4,A,i4.4,A)") &
         trim(adjustl(dump_dir)), "/dump_step", pf%state%step+1, &
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


  ! subroutine print_(pf, level_index)
  !   type(pf_pfasst_t), intent(inout) :: pf
  !   integer, intent(in) :: level_index

  !   class(hypre_vector_encap), pointer :: qhr
  !   integer :: j


  !   do j = 1, pf%levels(1)%nnodes-1
  !     qhr => cast_as_hypre_vector(pf%levels(1)%I(j))
  !     !> Get the solution at the end of this step
  !     print *, j
  !     call qhr%eprint()
  !   end do
  !     ! do j = 1, pf%levels(1)%nnodes-1
  !     !   qhr => cast_as_hypre_vector(pf%levels(1)%I(j))
  !     !   !> Get the solution at the end of this step
  !     !   print *, j
  !     !   call qhr%eprint()
  !     ! end do
  ! end subroutine print_level_vals

  subroutine print_level_valshelper(sweeper)
    class(pf_sweeper_t), target, intent(in) :: sweeper

    class(hypre_vector_encap), pointer :: qhr
    type(my_sweeper_t), pointer :: my_sweeper

    print *, "rhs"
    select type(sweeper)
    type is (my_sweeper_t)
      my_sweeper => sweeper
    end select
    qhr => cast_as_hypre_vector(my_sweeper%rhs)
    call qhr%eprint()
  end subroutine print_level_valshelper

  subroutine print_level_vals(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    class(hypre_vector_encap), pointer :: qhr

    print *, "Q(1)"
    qhr => cast_as_hypre_vector(pf%levels(level_index)%Q(1))
    call qhr%eprint()
    print *, "Q(2)"
    qhr => cast_as_hypre_vector(pf%levels(level_index)%Q(2))
    call qhr%eprint()
    print *, "F(1,1)"
    qhr => cast_as_hypre_vector(pf%levels(level_index)%F(1,1))
    call qhr%eprint()
    print *, "F(1,2)"
    qhr => cast_as_hypre_vector(pf%levels(level_index)%F(1,2))
    call qhr%eprint()
    print *, "F(2,1)"
    qhr => cast_as_hypre_vector(pf%levels(level_index)%F(2,1))
    call qhr%eprint()
    print *, "F(2,2)"
    qhr => cast_as_hypre_vector(pf%levels(level_index)%F(2,2))
    call qhr%eprint()

    call print_level_valshelper(pf%levels(level_index)%ulevel%sweeper)

  end subroutine print_level_vals


  subroutine pstep(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    print *,pf%state%step
  end subroutine pstep

end module hooks
