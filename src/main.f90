!
! This file is part of LIBPFASST.
!
!> Example of using LIBPFASST.
!!
!!  This program solves the linear hypre_vector equation
!!
!!    y'=lambda*y
!!
!>  The main program here just initializes mpi, calls the solver and then finalizes mpi
program main
  use pf_mod_mpi
  implicit none

  integer ::  ierr
  integer :: session

  !> Initialize MPI Session
  call mpi_session_init(MPI_INFO_NULL, MPI_ERRHANDLER_NULL, session, ierr)
  if (ierr /= 0) &
      stop "ERROR: Can't initialize MPI."

  !> Call the  solver 
  call run_pfasst(session)

  !> Close mpi
  call mpi_session_finalize(session, ierr)

contains
  !>  This subroutine set ups and calls libpfasst 
  subroutine run_pfasst(session)
    use pfasst            !< This module has include statements for the main pfasst routines
    use pf_my_sweeper     !< Local module for sweeper
    use pf_my_level       !< Local module for level
    use hooks             !< Local module for diagnostics and i/o
    use probin            !< Local module reading/parsing problem parameters
    use encap             !< Local module defining the encapsulation
    use pf_space_comm
    use pfasst_hypre
    use pf_mod_parareal
    use global_state

    implicit none

    !> argument
    integer, intent(in) :: session
    integer mcomm
    integer mgroup

    !>  Local variables
    type(pf_pfasst_t) :: pf        !<  the main pfasst structure
    type(pf_comm_t) :: comm      !<  the communicator (here it is mpi)
    type(hypre_vector_encap) :: y_0      !<  the initial condition
    type(hypre_vector_encap) :: y_end    !<  the solution at the final time
    type(hypre_vector_encap) :: y_other
    type(hypre_vector_encap) :: u
    class(pf_encap_t), allocatable :: y_0_base, y_end_base, u_base
    character(256) :: pf_fname   !<  file name for input of PFASST parameters
    integer, allocatable :: lev_shape(:,:)
    type(hypre_vector_factory) :: hvf
    type(my_sweeper_t) :: s_finest, s
    type(my_level_t) :: my_lev 

    integer :: l, l_finest   !  loop variable over levels
    integer :: n, m
    integer :: space_comm, time_comm, space_color, time_color
    integer :: level_index

    integer :: nproc, rank, error
    integer :: ierr
    real(pfdp) :: f
    integer :: nrows, ilower0, ilower1, iupper0, iupper1
    integer :: spacial_coarsen_flag
    type(mgrit_level_data), allocatable :: mg_ld(:)
    character(len=1000) :: fname

    integer :: tmp
    character(len=100) :: tmpstr

    character(len=MPI_MAX_PSET_NAME_LEN)  :: time_pset
    integer :: ntime

    ! TODO: time color and ntime

    !> Read problem parameters
    call probin_init(pf_fname)



    call mpi_group_from_session_pset(session, "mpi://WORLD", mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group from pset fail, error=',ierr)
    call mpi_comm_create_from_group(mgroup, "showcase", MPI_INFO_NULL, MPI_ERRORS_RETURN, mcomm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm create from group fail, error=',ierr)
    call mpi_group_free(mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierr)

    ! check size
    call mpi_comm_size(mcomm, nproc, error)
    call mpi_comm_rank(mcomm, rank,  error)

    print *,"nspace: ", nspace
    call create_pset_grid(session, "mpi://WORLD", mcomm, nspace, space_dim, &
         space_comm, space_color, time_comm, time_color, time_pset)

    !>  Set up communicator
    call pf_mpi_create(comm, time_comm)
    

    pf%use_rk_stepper = .false.
    pf%use_sdc_sweeper = .true.

    !>  Create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=pf_fname)

    spacial_coarsen_flag = 0
    call mpi_barrier(space_comm, ierr);
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi barrier fail, error=',ierr)
    call PfasstHypreInit(pf, mg_ld, lev_shape, space_comm, time_color, spacial_coarsen_flag)
    print *, "PfasstHypreInit done"
    print *,time_color,space_color,pf%rank

    ! hooks
    call pf_add_hook(pf, -1, PF_POST_BLOCK, echo_error)

    if (dump_values) then
        ! Save some global variables for the dump hook
        call set_global_str("dump_dir", dump_dir)
        call set_global_int("time_color", time_color)
        call set_global_int("space_color", space_color)

        ! TODO: what about spaces in dump_dir?
        call system('mkdir -p ' // trim(adjustl(dump_dir)))

        call pf_add_hook(pf, -1, PF_POST_BLOCK, dump_hook)
        !call pf_add_hook(pf, -1, PF_POST_PREDICTOR, dump_hook)
    end if

    !>  Output run parameters
    if (space_color == 0 .and. time_color == 0) then
        call print_loc_options(pf,un_opt=6)
    end if

    level_index = pf%nlevels
    !> Set the initial condition
    call hvf%create_single(y_0_base, level_index, lev_shape(pf%nlevels,:))
    call hvf%create_single(y_end_base, level_index, lev_shape(pf%nlevels,:))
    y_0 = cast_as_hypre_vector(y_0_base)
    y_end = cast_as_hypre_vector(y_end_base)
    call initial(y_0)

!    ! dump initial values
!    if (dump_values) then
!        write(fname, "(A,A,i5.5,A,i4.4,A,i4.4,A,i4.4,A)") &
!            trim(adjustl(dump_dir)), "/dump_step", 0, &
!            "_time", time_color, "_space", space_color, "_level", 2, ".csv"
!        call y_0%dump(fname)
!    end if
!
!    print *, "initial done"
!
!
!    !> Do the PFASST time stepping
!    call pf_pfasst_run(pf, y_0, dt, Tfin, nsteps, y_end)
!    if (pf%rank .eq. pf%comm%nproc-1) call y_end%eprint()
!
!    if (dump_values .and. space_color == ntime-1) then
!       write(fname, "(A, A,i1,A)") trim(adjustl(dump_dir)), "/final_dump_", time_color, ".csv"
!       call y_end%dump(fname)
!    end if

    !>  Wait for everyone to be done
    call mpi_barrier(mcomm, ierr)
    
    !>  Deallocate pfasst structure
    call pf_pfasst_destroy(pf)

    !> Free the communicator
    call mpi_comm_free(mcomm, ierr)

  end subroutine run_pfasst
end program
