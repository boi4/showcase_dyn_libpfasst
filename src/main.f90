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
  call mpi_session_init(MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, session, ierr)
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
    use pf_space_comm     !< Local module for communication
    use hooks             !< Local module for diagnostics and i/o
    use probin            !< Local module reading/parsing problem parameters
    use encap             !< Local module defining the encapsulation
    use pfasst_hypre
    use global_state

    implicit none

    !> argument
    integer, intent(in) :: session

    !>  Local variables
    type(pf_pfasst_t) :: pf              !<  the main pfasst structure
    type(pf_dynprocs_t) :: dynprocs      !<  the dynprocs object
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

    real(pfdp) :: nspace_real
    integer :: nproc, rank, error
    integer :: ierr
    real(pfdp) :: f
    integer :: nrows, ilower0, ilower1, iupper0, iupper1
    integer :: spacial_coarsen_flag
    character(len=1000) :: fname
    integer :: ntime
    character(len=MPI_MAX_PSET_NAME_LEN)  :: space_pset
    character(len=MPI_MAX_PSET_NAME_LEN)  :: time_pset
    logical :: is_dynamic
    logical :: premature_exit
    integer :: tmp_int
    character(len=*), parameter :: ESC = achar(27)
    character(len=9) :: COLORS(10) =[&
"[0;31m", &
"[0;32m", &
"[0;33m", &
"[0;34m", &
"[0;35m", &
"[0;36m", &
"[0;91m", &
"[0;92m", &
"[0;93m", &
"[0;94m" &
]

    !> Read problem parameters
    call probin_init(pf_fname)

    ! for now, nspace must be a perfect square
    if (space_dim .eq. 2) then
       nspace_real = sqrt(real(nspace))
       if (nspace_real .ne. nint(nspace_real)) then
          print'(a)', 'ERROR: create_simple_communicators: nspace must be perfect square.'
          stop
       end if
    end if


    ! split up processes into grid
    call create_pset_grid(session, "mpi://WORLD", nspace, space_dim, &
                          space_comm, space_color, space_pset, time_comm, time_color, time_pset)

    if (color_output) then
        write (*, '(a)', advance='no') ESC // COLORS(time_color+1)
    end if

    call mpi_comm_size(time_comm, ntime, ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi comm size fail, error=',ierr)

    call mpi_comm_free(time_comm, ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi comm free fail, error=',ierr)

    call mpi_barrier(space_comm, ierr);
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi barrier fail, error=',ierr)



    ! determine if dynamic
    call pf_dynprocs_check_dynamic(session, is_dynamic)

    !>  Set up communicator
    call pf_dynprocs_create(dynprocs, session, time_pset, "mpi://WORLD", space_pset)

    !>  Create the pfasst structure
    call pf_pfasst_create_dynamic(pf, dynprocs, fname=pf_fname)
    pf%use_rk_stepper = .false.
    pf%use_sdc_sweeper = .true.

    !> Update space_color if dynamic
    !> Only relevant for debugging
    if (is_dynamic) then
        space_color = pf%rank
    end if

    print *, "============================================================================="
    print *, "space_pset = ", trim(time_pset), ", time_pset = ", trim(space_pset)
    print *, "space_color = ", space_color, ", time_color = ", time_color
    print *, "============================================================================="

    !> some more initializations
    spacial_coarsen_flag = 0
    call mpi_barrier(space_comm, ierr);
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi barrier fail, error=',ierr)
    call PfasstHypreInit(pf, lev_shape, space_comm, spacial_coarsen_flag)
    ! print *, "PfasstHypreInit done"
    ! print *,time_color,space_color,pf%rank

    !> echo error after each hook
    call pf_add_hook(pf, -1, PF_POST_BLOCK, echo_error)

    if (dump_values) then
        dump_dir = trim(dump_dir)
        print *, "Setting globalvars"
        ! Save some global variables for the dump hook
        call set_global_str("dump_dir", dump_dir)
        call set_global_int("time_color", time_color)
        call set_global_int("space_color", space_color)

        ! fails for spaces in dump dir
        call system('mkdir -p ' // trim(adjustl(dump_dir)))

        call pf_add_hook(pf, -1, PF_POST_BLOCK, dump_hook)
    end if


    !> setup hook to control how to resize libpfasst
    call pf_add_hook(pf, -1, PF_PRE_POT_RESIZE, resize_decider)

    !>  Output run parameters
    if (.not. is_dynamic .and. space_color == 0 .and. time_color == 0) then
        call print_loc_options(pf,un_opt=6)
    end if


    !> Set the initial condition
    level_index = pf%nlevels
    call hvf%create_single(y_0_base, level_index, lev_shape(pf%nlevels,:))
    call hvf%create_single(y_end_base, level_index, lev_shape(pf%nlevels,:))
    y_0 = cast_as_hypre_vector(y_0_base)
    y_end = cast_as_hypre_vector(y_end_base)

    if (.not. is_dynamic) then
        call initial(y_0)
        ! dump initial values
        if (dump_values) then
            write(fname, "(A,i5.5,A,i4.4,A,i4.4,A,i4.4,A)") &
                trim(adjustl(dump_dir)) // trim("/dump_step"), 0, &
                "_time", time_color, "_space", space_color, "_level", 2, ".csv"
            print *, "Dumping to ", trim(fname)
            call y_0%dump(fname)
        end if
        print *, "initial done"
    end if


    !> Do the PFASST time stepping
    call pf_pfasst_run(pf, y_0, dt, Tfin, nsteps, y_end, &
                       join_existing=is_dynamic, premature_exit=premature_exit)

    if (.not. premature_exit) then
       !if (pf%state%step .eq. nsteps-1) call y_end%eprint()

       !>  Wait for everyone to be done
       call mpi_barrier(pf%comm%comm, ierr)
    end if

    !> free time communicator
    call mpi_comm_free(pf%comm%comm, ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi comm free fail, error=',ierr)

    !> free space communicator
    call mpi_comm_free(space_comm, ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi comm free fail, error=',ierr)

    !>  Deallocate pfasst structure
    call pf_pfasst_destroy(pf)

    call pf_dynprocs_destroy(dynprocs)

  end subroutine run_pfasst
end program
