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

  integer ::  ierror

  !> Initialize MPI
  call mpi_init(ierror)
  if (ierror /= 0) &
       stop "ERROR: Can't initialize MPI."

  !> Call the  solver 
  call run_pfasst()

  !> Close mpi
  call mpi_finalize(ierror)

contains
  !>  This subroutine set ups and calls libpfasst 
  subroutine run_pfasst()  
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
    real(pfdp) :: f
    integer :: nrows, ilower0, ilower1, iupper0, iupper1
    integer :: spacial_coarsen_flag
    type(mgrit_level_data), allocatable :: mg_ld(:)
    character(len=1000) :: fname

    integer :: tmp
    character(len=100) :: tmpstr

    ! check size
    call mpi_comm_size(MPI_COMM_WORLD, nproc, error)
    call mpi_comm_rank(MPI_COMM_WORLD, rank,  error)

    !> Read problem parameters
    call probin_init(pf_fname)

    n = num_grid_points * num_grid_points
    call create_simple_communicators(nspace, ntime, space_comm, time_comm, space_color, time_color, space_dim)

    !>  Set up communicator
    call pf_mpi_create(comm, time_comm)
    

    if ((solver_type .eq. 1) .or. (solver_type .eq. 2) .or. (solver_type .eq. 3)) then
       pf%use_rk_stepper = .true.
       pf%use_sdc_sweeper = .false.
    else
       pf%use_rk_stepper = .false.
       pf%use_sdc_sweeper = .true.
    end if

    !>  Create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=pf_fname)
    if ((solver_type .eq. 2) .and. (pf%nlevels .ne. 2)) then
       print *, 'ERROR: nlevels must be 2 for Parareal.'
       return
    end if

    spacial_coarsen_flag = 0
    call PfasstHypreInit(pf, mg_ld, lev_shape, space_color, time_color, spacial_coarsen_flag)
    print *, "PfasstHypreInit done"
    print *,time_color,space_color,pf%rank
    
    ! !>  Add some hooks for output
    ! if (solver_type .eq. 1) then
    !    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)
    ! else
    !    !call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)
    !    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)
    ! end if
    ! call pf_add_hook(pf, -1, PF_POST_BLOCK, print_sol)
    ! call pf_add_hook(pf, -1, PF_POST_PREDICTOR, ppred)
    ! call pf_add_hook(pf, -1, PF_POST_BLOCK, echo_error)

    if (dump_values) then
        ! Save some global variables for the dump hook
        call set_global_str("dump_dir", dump_dir)
        call set_global_int("time_color", time_color)
        call set_global_int("space_color", space_color)

        ! TODO: what about spaces in dump_dir?
        call system('mkdir -p ' // trim(adjustl(dump_dir)))

        call pf_add_hook(pf, -1, PF_POST_BLOCK, dump_hook)
    end if

    !>  Output run parameters
    call print_loc_options(pf,un_opt=6)

    level_index = pf%nlevels
    !> Set the initial condition
    call hvf%create_single(y_0_base, level_index, lev_shape(pf%nlevels,:))
    call hvf%create_single(y_end_base, level_index, lev_shape(pf%nlevels,:))
    y_0 = cast_as_hypre_vector(y_0_base)
    y_end = cast_as_hypre_vector(y_end_base)
    call initial(y_0)

    ! dump initial values
    if (dump_values) then
        write(fname, "(A,A,i5.5,A,i4.4,A,i4.4,A)") trim(adjustl(dump_dir)), "/dump_step", 0, "_time", time_color, "_space", space_color, ".csv"

        write(fname, "(A,A,i5.5,A,i4.4,A,i4.4,A,i1.1,A)") &
            trim(adjustl(dump_dir)), "/dump_step", 0, &
            "_time", time_color, "_space", space_color, "_level", 1, ".csv"
        call y_0%dump(fname)
    end if

    print *, "initial done"


    !> Do the PFASST time stepping
    if (solver_type .eq. 1) then
       call pf_MGRIT_run(pf, mg_ld, y_0, y_end)
    else if (solver_type .eq. 2) then
       call pf_parareal_run(pf, y_0, dt, Tfin, nsteps, y_end)
    else if (solver_type .eq. 3) then
       call initialize_results(pf)
       if (pf%save_timings > 0) call pf_start_timer(pf, T_TOTAL)
       call pf%levels(1)%ulevel%stepper%do_n_steps(pf, 1, T0, y_0, y_end, dt, nsteps)
       if (pf%save_timings > 0) call pf_stop_timer(pf, T_TOTAL)
       call pf_dump_stats(pf)
    else
       call pf_pfasst_run(pf, y_0, dt, Tfin, nsteps, y_end)
    end if
    ! if (pf%rank .eq. pf%comm%nproc-1) call y_end%eprint()
    ! call y_end%eprint()

!    if (space_color == ntime - 1) then
!       ! coarse grid solution
!       !write(filename, "(A,i1,A)") "data/final_cdump_", time_color, ".csv"
!       !y_other = cast_as_hypre_vector(pf%levels(1)%qend)
!       !call y_other%dump(filename)
!
!       write(filename, "(A,i1,A)") "data/final_dump_", time_color, ".csv"
!       call y_end%dump(filename)
!    end if

    !>  Wait for everyone to be done
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    
    !>  Deallocate pfasst structure
    call pf_pfasst_destroy(pf)

  end subroutine run_pfasst
end program
