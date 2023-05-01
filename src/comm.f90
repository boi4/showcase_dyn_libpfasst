module pf_space_comm
  use pf_mod_mpi
  use pfasst
  implicit none
contains

   subroutine psetop2str(psetop, str)
     integer, intent(in) :: psetop
     character(len=*), intent(out) :: str

     character(len=23) :: myStrings(11) = [ &
       "MPI_PSETOP_NULL        ", &
       "MPI_PSETOP_ADD         ", &
       "MPI_PSETOP_SUB         ", &
       "MPI_PSETOP_REPLACE     ", &
       "MPI_PSETOP_MALLEABLE   ", &
       "MPI_PSETOP_GROW        ", &
       "MPI_PSETOP_SHRINK      ", &
       "MPI_PSETOP_UNION       ", &
       "MPI_PSETOP_DIFFERENCE  ", &
       "MPI_PSETOP_INTERSECTION", &
       "MPI_PSETOP_MULTI       " &
     ]

     str = mystrings(psetop+1)
   end subroutine psetop2str

  subroutine create_pset_grid(session, base_pset, base_comm, nspace, space_dim, &
                              space_comm, space_color, time_comm, time_color, time_pset)
    integer, intent(in) :: session
    character(len=*) , intent(in) :: base_pset
    integer, intent(in) :: base_comm
    integer, intent(in) :: nspace
    integer, intent(in) :: space_dim

    integer, intent(out) :: space_comm
    integer, intent(out) :: space_color
    integer, intent(out) :: time_comm
    integer, intent(out) :: time_color
    character(len=MPI_MAX_PSET_NAME_LEN) , intent(out) :: time_pset

    character(len=MPI_MAX_PSET_NAME_LEN)  :: output_psets(nspace)
    character(len=MPI_MAX_PSET_NAME_LEN)  :: input_psets(1)
    character(len=:), allocatable  :: splitstr
    character(len=20)  :: tmpstr

    real(pfdp) :: nspace_real
    integer :: base_size, base_rank
    integer :: noutput
    integer :: mgroup
    integer :: info
    integer :: ierr
    integer :: ierr2
    integer :: nlen
    integer :: i
    integer :: ntime
    integer :: op
    logical :: contains_key

    ! for now, nspace must be a perfect square
    if (space_dim .eq. 2) then
       nspace_real = sqrt(real(nspace))
       if (nspace_real .ne. nint(nspace_real)) then
          print'(a)', 'ERROR: create_simple_communicators: nspace must be perfect square.'
          stop
       end if
    end if

    call mpi_comm_size(base_comm, base_size, ierr)
    call mpi_comm_rank(base_comm, base_rank, ierr)

    ntime = base_size / nspace

    if (base_rank == 0) then
       ! prepare split arguments
       call mpi_info_create(info, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)
       allocate(character(20*nspace)::splitstr)

       splitstr = ""
       do i = 1, nspace
           write(tmpstr,'(I0)') ntime
           splitstr = trim(splitstr)//trim(tmpstr)
           if (i < nspace) splitstr = trim(splitstr)//","
       end do

       call mpi_info_set(info, "mpi_part_sizes", splitstr, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info set fail, error=',ierr)
       noutput = nspace
       input_psets(1) = base_pset
       op = MPI_PSETOP_SPLIT

       ! do the split
       print *, "Splitting ", trim(input_psets(1)), " into ", trim(splitstr)
       call MPI_Session_dyn_v2a_psetop(session, op, input_psets, 1, output_psets, noutput, info, ierr)

       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi psetop split operation failed, error=',ierr)
    end if

    ! determine which pset we are part of
    noutput = nspace
    do i = 1,noutput
        call mpi_bcast(output_psets(i), MPI_MAX_PSET_NAME_LEN, MPI_CHAR, 0, base_comm, ierr)
        if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierr)
    end do

    do i = 1,noutput
        call MPI_Session_get_pset_info(session, output_psets(i), info, ierr)
        if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get pset info fail, error=',ierr)
        call MPI_Info_get(info, "mpi_included", 20, tmpstr, contains_key, ierr)
        if (ierr /=0 .OR. .not. contains_key) call pf_stop(__FILE__,__LINE__,'mpi info get fail, error=',ierr)

        if (tmpstr == "True") then
            time_color = i - 1
            time_pset = output_psets(i)
            exit
          end if
    end do
    call mpi_barrier(base_comm, ierr)

    ! create communicator from time_pset
    call mpi_group_from_session_pset(session, time_pset, mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group from pset fail, error=',ierr)
    call mpi_barrier(base_comm, ierr)
    call mpi_comm_create_from_group(mgroup, "showcase", MPI_INFO_NULL, MPI_ERRORS_RETURN, time_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm create from group fail, error=',ierr)
    call mpi_group_free(mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierr)

    call mpi_barrier(time_comm, ierr)
    call mpi_barrier(base_comm, ierr)

    ! create space communicator based on rank in time communicator
    call MPI_Comm_rank(time_comm, space_color, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm rank fail, error=',ierr)

    print *, "Splitting base_comm into space_comm using color = ", space_color, " and key = ", time_color
    call MPI_Comm_split(base_comm, space_color, time_color, space_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm split fail, error=',ierr)

    print *, "time color = ", time_color, " space color = ", space_color

    call mpi_barrier(time_comm, ierr)
    call mpi_barrier(base_comm, ierr)
  end subroutine create_pset_grid

end module pf_space_comm
