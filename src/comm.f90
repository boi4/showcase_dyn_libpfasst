module pf_space_comm
  use mpi
  use pf_mod_mpi
  use pfasst
  implicit none
contains

  !> This function will creat a grid of process sets as seen below
  !>
  !>    example ntime: 4, nspace: 3
  !>                         time psets (time color 0,1,2,3 respectively)
  !>         ^                v v v v
  !>         |              +---------+
  !>         |           >  | 8 9 1011|
  !>         |     space    +---------+
  !>         |           >  | 4 5 6 7 |
  !>         |     psets    +---------+
  !>         |           >>>>>0 1 2 3 |
  !>         |           ^  +---------+
  !>         |           ^
  !>         |           ^
  !>         |           this one is both rank 0 in space pset and time pset and rank 0 in base_comm
  !>   time  |           this one will do the splits and unions
  subroutine create_pset_grid(session, base_pset, nspace, space_dim, &
                              space_comm, space_color, space_pset, &
                              time_comm, time_color, time_pset)
    integer, intent(in) :: session
    character(len=*) , intent(in) :: base_pset
    integer, intent(in) :: nspace
    integer, intent(in) :: space_dim

    integer, intent(out) :: space_comm
    integer, intent(out) :: space_color
    character(len=MPI_MAX_PSET_NAME_LEN) , intent(out) :: space_pset
    integer, intent(out) :: time_comm
    integer, intent(out) :: time_color
    character(len=MPI_MAX_PSET_NAME_LEN) , intent(out) :: time_pset

    character(len=MPI_MAX_PSET_NAME_LEN), allocatable  :: split_split_psets(:)
    character(len=MPI_MAX_PSET_NAME_LEN), allocatable  :: time_psets(:)
    character(len=MPI_MAX_PSET_NAME_LEN), allocatable  :: space_psets(:)
    character(len=MPI_MAX_PSET_NAME_LEN), allocatable  :: pset_buf(:)
    character(len=MPI_MAX_PSET_NAME_LEN), allocatable  :: input_psets(:)
    character(len=MPI_MAX_PSET_NAME_LEN)  :: single_pset
    character(len=:), allocatable  :: splitstr
    character(len=20)  :: tmpstr

    integer :: base_comm
    integer :: base_group
    integer :: base_size, base_rank
    integer :: noutput
    integer :: mgroup
    integer :: info
    integer :: ierr
    integer :: ierr2
    integer :: nlen
    integer :: i
    integer :: j
    integer :: ntime
    integer :: op
    integer :: status(MPI_STATUS_SIZE)
    logical :: contains_key
    logical :: contains_me


    ! create base_comm from base_pset
    call mpi_group_from_session_pset(session, base_pset, base_group, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group from pset fail, error=',ierr)
    call mpi_comm_create_from_group(base_group, base_pset, MPI_INFO_NULL, MPI_ERRORS_RETURN, base_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm create from group fail, error=',ierr)
    call mpi_group_free(base_group, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierr)

    call mpi_comm_size(base_comm, base_size, ierr)
    call mpi_comm_rank(base_comm, base_rank, ierr)

    ntime = base_size / nspace


    allocate(time_psets(nspace))
    allocate(space_psets(ntime))

    ! PART 1: split up base_comm into space psets
    ! ===================================================
    if (base_rank == 0) then
       ! allocate additional stuff
       allocate(character(20*nspace*ntime)::splitstr)
       allocate(split_split_psets(ntime*nspace))
       allocate(input_psets(ntime))

       print *, "Splitting base communicator into time psets"
       ! prepare split arguments
       call mpi_info_create(info, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)

       splitstr = ""
       do i = 1, ntime
           write(tmpstr,'(I0)') nspace
           splitstr = trim(splitstr)//trim(tmpstr)
           if (i < ntime) splitstr = trim(splitstr)//","
       end do

       call mpi_info_set(info, "mpi_part_sizes", splitstr, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info set fail, error=',ierr)
       noutput = ntime
       op = MPI_PSETOP_SPLIT

       ! do the split
       print *, "Splitting ", trim(base_pset), " into ", trim(splitstr)
       print *, "base_size is ", base_size
       call mpi_session_dyn_v2a_psetop(session, op, base_pset, 1, space_psets, noutput, info, ierr)
       call mpi_info_free(info, ierr2)
       if (ierr /=0 .or. ierr2 /= 0) call pf_stop(__FILE__,__LINE__,'mpi psetop split operation failed, error=',ierr,ierr2)


       ! PART 2: Split up space psets into split split psets
       ! ===================================================
       ! prepare split arguments once for all splits
       call mpi_info_create(info, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)

       splitstr = ""
       do i = 1, nspace
          write(tmpstr,'(I0)') 1
          splitstr = trim(splitstr)//trim(tmpstr)
          if (i < nspace) splitstr = trim(splitstr)//","
       end do

       call mpi_info_set(info, "mpi_part_sizes", splitstr, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info set fail, error=',ierr)

       do i = 1, ntime
          noutput = nspace
          op = MPI_PSETOP_SPLIT
          print *, "Splitting ", trim(space_psets(i)), " into ", trim(splitstr)
          call mpi_session_dyn_v2a_psetop(session, op, space_psets(i), 1, split_split_psets((i-1)*nspace+1:i*nspace), noutput, info, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi psetop split operation failed, error=',ierr)
       end do
       call mpi_info_free(info, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)


       ! PART 3: Union the split split psets to get the time psets
       ! ===================================================

       ! prepare arguments
       call mpi_info_create(info, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)

       do i = 1, nspace
          ! see figure in above comment to see why these ranks are in the same space pset
          do j = 1, ntime
             input_psets(j) = split_split_psets((j-1)*nspace + i)
          end do
          op = MPI_PSETOP_UNION

          ! prepare arguments
          ! do psetop
          print *, "Unioning "
          do j = 1, ntime
              print *, trim(input_psets(j))
          end do
          noutput = 1
          call mpi_session_dyn_v2a_psetop(session, op, input_psets, ntime, time_psets(i), noutput, info, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi psetop union operation failed, error=',ierr)
          print *, "-> Unioned to ", trim(time_psets(i))
       end do
       call mpi_info_free(info, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)

       ! deallocate additional stuff
       deallocate(input_psets)
       deallocate(splitstr)
       deallocate(split_split_psets)
    end if


    ! Part 4: Broadcast the space & time psets to all processes
    ! ===================================================
    call mpi_bcast(time_psets, nspace*MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, 0, base_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierr)

    call mpi_bcast(space_psets, ntime*MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, 0, base_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierr)

    ! Part 5: Figure out time & space pset of this process
    ! ===================================================
    ! check figure at the top of this routine to see how the psets are distributed
    space_color = base_rank / nspace
    time_color  = MOD(base_rank,nspace)


    space_pset = space_psets(space_color+1)
    time_pset = time_psets(time_color+1)

    ! print *, "nspace=", nspace
    ! print *, "my space_color=", space_color
    ! print *, "my space_pset=", trim(space_pset)
    ! print *, "ntime=", ntime
    ! print *, "my time_color=", time_color
    ! print *, "my time_pset=", trim(time_pset)

    ! call print_my_psets(session)

    !! double check we are actualy part of these psets
    !call pf_dynprocs_pset_contains_me(session, space_pset, contains_me)
    !if (.not. contains_me) then
    !   call pf_stop(__FILE__,__LINE__,"Error: space pset does not contain me")
    !end if

    !call pf_dynprocs_pset_contains_me(session, time_pset, contains_me)
    !if (.not. contains_me) then
    !   call pf_stop(__FILE__,__LINE__,"Error: time pset does not contain me")
    !end if

    deallocate(time_psets)
    deallocate(space_psets)


    ! Part 6: Derive communicators from the psets
    ! ===================================================
    call mpi_group_from_session_pset(session, space_pset, mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group from pset fail, error=',ierr)
    call mpi_comm_create_from_group(mgroup, space_pset, MPI_INFO_NULL, MPI_ERRORS_RETURN, space_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm create from group fail, error=',ierr)
    call mpi_group_free(mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierr)


    call mpi_group_from_session_pset(session, time_pset, mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group from pset fail, error=',ierr)
    call mpi_comm_create_from_group(mgroup, time_pset, MPI_INFO_NULL, MPI_ERRORS_RETURN, time_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm create from group fail, error=',ierr)
    call mpi_group_free(mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierr)


    call mpi_barrier(time_comm, ierr)
    call mpi_barrier(space_comm, ierr)
    call mpi_barrier(base_comm, ierr)
  end subroutine create_pset_grid







  ! Debugging

  subroutine print_my_psets(session)
    integer, intent(in) :: session
    character(len=MPI_MAX_PSET_NAME_LEN) pset
    integer :: n
    integer :: len
    integer :: ierr
    integer :: i
    integer :: info
    character(len=20)  :: tmpstr
    logical :: contains_key

    call mpi_session_get_num_psets(session, MPI_INFO_NULL, n, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get num psets fail, error=',ierr)

    do i = 0, n-1
      len = MPI_MAX_PSET_NAME_LEN
      call mpi_session_get_nth_pset(session, MPI_INFO_NULL, i, len, pset, ierr)
      if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get pset name fail, error=',ierr)

      call mpi_session_get_pset_info(session, pset, info, ierr)
      if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get pset info fail, error=',ierr)
      call mpi_info_get(info, "mpi_included", 20, tmpstr, contains_key, ierr)
      if (ierr /=0 .OR. .not. contains_key) call pf_stop(__FILE__,__LINE__,'mpi info get fail, error=',ierr)

      if (tmpstr == "True") then
         print *, "pset ", i, ": ", trim(pset)
      else
         print *, "//pset ", i, ": ", trim(pset)
      end if
    end do
  end subroutine print_my_psets

 subroutine get_pset_size(session, global_comm, pset, size)
   integer, intent(in) :: session
   integer, intent(in) :: global_comm
   character(len=MPI_MAX_PSET_NAME_LEN), intent(in) :: pset
   integer, intent(out) :: size

   logical :: contains_me
   integer :: union_size
   integer :: indicator
   integer :: ierr
   call pf_dynprocs_pset_contains_me(session, pset, contains_me)
   if (contains_me) then
       indicator = 1
   else
       indicator = 0
   end if

   call mpi_allreduce(indicator, size, 1, MPI_INTEGER, MPI_SUM, global_comm, ierr)
   if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi allreduce fail, error=',ierr)
  end subroutine get_pset_size


end module pf_space_comm
