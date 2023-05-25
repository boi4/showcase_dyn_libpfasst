module pf_space_comm
  use mpi
  use pf_mod_mpi
  use pfasst
  implicit none
contains

   ! subroutine print_my_psets(session)
   !   integer, intent(in) :: session
   !   character(len=MPI_MAX_PSET_NAME_LEN) pset
   !   integer :: n
   !   integer :: len
   !   integer :: ierr
   !   integer :: i
   !   integer :: info
   !   character(len=20)  :: tmpstr
   !   logical :: contains_key

   !   call mpi_session_get_num_psets(session, MPI_INFO_NULL, n, ierr)
   !   if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get num psets fail, error=',ierr)

   !   do i = 0, n-1
   !     len = MPI_MAX_PSET_NAME_LEN
   !     call mpi_session_get_nth_pset(session, MPI_INFO_NULL, i, len, pset, ierr)
   !     if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get pset name fail, error=',ierr)

   !     call mpi_session_get_pset_info(session, pset, info, ierr)
   !     if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get pset info fail, error=',ierr)
   !     call mpi_info_get(info, "mpi_included", 20, tmpstr, contains_key, ierr)
   !     if (ierr /=0 .OR. .not. contains_key) call pf_stop(__FILE__,__LINE__,'mpi info get fail, error=',ierr)

   !     if (tmpstr == "True") then
   !        print *, "pset ", i, ": ", trim(pset)
   !     end if
   !   end do
   ! end subroutine print_my_psets


  ! TODO: add comment
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

    character(len=MPI_MAX_PSET_NAME_LEN), allocatable  :: output_psets(:)
    character(len=MPI_MAX_PSET_NAME_LEN), allocatable  :: input_psets(:)
    character(len=MPI_MAX_PSET_NAME_LEN), allocatable  :: pset_buf(:)
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

    ! allocate psets
    allocate(input_psets(1))
    allocate(output_psets(ntime))
    allocate(character(20*ntime)::splitstr)


    ! PART 1: split up base_comm into space psets
    ! ===================================================
    if (base_rank == 0) then
       print *, "Splitting base communicator into space psets"
       ! prepare split arguments
       call mpi_info_create(info, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)

       splitstr = ""
       do i = 1, ntime
           write(tmpstr,'(I0)') nspace
           splitstr = trim(splitstr)//trim(tmpstr)
           if (i < nspace) splitstr = trim(splitstr)//","
       end do

       call mpi_info_set(info, "mpi_part_sizes", splitstr, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info set fail, error=',ierr)
       noutput = ntime
       input_psets(1) = base_pset
       op = MPI_PSETOP_SPLIT

       ! do the split
       print *, "Splitting ", trim(input_psets(1)), " into ", trim(splitstr)
       call mpi_session_dyn_v2a_psetop(session, op, input_psets, 1, output_psets, noutput, info, ierr)
       call mpi_info_free(info, ierr2)
       if (ierr /=0 .or. ierr2 /= 0) call pf_stop(__FILE__,__LINE__,'mpi psetop split operation failed, error=',ierr,ierr2)
    end if

    ! determine which pset we are part of
    noutput = ntime
    call mpi_bcast(output_psets, MPI_MAX_PSET_NAME_LEN*noutput, MPI_CHAR, 0, base_comm, ierr)

    do i = 1,noutput
        call mpi_session_get_pset_info(session, output_psets(i), info, ierr)
        if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get pset info fail, error=',ierr)
        call mpi_info_get(info, "mpi_included", 20, tmpstr, contains_key, ierr)
        if (ierr /=0 .OR. .not. contains_key) call pf_stop(__FILE__,__LINE__,'mpi info get fail, error=',ierr)

        if (tmpstr == "True") then
           space_color = i - 1
           space_pset = output_psets(i)
           exit
        end if
    end do
    call mpi_barrier(base_comm, ierr)

    ! can free output_psets
    deallocate(input_psets)
    deallocate(output_psets)
    deallocate(splitstr)

    ! create communicator from space_pset
    call mpi_group_from_session_pset(session, space_pset, mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group from pset fail, error=',ierr)
    call mpi_comm_create_from_group(mgroup, space_pset, MPI_INFO_NULL, MPI_ERRORS_RETURN, space_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm create from group fail, error=',ierr)
    call mpi_group_free(mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierr)


    ! create space pset based on rank in time communicator
    call mpi_comm_rank(space_comm, time_color, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm rank fail, error=',ierr)

    ! sync up
    call mpi_barrier(space_comm, ierr)
    call mpi_barrier(base_comm, ierr)



    ! PART 2: create time psets
    ! ===================================================

    ! we need a temporary communicator so that all space rank 0's can communicate
    call mpi_comm_split(base_comm, time_color, space_color, time_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm split fail, error=',ierr)

    !
    !    example ntime: 4, nspace: 3
    !                         time psets (time color 0,1,2,3 respectively)
    !         ^                v v v v
    !         |              +---------+
    !         |           >  | 8 9 1011|
    !         |     space    +---------+
    !         |           >  | 4 5 6 7 |
    !         |     psets    +---------+
    !         |           >>>>>0 1 2 3 |
    !         |           ^  +---------+
    !         |           ^    ^ these have rank 0 in time psets and do the split in their respective space pset
    !         |           ^
    !         |           this one is both rank 0 in space pset and time pset and rank 0 in base_comm
    !   time  |           this one will do the unions
    !
    !    example ntime: 3, nspace: 4

    if (time_color == 0) then
       ! the "local leader" of each space pset does the split
       allocate(input_psets(1))
       allocate(output_psets(nspace))
       allocate(character(20*nspace)::splitstr)
       print *, "Splitting up my space pset"

       ! prepare split arguments
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
       noutput = nspace
       input_psets(1) = space_pset
       op = MPI_PSETOP_SPLIT

       ! do the split
       print *, "Splitting ", trim(input_psets(1)), " into ", trim(splitstr)
       call mpi_session_dyn_v2a_psetop(session, op, input_psets, 1, output_psets, noutput, info, ierr)
       call mpi_info_free(info, ierr2)
       if (ierr /=0 .or. ierr2 /= 0) call pf_stop(__FILE__,__LINE__,'mpi psetop split operation failed, error=',ierr,ierr2)


       ! prepare for gather
       allocate(pset_buf(nspace * ntime))

       print *, "Gathering all splits to base comm"
       ! gather all single-element splits to global comm
       call mpi_gather(output_psets, nspace*MPI_MAX_PSET_NAME_LEN, MPI_CHAR, &
                       pset_buf, nspace*MPI_MAX_PSET_NAME_LEN, MPI_CHAR, &
                       0, time_comm, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi gather fail, error=',ierr)
       print *, "Done gathering splits to base comm"

       deallocate(splitstr)
       deallocate(input_psets)
       deallocate(output_psets)

       if (space_color == 0) then
          allocate(input_psets(ntime))
          allocate(output_psets(nspace))

          print *, "Doing unions for each time comm"

          ! do the unions
          call mpi_info_create(info, ierr)
          do i = 1, nspace
             ! see figure in above comment to see why these ranks are in the same time pset
             do j = 1, ntime
                input_psets(j) = pset_buf((j-1)*nspace + i)
             end do
             op = MPI_PSETOP_UNION

             ! prepare arguments
             ! do psetop
             print *, "Unioning "
             do j = 1, ntime
                 print *, trim(input_psets(j))
             end do
             call mpi_session_dyn_v2a_psetop(session, op, input_psets, ntime, output_psets(i), noutput, info, ierr)
             if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi psetop union operation failed, error=',ierr)
             print *, "-> Unioned to ", trim(output_psets(i))
          end do
          call mpi_info_free(info, ierr)

          deallocate(input_psets)
       end if
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)

       deallocate(pset_buf)
    end if


    ! broadcast unions to all ranks in base comm
    if (base_rank /= 0) then
       allocate(output_psets(nspace))
    end if
    print *, "broadcasting unions to all ranks in base comm"
    call mpi_bcast(output_psets, nspace*MPI_MAX_PSET_NAME_LEN, MPI_CHAR, 0, base_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierr)

    ! get time pset
    time_pset = output_psets(time_color + 1)
    deallocate(output_psets)
    print *, "My time pset is ", trim(time_pset)


    ! assert we are part of time pset
    call mpi_session_get_pset_info(session, time_pset, info, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get pset info fail, error=',ierr)
    call mpi_info_get(info, "mpi_included", 20, tmpstr, contains_key, ierr)
    if (ierr /=0 .OR. .not. contains_key) call pf_stop(__FILE__,__LINE__,'mpi info get fail, error=',ierr)

    ! call print_my_psets(session)

    if (tmpstr /= "True") then
       call pf_stop(__FILE__,__LINE__,"ERROR: I am not part of time pset")
    end if

    ! re-derive time comm
    call mpi_comm_free(time_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm free fail, error=',ierr)

    call mpi_group_from_session_pset(session, time_pset, mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group from pset fail, error=',ierr)
    call mpi_comm_create_from_group(mgroup, time_pset, MPI_INFO_NULL, MPI_ERRORS_RETURN, time_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm create from group fail, error=',ierr)
    call mpi_group_free(mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierr)

    call mpi_barrier(space_comm, ierr)
    call mpi_barrier(time_comm, ierr)
    call mpi_barrier(base_comm, ierr)
  end subroutine create_pset_grid











  ! old version of the above that first splits into time psets, then space psets
  subroutine create_pset_grid_time_major(session, base_pset, nspace, space_dim, &
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

    character(len=MPI_MAX_PSET_NAME_LEN), allocatable  :: output_psets(:)
    character(len=MPI_MAX_PSET_NAME_LEN), allocatable  :: input_psets(:)
    character(len=MPI_MAX_PSET_NAME_LEN), allocatable  :: pset_buf(:)
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

    ! allocate psets
    allocate(input_psets(1))
    allocate(output_psets(nspace))
    allocate(character(20*nspace)::splitstr)


    ! PART 1: split up base_comm into time psets
    ! ===================================================
    if (base_rank == 0) then
       print *, "Splitting base communicator into time psets"
       ! prepare split arguments
       call mpi_info_create(info, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)

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
       call mpi_session_dyn_v2a_psetop(session, op, input_psets, 1, output_psets, noutput, info, ierr)
       call mpi_info_free(info, ierr2)
       if (ierr /=0 .or. ierr2 /= 0) call pf_stop(__FILE__,__LINE__,'mpi psetop split operation failed, error=',ierr,ierr2)
    end if

    ! determine which pset we are part of
    noutput = nspace
    call mpi_bcast(output_psets, MPI_MAX_PSET_NAME_LEN*noutput, MPI_CHAR, 0, base_comm, ierr)

    do i = 1,noutput
        call mpi_session_get_pset_info(session, output_psets(i), info, ierr)
        if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get pset info fail, error=',ierr)
        call mpi_info_get(info, "mpi_included", 20, tmpstr, contains_key, ierr)
        if (ierr /=0 .OR. .not. contains_key) call pf_stop(__FILE__,__LINE__,'mpi info get fail, error=',ierr)

        if (tmpstr == "True") then
           time_color = i - 1
           time_pset = output_psets(i)
           exit
        end if
    end do
    call mpi_barrier(base_comm, ierr)

    ! can free output_psets
    deallocate(input_psets)
    deallocate(output_psets)
    deallocate(splitstr)

    ! create communicator from time_pset
    call mpi_group_from_session_pset(session, time_pset, mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group from pset fail, error=',ierr)
    call mpi_comm_create_from_group(mgroup, time_pset, MPI_INFO_NULL, MPI_ERRORS_RETURN, time_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm create from group fail, error=',ierr)
    call mpi_group_free(mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierr)


    ! create space pset based on rank in time communicator
    call mpi_comm_rank(time_comm, space_color, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm rank fail, error=',ierr)

    ! sync up
    call mpi_barrier(time_comm, ierr)
    call mpi_barrier(base_comm, ierr)



    ! PART 2: create space psets
    ! ===================================================

    ! we need a temporary communicator so that all time rank 0's can communicate
    call mpi_comm_split(base_comm, space_color, time_color, space_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm split fail, error=',ierr)

    ! only "time layer 0" is communicating in this if:
    !                         time psets (time color 0,1,2,3 respectively) (each running pfasst)
    !         ^                v v v v
    !         |              +---------+
    !         |           >  | 2 5 8 11|
    !         |     space    +---------+
    !         |           >  | 1 4 7 10|
    !         |     psets    +---------+
    !         |           >  | 0 3 6 9 | < these have rank 0 in time psets and do the split in their respective time pset
    !         |              +---------+
    !         |                ^
    !         |                this one is both rank 0 in time pset and space pset and rank 0 in base_comm
    !    time |                this one will do the unions
    !
    !    example nspace: 4, ntime: 3

    print *, time_color, space_color
    if (space_color == 0) then
       ! the "local leader" of each time pset does the split
       allocate(input_psets(1))
       allocate(output_psets(ntime))
       allocate(character(20*ntime)::splitstr)
       print *, "Splitting up my time pset"

       ! prepare split arguments
       call mpi_info_create(info, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)

       splitstr = ""
       do i = 1, ntime
          write(tmpstr,'(I0)') 1
          splitstr = trim(splitstr)//trim(tmpstr)
          if (i < ntime) splitstr = trim(splitstr)//","
       end do

       call mpi_info_set(info, "mpi_part_sizes", splitstr, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info set fail, error=',ierr)
       noutput = ntime
       input_psets(1) = time_pset
       op = MPI_PSETOP_SPLIT

       ! do the split
       print *, "Splitting ", trim(input_psets(1)), " into ", trim(splitstr)
       call mpi_session_dyn_v2a_psetop(session, op, input_psets, 1, output_psets, noutput, info, ierr)
       call mpi_info_free(info, ierr2)
       if (ierr /=0 .or. ierr2 /= 0) call pf_stop(__FILE__,__LINE__,'mpi psetop split operation failed, error=',ierr,ierr2)


       ! prepare for gather
       allocate(pset_buf(ntime * nspace))

       print *, "Gathering all splits to base comm"
       ! gather all single-element splits to global comm
       call mpi_gather(output_psets, ntime*MPI_MAX_PSET_NAME_LEN, MPI_CHAR, &
                       pset_buf, ntime*MPI_MAX_PSET_NAME_LEN, MPI_CHAR, &
                       0, space_comm, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi gather fail, error=',ierr)
       print *, "Done gathering splits to base comm"

       deallocate(splitstr)
       deallocate(input_psets)
       deallocate(output_psets)

       if (time_color == 0) then
          allocate(input_psets(nspace))
          allocate(output_psets(ntime))

          print *, "Doing unions for each space comm"

          ! do the unions
          call mpi_info_create(info, ierr)
          do i = 1, ntime
             ! see figure in above comment to see why these ranks are in the same space pset
             do j = 1, nspace
                input_psets(j) = pset_buf((j-1)*ntime + i)
             end do
             op = MPI_PSETOP_UNION

             ! prepare arguments
             ! do psetop
             print *, "Unioning "
             do j = 1, nspace
                 print *, trim(input_psets(j))
             end do
             call mpi_session_dyn_v2a_psetop(session, op, input_psets, nspace, output_psets(i), noutput, info, ierr)
             if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi psetop union operation failed, error=',ierr)
             print *, "-> Unioned to ", trim(output_psets(i))
          end do
          call mpi_info_free(info, ierr)

          deallocate(input_psets)
       end if
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)

       deallocate(pset_buf)
    end if


    ! broadcast unions to all ranks in base comm
    if (base_rank /= 0) then
       allocate(output_psets(ntime))
    end if
    print *, "broadcasting unions to all ranks in base comm"
    call mpi_bcast(output_psets, ntime*MPI_MAX_PSET_NAME_LEN, MPI_CHAR, 0, base_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierr)

    ! get space pset
    space_pset = output_psets(space_color + 1)
    deallocate(output_psets)
    print *, "My space pset is ", trim(space_pset)


    ! assert we are part of space pset
    call mpi_session_get_pset_info(session, space_pset, info, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get pset info fail, error=',ierr)
    call mpi_info_get(info, "mpi_included", 20, tmpstr, contains_key, ierr)
    if (ierr /=0 .OR. .not. contains_key) call pf_stop(__FILE__,__LINE__,'mpi info get fail, error=',ierr)

    ! call print_my_psets(session)

    if (tmpstr /= "True") then
       call pf_stop(__FILE__,__LINE__,"ERROR: I am not part of space pset")
    end if

    ! re-derive space comm
    call mpi_comm_free(space_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm free fail, error=',ierr)

    call mpi_group_from_session_pset(session, space_pset, mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group from pset fail, error=',ierr)
    call mpi_comm_create_from_group(mgroup, space_pset, MPI_INFO_NULL, MPI_ERRORS_RETURN, space_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm create from group fail, error=',ierr)
    call mpi_group_free(mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierr)

    call mpi_barrier(time_comm, ierr)
    call mpi_barrier(space_comm, ierr)
    call mpi_barrier(base_comm, ierr)
  end subroutine create_pset_grid_time_major

end module pf_space_comm
