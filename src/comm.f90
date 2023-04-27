module pf_space_comm
  use pf_mod_mpi
  use pfasst
  implicit none
contains

   subroutine create_simple_communicators(base_comm, nspace, ntime, space_comm, time_comm, space_color, time_color, space_dim)
    integer, intent(in) :: base_comm
    integer, intent(in) :: nspace, ntime
    integer, intent(in) :: space_dim
    integer, intent(out) :: space_comm, time_comm
    integer, intent(out) :: space_color, time_color

    integer :: nproc, rank, error
    real(pfdp) :: nspace_real

    ! check size
    call mpi_comm_size(base_comm, nproc, error)
    call mpi_comm_rank(base_comm, rank,  error)

    if (nproc .ne. (nspace * ntime)) then
       print '(a)', 'ERROR: create_simple_communicators: processor number mismatch.'
       print '(a,i4,a,i4)', '       Expecting ', &
            nspace * ntime, ' MPI processors but received ', &
            nproc
       stop
    end if

    ! for now, nspace must be a perfect square
    if (space_dim .eq. 2) then
       nspace_real = sqrt(real(nspace))
       if (nspace_real .ne. nint(nspace_real)) then
          print'(a)', 'ERROR: create_simple_communicators: nspace must be perfect square.'
          stop
       end if
    end if

    if (ntime == 1) then
       time_color = rank
       space_color = 0
       space_comm = base_comm
       time_comm  = MPI_COMM_SELF
    else if (nspace == 1) then
       time_color = 0
       space_color = rank
       space_comm = MPI_COMM_SELF
       time_comm  = base_comm
    else
       ! split by color
       space_color = rank / nspace
       call mpi_comm_split(base_comm, space_color, rank, space_comm, error)

       time_color = mod(rank, nspace)
       call mpi_comm_split(base_comm, time_color, rank, time_comm, error)
    end if

  end subroutine create_simple_communicators

end module pf_space_comm
