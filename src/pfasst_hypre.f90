module pfasst_hypre
  use encap
  use pfasst
  use pf_my_level
  use probin
  use pf_my_sweeper
  use pf_mod_mgrit
  implicit none
contains

  subroutine PfasstHypreInit(pf, lev_shape, space_comm, spacial_coarsen_flag)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, allocatable, intent(inout) :: lev_shape(:,:)
    integer, intent(in) :: space_comm
    integer, intent(in) :: spacial_coarsen_flag
    
    type(pf_comm_t) :: comm
    type(my_sweeper_t) :: sw_finest, sw_lev
    type(my_level_t) :: my_lev

    integer :: l, l_finest   !  loop variable over levels
    integer :: n, m
    integer :: level_index

    real(pfdp) :: f
    integer :: nrows, ilower0, ilower1, iupper0, iupper1
    integer :: n_init, refine_factor
    logical :: setup_start_coarse_flag



    allocate(lev_shape(pf%nlevels,10))
    !> Loop over levels and set some level specific parameters
    do l = 1,pf%nlevels
       lev_shape(l,1) = num_grid_points
       lev_shape(l,2) = space_comm
       lev_shape(l,3) = space_dim
       lev_shape(l,4) = max_space_v_cycles
       lev_shape(l,10) = spacial_coarsen_flag
       !>  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)

       !>  Allocate the user specific data constructor
       allocate(hypre_vector_factory::pf%levels(l)%ulevel%factory)

       !>  Add the sweeper to the level
       allocate(my_sweeper_t::pf%levels(l)%ulevel%sweeper)

       call pf_level_set_size(pf, l, lev_shape(l,:), 0)
    end do

    l_finest = pf%nlevels
    if (spacial_coarsen_flag .eq. 1) then
       sw_finest = cast_as_my_sweeper_t(pf%levels(l_finest)%ulevel%sweeper)
       call HypreSolverInit(sw_finest%c_hypre_solver_ptr, &
                            l_finest, &
                            num_grid_points, &
                            space_comm, &
                            space_dim, &
                            max_space_v_cycles, &
                            pf%nlevels, &
                            spacial_coarsen_flag)
    end if

    do l = pf%nlevels,1,-1
       if (spacial_coarsen_flag .eq. 1) then
          sw_lev = cast_as_my_sweeper_t(pf%levels(l_finest)%ulevel%sweeper)
       else
          sw_lev = cast_as_my_sweeper_t(pf%levels(l)%ulevel%sweeper)
          call HypreSolverInit(sw_lev%c_hypre_solver_ptr, &
                               l, &
                               num_grid_points, &
                               space_comm, &
                               space_dim, &
                               max_space_v_cycles, &
                               pf%nlevels, &
                               spacial_coarsen_flag)
       end if

       n = sw_lev%get_nrows(l)
       ilower0 = sw_lev%get_extent(l, 0)
       ilower1 = sw_lev%get_extent(l, 1)
       iupper0 = sw_lev%get_extent(l, 2)
       iupper1 = sw_lev%get_extent(l, 3)

       lev_shape(l,5) = n
       lev_shape(l,6) = ilower0
       lev_shape(l,7) = ilower1
       lev_shape(l,8) = iupper0
       lev_shape(l,9) = iupper1

       call pf_level_set_size(pf, l, lev_shape(l,:), n)
    end do

    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)

    call sweeper_hypre_set_level_data(pf)
  end subroutine PfasstHypreInit

end module pfasst_hypre
