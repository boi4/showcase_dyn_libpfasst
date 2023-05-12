module global_state
   use, intrinsic :: iso_c_binding

   implicit none

   interface
      subroutine get_global_int_c(key, res) bind(c, name="get_global_int_c")
         use, intrinsic :: iso_c_binding, only: c_int, c_char
         implicit none
         character(c_char), dimension(*), intent(in) :: key
         integer(c_int), intent(out) :: res
      end subroutine get_global_int_c

      subroutine set_global_int_c(key, val) bind(c, name="set_global_int_c")
         use, intrinsic :: iso_c_binding, only: c_int, c_char
         character(c_char), dimension(*), intent(in) :: key
         integer(c_int), intent(in) :: val
      end subroutine set_global_int_c

      subroutine get_global_str_c(key, res) bind(c, name="get_global_str_c")
         use, intrinsic :: iso_c_binding, only: c_int, c_char
         character(c_char), dimension(*), intent(in) :: key
         character(c_char), dimension(*), intent(out) :: res
      end subroutine get_global_str_c

      subroutine set_global_str_c(key, val) bind(c, name="set_global_str_c")
         use, intrinsic :: iso_c_binding, only: c_int, c_char
         character(c_char), dimension(*), intent(in) :: key
         character(c_char), dimension(*), intent(in) :: val
      end subroutine set_global_str_c
   end interface

contains

   subroutine get_global_int(key, res)
      character(len=*), intent(in) :: key
      integer, intent(out) :: res

      character(c_char), dimension(len(key)+1) :: key_c

      ! Convert key to key_c
      key_c = transfer(trim(key)//char(0), key_c)

      call get_global_int_c(key_c, res)
   end subroutine get_global_int


   subroutine set_global_int(key, val)
      character(len=*), intent(in) :: key
      character(c_char), dimension(len(key)+1) :: key_c
      integer, intent(in) :: val

      ! Convert key to key_c
      key_c = transfer(trim(key)//char(0), key_c)

      call set_global_int_c(key_c, val)
   end subroutine set_global_int


   subroutine get_global_str(key, res)
      character(len=*), intent(in) :: key
      character(len=*), intent(out) :: res

      character(c_char), dimension(len(key)+1) :: key_c
      character(c_char), dimension(len(res)+1) :: res_c
      integer :: pos,i

      ! Convert key to key_c
      key_c = transfer(trim(key)//char(0), key_c)

      call get_global_str_c(key_c, res_c)

      res = transfer(res_c, res)

      pos = index(res, char(0))

      ! Pad res with spaces after and including the first null byte
      if (pos > 0) then
          do i = pos, len(res)
              res(i:i) = ' '
          end do
      endif

   end subroutine get_global_str


   subroutine set_global_str(key, val)
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: val

      character(c_char), dimension(len(key)+1) :: key_c
      character(c_char), dimension(len(val)+1) :: val_c

      ! Convert key to key_c
      key_c = transfer(trim(key)//char(0), key_c)
      ! Convert val to val_c
      val_c = transfer(trim(val)//char(0), val_c)

      call set_global_str_c(key_c, val_c)
   end subroutine set_global_str

end module global_state
