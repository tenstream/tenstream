@test
subroutine test_subtract_numbers()

   use hello
   use pfunit_mod

   implicit none

   real(8) :: res

   call subtract_numbers(1.0d0, 2.0d0, res)
   @assertEqual(res, -1.0d0)

end subroutine
