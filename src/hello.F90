module hello

   implicit none

   public add_numbers
   public subtract_numbers
   public multiply_numbers

   private

contains

   subroutine add_numbers(f1, f2, res)

      real(8), intent(in)  :: f1
      real(8), intent(in)  :: f2
      real(8), intent(out) :: res

      res = f1 + f2

   end subroutine

   subroutine subtract_numbers(f1, f2, res)

      real(8), intent(in)  :: f1
      real(8), intent(in)  :: f2
      real(8), intent(out) :: res

      res = f1 - f2

   end subroutine

   subroutine multiply_numbers(f1, f2, res)

      real(8), intent(in)  :: f1
      real(8), intent(in)  :: f2
      real(8), intent(out) :: res

      res = f1*f2

   end subroutine

end module
