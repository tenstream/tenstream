program main
      use m_eddington
      use mpi
      use m_data_parameters, only : mpiint,ireals,iintegers,one,zero, init_mpi_data_parameters,i1

      real(ireals) :: inp(4) ! tau, omega0, g, mu0
      real(ireals) :: out(7) ! a11,a12,a13,a23,a33,g1,g2

      call MPI_Init(ierr)
      call init_mpi_data_parameters(MPI_COMM_WORLD)

      print *,'Checking eddington coefficients'

      inp = [ 0.4824516550E-01, 0.5542391539, 0.4550637007, 1.00000000000]
      call calc(inp)

      inp = [0.2018013448, 0.3797843754, 0.4556422830, 1.00000000000]
      call calc(inp)

      inp = [ 0.1731484532,  0.6180083156, 0.4121485054, 1.00000000000 ]
      call calc(inp)

      inp = [ 0.1931012775E-04, 0.4384377003, -0.0000000000E+00,  0.7070999742 ]
      call calc(inp)

      inp = [ 4.895462513 , 0.3104626103E-05 , -0.0000000000E+00 , 0.49999997019767761  ] ! resonance case
      call calc(inp)

      inp = [ 3.2662250689525390E-011 , 0.99999171495417127 ,       0.0000000000000000  , 0.17364817766693041 ]
      call calc(inp)

      inp = [ 2.9317851124478626E-012,   1.0000000000000000  ,      0.0000000000000000 , 0.17364817766693041 ]
      call calc(inp)

      inp = [1.93321303E-10,  0.999984443 ,      2.22044605E-16,  0.17364817766693041 ]
      call calc(inp)

      inp = [7.89528581E-11,  0.999988437  ,     2.22044605E-16 , 0.17364817766693041 ]
      call calc(inp)

      inp = [ 1.3865453490508738E-011 , 0.99999499320987351 , 0.0000000000000000 , 0.17364817766693041 ]
      call calc(inp)

      inp = [ 113.59224626216431 , 2.7005225550306174E-008 , 0.0000000000000000  , 0.17364817766693041 ]
      call calc(inp)

      contains
        subroutine calc(inp)
            real(ireals) :: inp(4)
            real(ireals) :: out1(7),out2(7)

            print *,'inp ::',inp
            call eddington_coeff_zdun(inp(1),inp(2),inp(3),inp(4),out2(1),out2(2),out2(3),out2(4),out2(5),out2(6),out2(7))
            print *,'zdun ::',out2
!            call eddington_coeff_fab(inp(1),inp(2),inp(3),inp(4),out1(1),out1(2),out1(3),out1(4),out1(5),out1(6),out1(7))
            print *,'fab  ::',out1
            print *,''
        end subroutine

end program
