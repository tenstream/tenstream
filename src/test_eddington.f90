program main
      use m_eddington
      use mpi
      use m_data_parameters, only : mpiint,ireals,iintegers,one,zero, init_mpi_data_parameters,i1

      real(ireals) :: inp(4) ! tau, omega0, g, mu0
      real(ireals) :: out(7) ! a11,a12,a13,a23,a33,g1,g2

      call MPI_Init(ierr)
      call init_mpi_data_parameters(MPI_COMM_WORLD)

      print *,'Checking eddington coefficients'

      inp = [ 3.2662250689525390E-011 , 0.99999171495417127 ,       0.0000000000000000  , 0.17364817766693041 ]
      call eddington_coeff_fab(inp(1),inp(2),inp(3),inp(4),out(1),out(2),out(3),out(4),out(5),out(6),out(7))
      print *,inp,'::',out

      inp = [ 2.9317851124478626E-012,   1.0000000000000000  ,      0.0000000000000000 , 0.17364817766693041 ]
      call eddington_coeff_fab(inp(1),inp(2),inp(3),inp(4),out(1),out(2),out(3),out(4),out(5),out(6),out(7))
      print *,inp,'::',out


      inp = [1.93321303E-10,  0.999984443 ,      2.22044605E-16,  0.17364817766693041 ]
      call eddington_coeff_fab(inp(1),inp(2),inp(3),inp(4),out(1),out(2),out(3),out(4),out(5),out(6),out(7))
      print *,inp,'::',out

      inp = [7.89528581E-11,  0.999988437  ,     2.22044605E-16 , 0.17364817766693041 ]
      call eddington_coeff_fab(inp(1),inp(2),inp(3),inp(4),out(1),out(2),out(3),out(4),out(5),out(6),out(7))
      print *,inp,'::',out
      
      inp = [ 1.3865453490508738E-011 , 0.99999499320987351 , 0.0000000000000000 , 0.17364817766693041 ]
      call eddington_coeff_fab(inp(1),inp(2),inp(3),inp(4),out(1),out(2),out(3),out(4),out(5),out(6),out(7))
      print *,inp,'::',out

      inp = [ 113.59224626216431 , 2.7005225550306174E-008 , 0.0000000000000000  , 0.17364817766693041 ]
      call eddington_coeff_fab(inp(1),inp(2),inp(3),inp(4),out(1),out(2),out(3),out(4),out(5),out(6),out(7))
      print *,inp,'::',out



end program
