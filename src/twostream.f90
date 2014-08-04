module m_twostream
      use m_data_parameters, only: ireals,iintegers,zero,one
      use m_eddington, only: eddington_coeff_fab
      implicit none

      private
      public delta_eddington_twostream

    contains

      subroutine delta_eddington_twostream(dtau,w0,g,mu0,incSolar,albedo, S,Edn,Eup)
        real(ireals),intent(in) :: dtau(:),w0(:),g(:),albedo,mu0,incSolar
        real(ireals),dimension(size(dtau)+1),intent(out):: S,Edn,Eup 

        real(ireals),dimension(size(dtau)) :: a11,a12,a13,a23,a33

        integer(iintegers) :: i,j,k,ke,ke1,bi
        real(ireals) :: R,T
        real(ireals),allocatable :: AB (:,:)
        real(ireals),allocatable :: B (:,:)
        integer,allocatable :: IPIV(:)
        integer :: N, KLU,  KL, KU, NRHS, LDAB, LDB, INFO

        ke = size(dtau)
        ke1 = ke+1
        N = 2*(ke1)
        KL =2
        KU =2
        NRHS=1
        LDAB = 2*KL+KU+1
        LDB  = N
        KLU = KL+KU+1

        do k=1,ke
          call eddington_coeff_fab (dtau(k), w0(k),g(k), mu0,a11(k),a12(k),a13(k),a23(k),a33(k))
!          print *,'eddington',k,' :: ',dtau(k), w0(k),g(k), mu0,'::',a11(k),a12(k),a13(k),a23(k),a33(k)
        enddo

        S(1) = incSolar *mu0 ! irradiance on tilted plane
        do k=1,ke
          S(k+1) = S(k) * a33(k)
        enddo

        allocate( IPIV(N) )
        allocate( AB (LDAB,N) )
        allocate( B (LDB,NRHS)   )
        AB   = zero
        B    = zero
        IPIV = 0

        ! Setup src vector
        do k=1,ke
          B(2*k-1,1) = S(k) * a13(k) ! Eup 
          B(2*k+2,1) = S(k) * a23(k) ! Edn 
        enddo
        B(2,1) = zero ! no Edn at TOA
        B(2*ke1-1,1) = S(ke1) * albedo 

        !diagonal entries
        do i=1,N
          j=i
          bi= KLU+i-j 
          AB( bi,j ) = one
        enddo

        do k=1,ke
          T = a11(k)
          R = a12(k)

        ! setting Eup coeffs
          i=2*k-1 ; j=i+2 ; bi= KLU+i-j ; AB(bi,j) = -T
          i=2*k-1 ; j=i+1 ; bi= KLU+i-j ; AB(bi,j) = -R

        ! setting Edn coeffs
          i=2*(k+1) ; j=i-2 ; bi= KLU+i-j ; AB(bi,j) = -T
          i=2*(k+1) ; j=i-1 ; bi= KLU+i-j ; AB(bi,j) = -R
        enddo
        i=2*ke1-1 ; j=i+1 ; bi= KLU+i-j ; AB(bi,j) = -albedo ! Eup at surface is Edn*albedo

        do k=1,ke
          if(any(isnan( [a11(k),a12(k),a13(k),a23(k),a33(k)] )) ) then
            print *,'eddington coefficients',k,' source',B(2*k-1,1),B(2*k,1), 'eddington',a11(k),a12(k),' :: ',a13(k),a23(k), ' :: ',a33(k), ' :: ',dtau(k),w0(k),g(k),'S=',S(k)
            call exit()
          endif
        enddo
        
        INFO=-1
        if(kind(one).eq.kind(real(one)) ) then !single_precision
          call SGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
        else if(kind(one).eq.kind(dble(one)) ) then !double_precision
          call DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
        else
          print *,'Dont know which LAPACK routine to call for real kind',kind(one)
          call exit(-5)
        endif

        if(INFO.ne.0) then
          print *,'INFO',INFO
          stop 'Error in twostream calculation - lapack returned Error'
        endif

        ! retrieve result from solver
        do k=1,ke1
          Eup(k) = B(2*k-1,NRHS) ! Eup
          Edn(k) = B(2*k,NRHS) ! Edn
          if(any(isnan( [Eup(k),Edn(k)] )) ) &
              print *,'setting value for Eup,Edn',k,' Source',B(2*k-1,1),B(2*k,1),'Eup/dn', Eup(k),Edn(k),'IPIV',IPIV(2*k-1:2*k)
        enddo

!        S = S*mu0

!        call exit()

        end subroutine


end module
