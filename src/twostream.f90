module m_twostream
      use data_parameters, only: ireals,iintegers,zero,one
      use eddington, only: rodents
      implicit none

      private
      public twostream_dir,twostream_diff

    contains

      subroutine twostream_dir(dtau,w0,g,mu0,incSolar,albedo, S)
        real(ireals),intent(in) :: dtau(:),w0(:),g(:),albedo,mu0,incSolar
        real(ireals),intent(out):: S(size(dtau)+1 ,3) ! direct + 2 diffuse src terms ( Eup,Edn ); direct radiation is irradiance through plane

        real(ireals) :: eddington(5,size(dtau))

        integer :: k,ke

        ke = size(dtau)

        do k=1,ke
          call rodents(dtau(k),w0(k),g(k),mu0, eddington(:,k) )
        enddo

        S(1,1) = incSolar
        do k=1,ke
          S(k+1,1) = S(k,1) * eddington(5,k)
        enddo

        do k=1,ke
          S(k  ,2) = S(k,1) * eddington(3,k)
          S(k+1,3) = S(k,1) * eddington(4,k)
        enddo

        ! Eup at surface also uses albedo *direct radiation - the term albedo*diffuse is done in the matrix solver part
        S(ke+1 ,2) = S(ke+1,1) * albedo

        S(1,3) = zero ! no contribution to Edn at TOA

        S(:,1) = S(:,1) * mu0

      end subroutine


      subroutine twostream_diff(dtau,w0,g,S,albedo, E)
          real(ireals),intent(in) :: dtau(:),w0(:),g(:),S(:,:),albedo
          real(ireals),intent(out) :: E(:,:) ! dim should be Nz+1, 2 -> 2 for Eup,Edn

          real(ireals),allocatable :: AB (:,:)

          real(ireals),allocatable :: B (:,:)
          real(ireals) :: eddington(5,size(dtau)), R,T

          integer :: N, KL, KU, NRHS, LDAB, LDB, INFO, m
          integer,allocatable :: IPIV(:)

          integer :: i,j,k,ke,bi

          ke = size(dtau)
          do k=1,ke
            call rodents(dtau(k),w0(k),g(k),one, eddington(:,k) )
          enddo

          N = 2*(ke+1)
          KL =2
          KU =2
          NRHS=1
          LDAB = 2*KL+KU+1
          LDB  = N

          allocate( IPIV(N) )

          allocate( AB (LDAB,N) )
          allocate( B (LDB,NRHS)   )

          ! Setup src vector
          B = zero
          do k=1,ke+1
            B(2*(k-1)+1,1) = -S(k,2) ! Eup 
            B(2*(k-1)+2,1) = -S(k,3) ! Edn 
          enddo

          ! Set self(diagonal) values to -1
          AB( :,: ) = zero
          do j=1,N
            i=j
            bi= KL+KU+1+i-j 
            AB( bi,j ) = -one
          enddo

          ! BC at grnd
           j=N-1 ; i=j+1 ; bi= KL+KU+1+i-j 
           AB( bi,j) = albedo !Eup is Edn*albedo


          ! in between
          do k=1,ke
            T = eddington(1,k)
            R = eddington(2,k)

            i=2*k-1 ; j=i+2 ; bi= KL+KU+1+i-j 
            AB(bi,j) = T
            
            i=2*k-1 ; j=i+1 ; bi= KL+KU+1+i-j 
            AB(bi,j) = R
          enddo

          !Edn at TOA = 0
          do k=2,ke+1
            T = eddington(1,k-1)
            R = eddington(2,k-1)
            !Edn

            i=2*k ; j=i-2 ; bi= KL+KU+1+i-j 
            AB(bi,j) = T
            
            i=2*k ; j=i-1 ; bi= KL+KU+1+i-j 
            AB(bi,j) = R
          enddo

          if(kind(one).eq.kind(real(one)) ) then !single_precision
            call SGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
          else
            call DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
          endif

          if(INFO.ne.0) then
            print *,'INFO',INFO
            stop 'Error in twostream calculation - lapack returned Error'
          endif

!          do i=1,N
!            print *,'x',i,'::',B(i,:)/70**2,'IPIV',IPIV(i)
!          enddo
!          print *,'INFO',INFO
          E(1,2) = zero ! Edn at TOA is 0

          do k=1,ke+1
            E(k,1) = B(2*k-1,NRHS) ! Eup
            E(k,2) = B(2*k,NRHS) ! Edn
!            print *,'setting value for Eup,Edn',k,' indices',2*(k-1)+1,2*(k-1)+2, E(k,1)/70**2,E(k,2)/70**2
          enddo

!         call exit()

        end subroutine


end module
