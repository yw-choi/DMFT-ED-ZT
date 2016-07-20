module numeric_utils

contains
    ! f =                      b(1)*b(1)
    !         ---------------------------------------------
    !                                 b(2)*b(2)  
    !          z - a(1) - ---------------------------------
    !                                          b(3)*b(3)
    !                         z - a(2)  -  ----------------
    !                                       z - a(3) - ... 
    double complex function continued_fraction_p(z,n,a,b) result(f)
        double complex, intent(in) :: z
        integer, intent(in) :: n
        double precision, intent(in) :: a(n), b(n)

        integer :: i 

        f = b(n)*b(n)/(z-a(n))
        do i = n-1, 2, -1
            f = b(i)*b(i)/(z-a(i)-f)
        enddo
        f = b(1)*b(1)/(z-a(1)-f)
    end function continued_fraction_p

    ! f =                      b(1)*b(1) 
    !         ---------------------------------------------
    !                                 b(2)*b(2)  
    !          z + a(1) - ---------------------------------
    !                                          b(3)*b(3)
    !                         z + a(2)  -  ----------------
    !                                       z + a(3) - ... 
    double complex function continued_fraction_m(z,n,a,b) result(f)
        double complex, intent(in) :: z
        integer, intent(in) :: n
        double precision, intent(in) :: a(n), b(n)

        integer :: i 

        f = b(n)*b(n)/(z+a(n))
        do i = n-1, 2, -1
            f = b(i)*b(i)/(z+a(i)-f)
        enddo
        f = b(1)*b(1)/(z+a(1)-f)
    end function continued_fraction_m

    double precision function boltzmann_factor(beta,de) result(prob)
        double precision, intent(in) :: beta, de

        prob = exp(-beta*de)

    end function boltzmann_factor

    ! Calculates the inverse matrix of a complex NxN matrix. 
    ! subroutine complex_inverse(N,A,INFO)
    !     double complex :: A(N,N), WORK(N)
    !     integer :: INFO, N, IPIV(N)

    !     call ZGETRF(N,N,A,N,IPIV,INFO)

    !     if (INFO.ne.0) then
    !         return
    !     endif

    !     call ZGETRI(N,A,N,IPIV,WORK,N,INFO)

    ! end subroutine complex_inverse

    subroutine cinv( A, N, MAXN, Ainv )

        integer*4 N, MAXN
        double complex A(MAXN,MAXN), Ainv(MAXN,MAXN)

        ! Inputs
        !   A       Matrix A to be inverted
        !   N       Elements used in matrix A (N by N)
        !  MAXN     Matrix dimenstions as A(MAXN,MAXN)
        ! Outputs
        !  Ainv     Inverse of matrix A

        integer*4 MAXMAXN
        parameter( MAXMAXN = 200 )
        integer*4 i, j, k, index(MAXMAXN), jPivot, indexJ
        double precision :: scale(MAXMAXN), scaleMax, ratio, ratioMax
        double complex :: AA(MAXMAXN,MAXMAXN), B(MAXMAXN,MAXMAXN), coeff, sum

        if( MAXN .gt. MAXMAXN ) then
            write(*,*) 'ERROR in cinv: Matrix too large'
            stop
        endif

        !* Matrix B is initialized to the identity matrix
        do i=1,N
        do j=1,N
        AA(i,j) = A(i,j)  ! Copy matrix so as not to overwrite
        B(i,j) = 0.0
        enddo
        B(i,i) = 1.0
        enddo

        !* Set scale factor, scale(i) = max( |a(i,j)| ), for each row
        do i=1,N
        index(i) = i     ! Initialize row index list
        scaleMax = 0.0
        do j=1,N
        if( abs(AA(i,j)) .gt. scaleMax ) then
            scaleMax = abs(AA(i,j))
        endif
        enddo
        scale(i) = scaleMax
        enddo

        !* Loop over rows k = 1, ..., (N-1)
        do k=1,(N-1)
        !* Select pivot row from max( |a(j,k)/s(j)| )
        ratiomax = 0.0
        jPivot = k
        do i=k,N
        ratio = abs(AA(index(i),k))/scale(index(i))
        if( ratio .gt. ratiomax ) then
            jPivot=i
            ratiomax = ratio
        endif
        enddo
        !* Perform pivoting using row index list
        indexJ = index(k)
        if( jPivot .ne. k ) then     ! Pivot
            indexJ = index(jPivot)
            index(jPivot) = index(k)   ! Swap index jPivot and k
            index(k) = indexJ
        endif
        !* Perform forward elimination
        do i=k+1,N
        coeff = AA(index(i),k)/AA(indexJ,k)
        do j=k+1,N
        AA(index(i),j) = AA(index(i),j) - coeff*AA(indexJ,j)
        enddo
        AA(index(i),k) = coeff
        do j=1,N
        B(index(i),j) = B(index(i),j) - AA(index(i),k)*B(indexJ,j)
        enddo
        enddo
        enddo

        !* Perform backsubstitution
        do k=1,N
        Ainv(N,k) = B(index(N),k)/AA(index(N),N)
        do i=N-1,1,-1
        sum = B(index(i),k)
        do j=i+1,N
        sum = sum - AA(index(i),j)*Ainv(j,k)
        enddo
        Ainv(i,k) = sum/AA(index(i),i)
        enddo
        enddo

        return
    end subroutine cinv

    double precision function mpi_dot_product(A,B,n) result(dab)
        use mpi

        integer n
        double precision :: A(n), B(n), dab_tmp

        dab_tmp = sum(A(1:n)*B(1:n)) 
        call mpi_allreduce(dab_tmp,dab,1,mpi_double_precision,&
            mpi_sum,comm,ierr)

        return
    end function mpi_dot_product

    double precision function mpi_norm(A,n) result(dab)
        use mpi

        integer n
        double precision :: A(n)
        dab = sqrt(mpi_dot_product(A,A,n))
        return    
    end function mpi_norm

    integer function ifact(n)

        integer:: n,i

        ifact = 1
        do i = 1, n
        ifact = ifact*i
        enddo

        return
    end function ifact

    integer function iP(n,m)

        integer:: n,m,i

        iP = 1

        do i = n-m+1, n
        iP = iP*i
        enddo 

        return
    end function iP


    integer function iCom(n,m) 
        integer:: n, m, k

        if((n.eq.m).or.(m.eq.0)) then
            iCom = 1
        else
            if(m.gt.n/2) then
                k = n-m
                iCom = iP(n,k)/ifact(k)
            elseif(m.le.n/2) then
                iCom = iP(n,m)/ifact(m)
            endif 
        endif

        return
    end function icom

    SUBROUTINE SORT(N,ARRIN,INDX)

!     SORTS AN ARRAY BY THE HEAPSORT METHOD
!     W. H. PREUSS ET AL. NUMERICAL RECIPES
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     this routine is not useful for samll size of n.
!     in particular n=1 lead to operation error.(arrin(0) is not defined.)
!     if n is greater than ~ 20 this routine is more fast.
!
      DIMENSION ARRIN(N),INDX(N)
!     DIMENSION ARRIN(1),INDX(1)
      IF (N .EQ. 1) THEN
        INDX(1) = 1
        RETURN
      ENDIF
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
    END SUBROUTINE SORT
end module numeric_utils
