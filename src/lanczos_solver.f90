module lanczos_solver
    use mpi
    use ed_basis, only: basis_t, generate_basis
    use dmft_params, only: norb, nbath, maxnstep
    use timer, only: timestamp
    use numeric_utils, only: mpi_norm, mpi_dot_product
    use impurity_hamiltonian, only: hx

    implicit none

    public :: &
        diag_lanczos, &
        lanczos_iteration

    private
contains

    ! Diagonalize a sector of the impurity Hamiltonian.
    ! basis has the information of the relevant sector.
    subroutine diag_lanczos(em, ek, vmk, basis, eigval, eigvec)
        double precision, intent(in) :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)
        type(basis_t), intent(in) :: basis

        double precision, intent(out) :: eigval
        double precision, intent(out) :: eigvec(basis%nloc)

        integer :: nstep
        double precision, allocatable :: &
            a(:), b(:), v_init(:), &
            ev(:), lanczos_v(:,:), coeff(:)

        integer :: i
        double precision :: t1, t2, r, residual
        
        allocate(v_init(basis%nloc))
        allocate(a(maxnstep),b(maxnstep))
        
        ! initial random v
        call random_seed
        do i=1,basis%nloc
            call random_number(r)
            v_init(i) = r
        enddo

        t1 = mpi_wtime(mpierr)
        call lanczos_iteration( em, ek, vmk, basis, v_init, &
                               nstep, a, b )

        t2 = mpi_wtime(mpierr)

        if (master) then
            print *, "nstep = ", nstep
            print *, "lanczos iteration time = ", (t2-t1), " sec."
        endif

        allocate(ev(nstep),lanczos_v(nstep,nstep))
        allocate(coeff(nstep))

        if (master) then
            t1 = mpi_wtime(mpierr)
            call lanczos_diagonalize( nstep, a, b, ev, lanczos_v )
            t2 = mpi_wtime(mpierr)

            eigval = ev(1)
            coeff = lanczos_v(:,1)
            
            print *, "lanczos diagonalization time = ", (t2-t1), " sec."
            print *, "lanczos ground state energy = ", eigval
        endif

        call mpi_bcast( coeff, nstep, mpi_double_precision, 0, comm, mpierr )
        call mpi_bcast( eigval, 1, mpi_double_precision, 0, comm, mpierr )

        eigvec = 0.d0

        t1 = mpi_wtime(mpierr)
        call lanczos_ground_state( em, ek, vmk, basis, v_init, &
                                   nstep, a, b, coeff, eigval, eigvec, residual)
        t2 = mpi_wtime(mpierr)

        if (master) then
            print *, "lanczos eigenvector calculation time = ", (t2-t1), " sec."
            print *, "residual = ", residual
        endif

        deallocate(v_init,a,b,ev,lanczos_v,coeff)
    end subroutine diag_lanczos

    subroutine lanczos_iteration( em, ek, vmk, basis, v_init, &
                                 nstep, a, b )

        double precision, intent(in) :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)
        type(basis_t), intent(in) :: basis
        double precision, intent(in) :: v_init(basis%nloc)

        integer, intent(out) :: nstep
        double precision, intent(out) :: a(maxnstep), b(maxnstep)

        ! temporary lanczos vectors
        double precision, allocatable :: v(:), w(:), x_all(:)
        double precision :: norm_v, t  
        
        integer :: i, j, ierr, k, nloc
        character(len=100) :: msg

        nloc = basis%nloc

        allocate(v(nloc), w(nloc), x_all(basis%ntot))
        
        a = 0.0D0
        b = 0.0D0

        ! Lanczos steps
        ! ref: G. Golub, Matrix Computations, 4th ed., p.562 (2013)

        ! if (master) then
        !     write(msg,"(a,I4)") "eigenvalue Lanczos iteration ", 1
        !     call timestamp(msg)
        ! endif

        ! normalize the initial vector
        norm_v = mpi_norm( v_init, nloc)

        w = v_init/norm_v

        call hx( em, ek, vmk, basis, w, v, x_all, .false. )

        a(1) = mpi_dot_product( w, v, nloc )
        b(1) = 0.d0

        if (maxnstep==1) then
            nstep = 1
            deallocate(v, w, x_all)
            return
        endif

        call daxpy( nloc, -a(1), w, 1, v, 1 )
        b(2) = mpi_norm( v, nloc )

        k = 2
        do while (1)
            ! if (master) then
            !     write(msg,"(a,I4)") "eigenvalue Lanczos iteration", k
            !     call timestamp(msg)
            ! endif
            do i=1,nloc
                t = w(i)
                w(i) = v(i)/b(k)
                v(i) = -b(k)*t
            enddo

            call hx( em, ek, vmk, basis, w, v, x_all, .true. )
            a(k) = mpi_dot_product( w, v, nloc ) 

            if (k>=maxnstep.or.b(k)<1.d-16) then
                nstep = k
                exit
            endif

            call daxpy( nloc, -a(k), w, 1, v, 1 )
            b(k+1) = mpi_norm( v, nloc )

            k = k+1
        enddo

        deallocate(v,w,x_all)
    end subroutine lanczos_iteration

    subroutine lanczos_diagonalize( nstep, a, b, ev, v )
        include 'mkl_lapack.fi'
        integer, intent(in) :: nstep
        double precision, intent(in) :: a(nstep), b(nstep)
        double precision, intent(out) :: ev(nstep), v(nstep,nstep)

        double precision :: d(nstep), e(nstep-1), abstol, &
                            work(5*nstep)
        integer :: il, iu, m, iwork(5*nstep), ifail(nstep), &
                   info
        double precision :: w(nstep), z(nstep,nstep)

        integer :: i

        ! copy a
        d = a
        e(1:nstep-1) = b(2:nstep)

        il = 1                 ! lowest eigenvalue index
        iu = nstep             ! highest eigenvalue index
        abstol = 2*DLAMCH('S') ! recommended tolerance

        call dstevx( 'V', 'A', nstep, d, e, 0.d0, 0.d0, il, iu, abstol, &
                     m, w, z, nstep, work, iwork, ifail, info )

        if (info/=0) then
            write(*,*) "dstevx error. info = ", info
            if (info>0) then
                write(*,*) "ifail = "
                write(*,*), ifail
            endif
            stop
        endif

        ev(1:m) = w(1:m)
        v(1:nstep,1:m) = z(1:nstep,1:m) 
    end subroutine lanczos_diagonalize

    subroutine lanczos_ground_state( em, ek, vmk, basis, v_init, &
                                nstep, a, b, coeff, eigval, eigvec, residual)
        double precision, intent(in) :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)
        type(basis_t), intent(in) :: basis
        integer, intent(in) :: nstep

        double precision, intent(in) :: &
            v_init(basis%nloc), &   ! starting vector for lanczos iteration
            eigval,             &  
            a(nstep),           &   ! diagonal matrix element in lanczos basis
            b(nstep),           &   ! off-diagonal matrix element in lanczos basis
            coeff(nstep)            ! ground eigenvector of the lanczos matrix 

        double precision, intent(out) :: &
            eigvec(basis%nloc), residual

        ! temporary lanczos vectors
        double precision, allocatable :: v(:), w(:), x_all(:)
        double precision :: norm_v, t
        
        integer :: i, j, ierr, k, nloc
        character(len=100) :: msg

        nloc = basis%nloc

        allocate(v(nloc), w(nloc), x_all(basis%ntot))

        ! if (master) then
        !     write(msg,"(a,I4)") "eigenvector Lanczos iteration", 1
        !     call timestamp(msg)
        ! endif

        ! normalize the initial vector
        norm_v = mpi_norm( v_init, nloc)

        w = v_init/norm_v ! v1

        ! eigvec = coeff(1)*w
        eigvec = 0.d0
        call daxpy( nloc, coeff(1), w, 1, eigvec, 1) 

        call hx( em, ek, vmk, basis, w, v, x_all, .false. )

        if (nstep==1) then
            deallocate(v, w, x_all)
            return
        endif

        call daxpy( nloc, -a(1), w, 1, v, 1 )

        k = 2
        do while (1)
            ! if (master) then
            !     write(msg,"(a,I4)") "eigenvector Lanczos iteration", k
            !     call timestamp(msg)
            ! endif
            do i=1,nloc
                t = w(i)
                w(i) = v(i)/b(k)
                v(i) = -b(k)*t
            enddo
            ! w = |v_k>
            ! eigvec = eigvec + coeff(k)*w
            call daxpy( nloc, coeff(k), w, 1, eigvec, 1 ) 

            call hx( em, ek, vmk, basis, w, v, x_all, .true. )

            if (k>=nstep.or.b(k)<1.d-16) then
                exit
            endif
            
            call daxpy( nloc, -a(k), w, 1, v, 1 )
            k = k+1
        enddo

        ! test eigenvector residual | H|gs>-E0|gs> |^2
        call hx( em, ek, vmk, basis, eigvec, v, x_all, .false. )
        call daxpy( nloc, -eigval, eigvec, 1, v, 1)

        residual = mpi_dot_product( v, v, nloc )        

        deallocate(v,w,x_all)
    end subroutine lanczos_ground_state
end module lanczos_solver
