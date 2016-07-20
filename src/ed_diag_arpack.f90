module ed_diag_arpack
!==============================================================================
! Parallel diagonalization using Parallel ARPACK library.
!==============================================================================

    use mpi
    use dmft_params, only: norb, beta
    use ed_params, only: nev, nsite, nbath, sectors, nsector, &
                         PROB_THRESHOLD, print_arpack_stat
    use ed_basis, only: generate_basis, basis_t, ed_basis_get, dealloc_basis
    use ed_hamiltonian, only: multiply_H
    use numeric_utils, only: boltzmann_factor, sort
    use utils, only: die
    use timer, only: t1_diag_loop, t2_diag_loop, print_elapsed_time
    use alloc, only: re_alloc, de_alloc, alloc_count
    use ed_eigpair, only: eigpair_t, eigpairs, nev_calc

    implicit none

    public :: diag_arpack

    private

contains

    subroutine diag_arpack(ia)
        integer, intent(in) :: ia

        type(basis_t) :: basis
        integer :: ne_up, ne_down, nh, nloc, isector, iev, ind(nev*nsector),i, &
                   multcount
        double precision :: eigval(nev*nsector), prob(nev*nsector)
        type(eigpair_t) :: eigpairs_all(nev, nsector)

        double precision :: ground, Z

        if (master) then
            write(*,*) "ARPACK Diagonalization"
        endif

        do isector=1,nsector
            ne_up = sectors(isector,1)
            ne_down = sectors(isector,2)
            nh = sectors(isector,3)
            if (master) then
                write(*,"(a,I3)") "    diagonalizing sector ",isector
                write(*,"(a,I2,a2,I2,a2,I,a)")"    (ne_up, ne_dn, dim) = (",&
                                    ne_up,", ",ne_down,", ",nh,")" 
            endif

            call generate_basis(ne_up, ne_down, basis)

            t1_diag_loop = mpi_wtime(mpierr)    
            call diagonalization(ia, isector, basis, eigpairs_all(:, isector),&
                                multcount)
            t2_diag_loop = mpi_wtime(mpierr)    

            if (master) then
                write(*,*) "    matrix-vector multiplication count = ", multcount
                call print_elapsed_time("    diagonalization time", &
                                        t1_diag_loop,t2_diag_loop)
            endif

            do iev=1, nev
                eigval((isector-1)*nev+iev) = eigpairs_all(iev,isector)%val
            enddo

            call dealloc_basis(basis)
        enddo

        ! sort all the eigenvalues
        call sort(nev*nsector,eigval,ind)

        ! boltzman factors and partition function
        Z = 1.d0
        prob(1) = 1.d0
        do i=2,nev*nsector
            prob(i) = exp(-beta*(eigval(i)-eigval(1)))
            Z = Z + prob(i)
        enddo
        prob = prob / Z

        ! how many states have significant boltzman factor?
        nev_calc = 0
        do iev=1,nev
            if (prob(ind(iev)) > PROB_THRESHOLD) then
                nev_calc = nev_calc + 1
            endif
        enddo

        ! take lowest nev_calc eigenvalues
        do i=1,nev_calc
            isector = (ind(i)-1)/nev+1
            iev = mod(ind(i)-1,nev)+1
            eigpairs(i)%sector =  eigpairs_all(iev,isector)%sector
            eigpairs(i)%level  =  eigpairs_all(iev,isector)%level
            eigpairs(i)%val    =  eigpairs_all(iev,isector)%val
            eigpairs(i)%prob   =  prob(ind(i))
            eigpairs(i)%nloc   =  eigpairs_all(iev,isector)%nloc
            ! This is NOT copying the array.
            ! Only the pointer is set to eigpairs
            eigpairs(i)%vec    => eigpairs_all(iev,isector)%vec
        enddo

        ! and free the rest
        ! do i=nev_calc+1,nev*nsector
        !     isector = (ind(i)-1)/nev+1
        !     iev = mod(ind(i)-1,nev)+1
        !     call de_alloc(eigpairs_all(iev,isector)%vec, &
        !         "ed_diag_arpack", "eigvec_all")
        ! enddo

        if (master) then
            write(6,*)
            write(6,"(a,ES10.3,a,I5)") " Number of eigenvalues (with prob > ", &
                                       PROB_THRESHOLD,") = ", nev_calc
            write(6,*)
            write(6,"(a)") " Eigenvalue          Prob           Sector   Level"
            do iev=1,nev_calc
                write(6,"(1x,F16.10,4x,ES13.5,2I8)") &
                    eigpairs(iev)%val, eigpairs(iev)%prob,&
                    eigpairs(iev)%sector, eigpairs(iev)%level
            enddo
            write(6,*)
        endif

    end subroutine diag_arpack

    subroutine diagonalization(ia, isector, basis, eig, multcount)
        integer, intent(in) :: isector, ia
        type(basis_t), intent(in) :: basis
        type(eigpair_t), intent(out) :: eig(nev)
        integer, intent(out) :: multcount

        character, parameter :: bmat = 'I'
        character(len=2), parameter :: which = 'SA'

        double precision, allocatable :: ax(:), workd(:), resid(:), &
                                         x_all(:), v(:,:)
        double precision :: workl(2*nev*(2*nev+8)), d(2*nev,2), sigma, tol

        logical :: select(2*nev)
        integer :: iparam(11), ipntr(11), lworkl, info, ido, nconv, ncv, &
                   maxitr, mode, ishfts, ldv

        double precision :: pdnorm2
        external :: pdnorm2, daxpy

        integer :: i, j, ierr

        allocate( v(basis%nloc,2*nev) )
        allocate( workd(3*basis%nloc) )
        allocate( resid(basis%nloc) )
        allocate( x_all(basis%ntot) )
        call alloc_count(2*nev*basis%nloc, 'D', 'ed_diag_arpack', 'v')
        call alloc_count(3*basis%nloc, 'D', 'ed_diag_arpack', 'workd')
        call alloc_count(basis%nloc, 'D', 'ed_diag_arpack', 'resid')
        call alloc_count(basis%ntot, 'D', 'ed_diag_arpack', 'x_all')

        ldv = basis%nloc
        ncv = 2*nev
        lworkl = ncv*(ncv+8)

        tol = 0.0

        ishfts = 1
        maxitr = 20000
        mode   = 1

        iparam(1) = ishfts
        iparam(3) = maxitr
        iparam(7) = mode

        info = 0
        ido = 0
        multcount = 0
        do
            call pdsaupd( comm, ido, bmat, basis%nloc, which, nev, tol, resid, &
                ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
            if (ido .eq. -1 .or. ido .eq. 1) then
                call multiply_H(ia, basis, workd(ipntr(1)), workd(ipntr(2)), &
                                x_all)
                multcount = multcount+1
            else
                exit
            endif
        enddo

        if ( info .lt. 0 ) then
            print *, ' Error with pdsaupd, info = ', info
            print *, iparam(5)
            call die("diagonalization", "PDSAUPD ERROR")
        else
            call pdseupd(comm, .true., 'All', select, d, v, ldv, sigma, &
                bmat, basis%nloc, which, nev, tol, resid, ncv, v, ldv, &
                iparam, ipntr, workd, workl, lworkl, ierr )
            if ( ierr .ne. 0) then
                print *, ' Error with pdseupd, info = ', ierr
                call die("diagonalization", "pdseupd ERROR")
            else
                nconv =  iparam(5)

                if (print_arpack_stat) then
                    allocate(ax(basis%nloc))
                    call alloc_count(basis%nloc, 'D', &
                        'ed_diag_arpack', 'ax')

                    do j=1, nconv
                        call multiply_H(ia, basis, v(:,j), ax, x_all)
                        call daxpy(basis%nloc, -d(j,1), v(1,j), 1, ax, 1)
                        d(j,2) = pdnorm2( comm, basis%nloc, ax, 1 )
                    enddo

                    call pdmout(comm, 6, nconv, 2, d, ncv, -6, &
                        'Ritz values and direct residuals')

                    deallocate(ax)
                    call alloc_count(-basis%nloc, 'D', &
                        'ed_diag_arpack', 'ax')
                endif 
            end if
        endif

        if (nconv /= nev) then
            call die("ed_diag_arpack", "diagonalization not converged properly.")
            return
        endif

        do i=1,nconv
            eig(i)%sector = isector
            eig(i)%level  = i
            eig(i)%val    = d(i,1)

            nullify(eig(i)%vec)
            call re_alloc(eig(i)%vec, 1, basis%nloc, &
                            "ed_diag_arpack", "eigvec_all")
            eig(i)%vec    = v(:,i)
            eig(i)%nloc   = basis%nloc
        enddo

        deallocate(v,resid,workd,x_all)
        call alloc_count(-2*nev*basis%nloc, 'D', 'ed_diag_arpack', 'v')
        call alloc_count(-3*basis%nloc, 'D', 'ed_diag_arpack', 'workd')
        call alloc_count(-basis%nloc, 'D', 'ed_diag_arpack', 'resid')
        call alloc_count(-basis%ntot, 'D', 'ed_diag_arpack', 'x_all')
    end subroutine diagonalization
end module ed_diag_arpack
