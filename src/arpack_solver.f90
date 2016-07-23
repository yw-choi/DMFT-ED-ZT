module arpack_solver
    use mpi
    use utils, only: die
    use ed_basis, only: basis_t, generate_basis
    use dmft_params, only: norb, nbath, maxnstep
    use timer, only: timestamp
    use numeric_utils, only: mpi_norm, mpi_dot_product
    use impurity_hamiltonian, only: hx

    implicit none

    public :: diag_arpack

    private
contains

    ! Diagonalize a sector of the impurity Hamiltonian.
    ! basis has the information of the relevant sector.
    subroutine diag_arpack(em, ek, vmk, basis, eigval, eigvec)
        double precision, intent(in) :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)
        type(basis_t), intent(in) :: basis

        double precision, intent(out) :: eigval
        double precision, intent(out) :: eigvec(basis%nloc)

        integer, parameter :: nev = 5
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

        integer :: i, j, ierr, multcount

        allocate( v(basis%nloc,2*nev) )
        allocate( workd(3*basis%nloc) )
        allocate( resid(basis%nloc) )
        allocate( x_all(basis%ntot) )

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
                call hx( em, ek, vmk, basis, &
                         workd(ipntr(1)), workd(ipntr(2)), x_all, .false. )
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

                allocate(ax(basis%nloc))

                do j=1, nconv
                    call hx( em, ek, vmk, basis, &
                             v(:,j), ax, x_all, .false. )
                    call daxpy(basis%nloc, -d(j,1), v(1,j), 1, ax, 1)
                    d(j,2) = pdnorm2( comm, basis%nloc, ax, 1 )
                enddo

                call pdmout(comm, 6, nconv, 2, d, ncv, -6, &
                    'Ritz values and direct residuals')

                deallocate(ax)
            end if
        endif

        if (nconv /= nev) then
            call die("ed_diag_arpack", "diagonalization not converged properly.")
            return
        endif

        eigval = d(1,1)
        eigvec = v(:,1)

        deallocate(v,resid,workd,x_all)
    end subroutine diag_arpack
end module arpack_solver
