module debug
    use mpi
    use ed_basis, only: basis_t, ed_basis_get
    use dmft_params, only: nsite, norb, nbath, nspin, U, Up, Jex, Jp, mu
    use timer, only: timestamp
    use numeric_utils, only: mpi_norm, mpi_dot_product
    use impurity_hamiltonian, only: hx

    implicit none

contains

    subroutine dump_hamiltonian(em, ek, vmk, basis)
        double precision, intent(in) :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)
        type(basis_t), intent(in) :: basis

        double precision, allocatable :: x(:), y(:), x_all(:)
        integer :: irow, icol
        integer :: b

        if (nprocs>1) then
            stop "only serial run is allowed"
        endif

        open(unit=137, file="h_imp.dump", status="replace")
        allocate( x(basis%ntot), y(basis%ntot), x_all(basis%ntot) )

        do irow=1,basis%ntot
            x = 0.d0
            x(irow) = 1.d0

            call hx( em, ek, vmk, basis, x, y, x_all, .false. )

            do icol=1,basis%ntot
                write(137,"(F6.3,1x)",advance="no") y(icol)
            enddo
            write(137,*)
        enddo
        close(137)
        open(unit=137, file="basis.dump", status="replace")
        do irow=1,basis%ntot
            b = ed_basis_get(basis, irow)
            write(137,"(I5,2x)",advance="no") irow
            do icol=1,2*nsite
                if (BTEST(b,icol-1)) then
                    write(137,"(I1)",advance="no") 1
                else
                    write(137,"(I1)",advance="no") 0
                endif
            enddo
            write(137,*)
        enddo
        close(137)
    end subroutine dump_hamiltonian

end module debug
