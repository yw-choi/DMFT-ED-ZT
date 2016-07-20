module ed_operator
    use dmft_params, only: norb
    use ed_params, only: nbath, nsite, kind_basis
    use mpi
    use ed_basis, only: basis_t, ed_basis_idx, ed_basis_get

    implicit none

contains

    ! Apply creation/destruction operator to a state vector.
    ! pm = 1 : create
    ! pm = 2 : destroy
    subroutine apply_c(basis, vec, vec_all, pm, iorb, ispin, basis_out, vec_out)
        type(basis_t), intent(in) :: basis
        double precision, intent(in) :: vec(basis%nloc)
        integer, intent(in) :: pm, iorb, ispin
        type(basis_t), intent(in) :: basis_out
        double precision, intent(out) :: vec_out(basis_out%nloc), &
                                         vec_all(basis%ntot) 
        integer(kind=kind_basis) :: basis_i, basis_j
        integer :: i,j,sgntot, isite

        vec_out = 0.0D0

        call mpi_allgatherv(vec,basis%nloc,mpi_double_precision,vec_all,&
            basis%nlocals,basis%offsets,mpi_double_precision,comm,mpierr)

        do i=1,basis_out%nloc
            basis_i = ed_basis_get(basis_out,i)
            isite = (ispin-1)*nsite+iorb
            sgntot = sgn(basis_i,1,isite-1)
            ! transposed. 1=creation, 2=destruction
            if (pm.eq.1) then
                if (.not.BTEST(basis_i,isite-1)) then
                    sgntot = 0
                endif
                basis_j = IBCLR(basis_i, isite-1)
            else
                if (BTEST(basis_i,isite-1)) then
                    sgntot = 0
                endif
                basis_j = IBSET(basis_i, isite-1)
            endif

            if (sgntot.eq.0) then
                cycle
            endif

            j = ed_basis_idx(basis, basis_j)
            vec_out(i) = vec_out(i) + vec_all(j)*sgntot
        enddo
    end subroutine apply_c

    integer function sgn(basis_in,i,j)
        integer(kind=kind_basis) :: basis_in
        integer :: i,j, k, sgnsum

        sgnsum = 0
        do k=i,j
            if (BTEST(basis_in,k-1)) then
                sgnsum = sgnsum + 1
            endif
        enddo

        if (mod(sgnsum,2)) then
            sgn = -1
        else
            sgn = +1
        endif
        return
    end function sgn
end module ed_operator
