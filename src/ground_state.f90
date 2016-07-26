subroutine ground_state(em, ek, vmk, gs)
    
    use mpi
    use utils, only: die
    use eigpair, only: eigpair_t
    use ed_basis, only: basis_t, generate_basis, dealloc_basis
    use dmft_params, only: norb, nbath, nsite, nspin, nsector, sectors, &
                            diag_solver
    use timer, only: t1_diag_loop, t2_diag_loop, print_elapsed_time
    use lanczos_solver, only: diag_lanczos
    use arpack_solver, only: diag_arpack

    use debug, only: dump_hamiltonian

    double precision, intent(in) :: &
        em(norb,2), &
        ek(norb,2), &
        vmk(norb,nbath,2)

    type(eigpair_t), intent(out) :: gs

    integer :: isector, ne_up, ne_down, nh, gs_sector 

    type(basis_t) :: basis
    double precision :: eigval_tmp
    double precision, allocatable :: eigvec_tmp(:)

    gs%val = huge(1.d0)

    do isector=1,nsector
        ne_up = sectors(isector,1)
        ne_down = sectors(isector,2)
        nh = sectors(isector,3)

        if (master) then
            write(*,"(1x,a,I3)") "diagonalizing sector ",isector
            write(*,"(1x,a,I2,a2,I2,a2,I,a)") "(ne_up, ne_dn, dim) = (",&
                                ne_up,", ",ne_down,", ",nh,")" 
        endif
        
        call generate_basis(ne_up, ne_down, basis)
        
        ! DEBUG
        ! call dump_hamiltonian( em, ek, vmk, basis )
        ! stop

        allocate(eigvec_tmp(basis%nloc))
        
        ! find the ground state in the given sector
        t1_diag_loop = mpi_wtime(mpierr)    

        call mpi_barrier(comm, mpierr)
        select case(diag_solver)
            case (1)
                call diag_lanczos(em, ek, vmk, basis, eigval_tmp, eigvec_tmp)
            case (2)
                call diag_arpack(em, ek, vmk, basis, eigval_tmp, eigvec_tmp)
            case default
                call die("ground_state", "invalid diag solver.")
        end select

        if (eigval_tmp < gs%val) then
            gs%isector = isector
            gs%val = eigval_tmp
            if (allocated(gs%vec)) deallocate(gs%vec)
            allocate(gs%vec(basis%nloc))
            gs%vec = eigvec_tmp
        endif

        deallocate(eigvec_tmp)
        call dealloc_basis(basis)
        t2_diag_loop = mpi_wtime(mpierr)
    enddo

    if (master) then
        write(*,"(1x,a,I2,a2,I2,a2,I,a)") &
            "Ground state is at sector (ne_up, ne_dn, dim) = (",&
                            basis%ne_up,", ",basis%ne_down,", ",basis%ntot,")" 
    endif

    call mpi_barrier(comm, mpierr)
end subroutine ground_state
