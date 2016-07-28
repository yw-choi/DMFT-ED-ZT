subroutine post_process( em, ek, vmk, D_cl, G_cl, G_loc, gs )
    use mpi
    use dmft_params, only: nspin, norb, nbath, nsite, nw, sectors
    use dmft_grid, only: nwloc, omega
    use eigpair, only: eigpair_t
    use ed_basis, only: basis_t, ed_basis_get, dealloc_basis, generate_basis

    double precision, intent(in) :: &
        em(norb,2),  &
        ek(nbath,2), &
        vmk(norb,nbath,2)

    double complex, intent(in) :: &
        D_cl(nwloc,norb,nspin), &
        G_cl(nwloc,norb,nspin), &
        G_loc(nwloc,norb,nspin)

    type(eigpair_t), intent(in) :: gs

    type(basis_t) :: basis

    double complex, allocatable :: &
        Sigma(:,:,:)

    double precision, allocatable :: &
        occ_loc(:,:), &
        occ(:,:)

    double precision, allocatable :: &
        docc_loc(:), &
        docc(:)
    double precision :: zq

    integer :: i, n_ms, ib

    allocate(Sigma(nwloc,norb,nspin))
    allocate(occ_loc(norb,nspin), occ(norb,nspin))
    allocate(docc_loc(norb), docc(norb))

    call dump_sigma( em, D_cl, G_cl, Sigma )

    ! approximate Quasiparticle weight 
    if (master) then
        write(*,*) 
        write(*,*) "Quasiparticle weight"
        write(*,*) 
        write(*,*) "    Z(iorb,ispin)"
        write(*,*) "    ---------------------------"
        do ispin=1,nspin
            do iorb=1,norb
                zq = aimag(sigma(1,iorb,ispin))/omega(1)
                zq = 1.d0/(1.d0-zq)
                write(*,"(a,I2,a,I2,a,F12.6)") &
                    "     Z(",iorb,",",ispin,") = ",zq
            enddo
        enddo
        write(*,*) 
    endif

    call mpi_barrier(comm, mpierr)

    ! Occupancy
    call generate_basis( sectors(gs%isector,1), sectors(gs%isector,2), basis ) 
    do ispin=1,nspin
        do iorb=1,norb
            occ_loc(iorb,ispin) = 0.d0
            do i=1,basis%nloc
                ib = ed_basis_get(basis, i)
                if (BTEST(ib,(ispin-1)*nsite+iorb-1)) then
                    occ_loc(iorb,ispin) = occ_loc(iorb,ispin) &
                                              + gs%vec(i)*gs%vec(i)
                endif
            enddo
        enddo
    enddo

    call mpi_allreduce(occ_loc,occ,norb*nspin,&
        mpi_double_precision,mpi_sum,comm,mpierr)

    if (master) then
        write(*,*) 
        write(*,*) "Particle occupancy"
        do ispin=1,nspin
            do iorb=1,norb
                write(*,"(1x,A,I2,A,I2,A,F6.3)") &
                    "n(",iorb,",",ispin,") = ", &
                    occ(iorb,ispin)
            enddo
        enddo
        write(*,"(1x,A,F6.3)") "Total    = ", sum(occ)
    endif

    ! double occupancy
    do iorb=1,norb
        docc_loc(iorb) = 0.d0
        do i=1,basis%nloc
            ib = ed_basis_get(basis, i)
            if (BTEST(ib,iorb-1).and.BTEST(ib,nsite+iorb-1)) then
                docc_loc(iorb) = docc_loc(iorb) + gs%vec(i)*gs%vec(i)
            endif
        enddo
    enddo
    call mpi_allreduce(docc_loc,docc,norb,&
        mpi_double_precision,mpi_sum,comm,mpierr)
    if (master) then
        write(*,*) 
        write(*,*) "Double occupancy"
        do iorb=1,norb
            write(*,"(1x,A,I2,A,F6.3)") &
                "docc(",iorb,") =  ", docc(iorb)
        enddo
    endif

    call dealloc_basis(basis)
    deallocate(occ_loc, occ)
    deallocate(docc_loc, docc)
    deallocate(Sigma)
    call mpi_barrier(comm, mpierr)
end subroutine

subroutine dump_sigma( em, D_cl, G_cl, Sigma )
    use dmft_params, only: norb, nbath, nspin, nw, mu
    use dmft_grid, only: nwloc, omega, nw_procs, nw_offsets
    use io_units
    use mpi

    double precision, intent(in) :: &
        em(norb,2)
    double complex, intent(in) :: &
        D_cl(nwloc,norb,nspin), &
        G_cl(nwloc,norb,nspin)
    double complex, intent(out) :: &
        Sigma(nwloc,norb,nspin)

    double precision, allocatable :: &
        omega_all(:)

    double complex, allocatable :: &
        Sigma_all(:)

    integer :: iw, iorb, ispin, ibath

    character(len=100) :: fn

    ! dump em, ek, vmk
    if (master) then
        write(*,*) "Writing the self-energy ..."
    endif

    ! dump D_cl, G_cl, G_loc
    allocate(omega_all(nw), Sigma_all(nw))
    call mpi_gatherv(omega, nwloc, mpi_double_precision, &
        omega_all, nw_procs, nw_offsets, mpi_double_precision, &
        0, comm, mpierr)

    do ispin=1,nspin
        do iorb=1,norb
            Sigma = 0.d0
            do iw=1,nwloc
                Sigma(iw,iorb,ispin) = cmplx(0.d0,omega(iw))+mu&
                    -em(iorb,ispin)-D_cl(iw,iorb,ispin)-1.d0/G_cl(iw,iorb,ispin)
            enddo
            call mpi_gatherv(Sigma(:,iorb,ispin), nwloc, mpi_double_complex, &
                Sigma_all, nw_procs, nw_offsets, mpi_double_complex, &
                0, comm, mpierr)

            if (master) then
                write(fn,"(A9,A1,I1,A1,I1)") FN_SIGMA, ".", iorb, ".", ispin
                open(unit=IO_SIGMA, file=fn, status="replace")
                write(IO_SIGMA,*) nw, iorb, ispin

                do iw=1,nw
                    write(IO_SIGMA,"(7F12.8)") omega_all(iw), &
                        real(Sigma_all(iw)), &
                        aimag(Sigma_all(iw))
                enddo
                close(IO_SIGMA)
            endif
        enddo
    enddo

    deallocate(omega_all, Sigma_all)

    call mpi_barrier(comm, mpierr)
end subroutine dump_sigma
