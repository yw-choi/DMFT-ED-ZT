subroutine dump_data( em, ek, vmk, D_cl, G_cl, G_loc )
    use dmft_params, only: norb, nbath, nspin, nw
    use dmft_grid, only: nwloc, omega, nw_procs, nw_offsets
    use io_units
    use mpi

    double precision, intent(in) :: &
        em(norb,2), &
        ek(nbath,2), &
        vmk(norb,nbath,2) 
    double complex, intent(in) :: &
        D_cl(nwloc,norb,nspin), &
        G_cl(nwloc,norb,nspin), &
        G_loc(nwloc,norb,nspin)

    double precision, allocatable :: &
        omega_all(:)

    double complex, allocatable :: &
        D_cl_all(:), &
        G_cl_all(:), &
        G_loc_all(:)

    integer :: iw, iorb, ispin, ibath

    character(len=100) :: fn

    ! dump em, ek, vmk
    if (master) then
        write(*,*) "Writing DMFT data ..."
        open(unit=IO_H_PARAMS, file=FN_H_PARAMS, status="replace")
        write(IO_H_PARAMS,*) norb, nbath, nspin
        do ispin=1,nspin
            do iorb=1,norb
                write(IO_H_PARAMS,*) em(iorb,ispin)
            enddo
            do ibath=1,nbath
                write(IO_H_PARAMS,*) ek(ibath,ispin)
            enddo
            do ibath=1,nbath
                do iorb=1,norb
                    write(IO_H_PARAMS,*) vmk(iorb,ibath,ispin)
                enddo
            enddo
        enddo
        close(IO_H_PARAMS)
    endif

    call mpi_barrier(comm, mpierr)

    ! dump D_cl, G_cl, G_loc
    allocate(omega_all(nw), D_cl_all(nw), G_cl_all(nw), G_loc_all(nw))
    call mpi_gatherv(omega, nwloc, mpi_double_precision, &
        omega_all, nw_procs, nw_offsets, mpi_double_precision, &
        0, comm, mpierr)

    do ispin=1,nspin
        do iorb=1,norb
            call mpi_gatherv(D_cl(:,iorb,ispin), nwloc, mpi_double_complex, &
                D_cl_all, nw_procs, nw_offsets, mpi_double_complex, &
                0, comm, mpierr)
            call mpi_gatherv(G_cl(:,iorb,ispin), nwloc, mpi_double_complex, &
                G_cl_all, nw_procs, nw_offsets, mpi_double_complex, &
                0, comm, mpierr)
            call mpi_gatherv(G_loc(:,iorb,ispin), nwloc, mpi_double_complex, &
                G_loc_all, nw_procs, nw_offsets, mpi_double_complex, &
                0, comm, mpierr)

            if (master) then
                write(fn,"(A10,A1,I1,A1,I1)") FN_GR, ".", iorb, ".", ispin
                open(unit=IO_GR, file=fn, status="replace")
                write(IO_GR,*) nw, iorb, ispin

                do iw=1,nw
                    write(IO_GR,"(7F12.8)") omega_all(iw), &
                        real(D_cl_all(iw)), &
                        aimag(D_cl_all(iw)), &
                        real(G_cl_all(iw)), &
                        aimag(G_cl_all(iw)), &
                        real(G_loc_all(iw)), &
                        aimag(G_loc_all(iw))
                enddo
                close(IO_GR)
            endif
        enddo
    enddo

    deallocate(omega_all, G_cl_all, D_cl_all, G_loc_all)

    call mpi_barrier(comm, mpierr)
end subroutine dump_data

