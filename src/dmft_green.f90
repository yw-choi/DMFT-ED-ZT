module dmft_green
    use dmft_params
    use dmft_lattice, only: dos, dosw, hk
    use mpi
    use dmft_grid
    use io_units
    use numeric_utils, only: cinv
    use alloc, only: alloc_count

    implicit none

    public :: &
        dmft_green_init, &
        update_local_green_ftn, &
        update_weiss_ftn, &
        dump_green, &
        find_density

    double complex, allocatable, public :: &
        G_prev(:,:,:),   & ! G_prev(nwloc,norb,nspin)
                           ! local Green's function of the previous step
        G(:,:,:),        & ! G(nwloc,norb,nspin)
                           ! local Green's function of the current step
        G0(:,:,:),       & ! G0(nwloc,norb,nspin)
                           ! Weiss field
        Sigma(:,:,:)       ! Sigma(nwloc,norb,nspin)
                           ! self energy

    double precision, allocatable, public :: &
        density(:,:)       ! density(norb,nspin)  particle density

    character(len=100), parameter :: FN_GR_SAVE = "green.save"

    private 
contains

    subroutine dmft_green_init
        logical :: found

        allocate(G_prev(nwloc,norb,nspin))
        allocate(G(nwloc,norb,nspin))
        allocate(G0(nwloc,norb,nspin))
        allocate(Sigma(nwloc,norb,nspin))

        allocate(density(norb,nspin,na))

        G_prev = cmplx(0.0d0,0.0d0)
        G      = cmplx(0.0d0,0.0d0)
        G0     = cmplx(0.0d0,0.0d0)
        Sigma  = cmplx(0.0d0,0.0d0)

        found = .false.
        if (read_gf_from_file) then
            call import_green_ftn(found)
        endif
        
        if (.not.found) then
            if (master) then
                write(*,*) "Initial Green's function from uncorrelated lattice."
            endif
            ! initial green function and weiss field
            ! from uncorrelated lattice green function
            call update_local_green_ftn
            call update_weiss_ftn
        end if

    end subroutine dmft_green_init

    ! Calculates new local green function from the given self-energy,
    ! by summing the lattice green function over k.
    subroutine update_local_green_ftn
        integer :: ik, iw, ispin, ia, iorb
        ! Lattice Green's function at (ispin,iw,ik)
        double complex :: gk(norb,nspin,na), gf

        if (master) then
            write(*,*) "Updating the local Green's function..."
        endif

        ! Keep the previous Local Green Function
        G_prev = G

        ! calculate new local Green's function from the lattice Green's function
        G = cmplx(0.0d0,0.0d0)
        if (tbham==0) then
            ! general tight-binding hamiltonian.
            ! obtain the lattice green's function through matrix inversion.
            do iw=1,nwloc
                do ik=1,nk
                    call gkdiag(iw, ik, sigma, gk)
                    G(iw,:,:,:) = G(iw,:,:,:) + gk(:,:,:)
                enddo
            enddo
            G = G / nk
        else if (tbham==1) then
            ! degenerate bands on 2d square lattice
            ! tight-binding hamiltonian is diagonal. no need of matrix inv.
            do iw=1,nwloc
                do ia=1,na
                    do ispin=1,nspin
                        do iorb=1,norb
                            do ik=1,nk
                                gf = cmplx(0.d0,omega(iw))+mu &
                                    -hk(ik,iorb,iorb,ispin,ia,ia) &
                                    -sigma(iw,iorb,ispin,ia)
                                gf = 1.d0/gf
                                G(iw,iorb,ispin,ia) = G(iw,iorb,ispin,ia)+gf
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            G = G / nk
        else if (tbham==2) then
            ! Bethe lattice.
            ! rather than sum over k, integrate over E using dos
            do iw=1,nwloc
                do ispin=1,nspin
                    do iorb=1,norb
                        call green_bethe(iw,ispin,iorb,gf) 
                        G(iw,iorb,ispin,1) = gf
                    enddo
                enddo
            enddo
        else
            ! other cases should be prevented in the initialization step
            stop "This should not happen"
        endif

    end subroutine update_local_green_ftn

    subroutine update_weiss_ftn

        if (master) then
            write(*,*) "Updating the Weiss field..."
        endif

        G0 = 1.d0/(1.d0/G+Sigma)

    end subroutine update_weiss_ftn

    subroutine dump_green
        ! dump current green ftn, weiss ftn, self energy
        integer :: iw, iorb, ia, ispin
        character(len=100) :: fn
        double precision :: zq, omega_all(nw)
        double complex :: G_all(nw), G0_all(nw), Sigma_all(nw)

        if (master) then
            open(unit=IO_GR_DATA, file=FN_GR_SAVE, status="replace")
            write(IO_GR_DATA, *) na,nspin,norb,nw
        endif

        do ia=1,na
            do ispin=1,nspin
                do iorb=1,norb
                    call mpi_allgatherv(omega(:), &
                                        nwloc, &
                                        mpi_double_precision, &
                                        omega_all, &
                                        nw_procs, &
                                        nw_offsets, &
                                        mpi_double_precision, &
                                        comm, &
                                        mpierr)
                    call mpi_allgatherv(G(:,iorb,ispin,ia), &
                                        nwloc, &
                                        mpi_double_complex, &
                                        G_all, &
                                        nw_procs, &
                                        nw_offsets, &
                                        mpi_double_complex, &
                                        comm, &
                                        mpierr)
                    call mpi_allgatherv(G0(:,iorb,ispin,ia), &
                                        nwloc, &
                                        mpi_double_complex, &
                                        G0_all, &
                                        nw_procs, &
                                        nw_offsets, &
                                        mpi_double_complex, &
                                        comm, &
                                        mpierr)
                    call mpi_allgatherv(Sigma(:,iorb,ispin,ia), &
                                        nwloc, &
                                        mpi_double_complex, &
                                        Sigma_all, &
                                        nw_procs, &
                                        nw_offsets, &
                                        mpi_double_complex, &
                                        comm, &
                                        mpierr)

                    if (master) then
                        write(IO_GR_DATA,*) ia, ispin, iorb
                        do iw=1,nw
                            write(IO_GR_DATA,"(7E20.10)") omega_all(iw), &
                                real(G_all(iw)), &
                                aimag(G_all(iw)), &
                                real(G0_all(iw)), &
                                aimag(G0_all(iw)), &
                                real(Sigma_all(iw)), &
                                aimag(Sigma_all(iw))
                        enddo
                    endif

                    call mpi_barrier(comm,mpierr)
                enddo
            enddo
        enddo

        if (master) then
            close(IO_GR_DATA)
        endif

        call mpi_barrier(comm,mpierr)
    end subroutine dump_green

    ! diagonal element of G_lat(iw,ik)
    subroutine gkdiag(iw,ik,sigma,gk)
        integer, intent(in) :: iw, ik
        double complex, intent(in) :: sigma(nwloc, norb, nspin, na)
        double complex, intent(out) :: gk(norb, nspin, na)

        double complex :: invGk(na*norb*nspin, na*norb*nspin),&
            Gkmat(na*norb*nspin,na*norb*nspin)
        integer :: i,j, ia,iorb, ja,jorb, ispin

        ! invGk = iw + mu - Hk - Sigma
        invGk = cmplx(0.0d0,0.0d0)

        do ia=1,na
            do ispin=1,nspin
                do iorb=1,norb
                    i = (ia-1)*norb*nspin+(ispin-1)*norb+iorb
                    invGk(i,i) = invGk(i,i) + &
                        cmplx(0.0d0,omega(iw)) + mu - Sigma(iw,iorb,ispin,ia)

                    do ja=1,na
                        do jorb=1,norb
                            j = (ja-1)*norb*nspin+(ispin-1)*norb+jorb
                            invGk(j,i) = invGk(j,i) - Hk(ik,iorb,jorb,ispin,ia,ja)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! matrix inversion
        call cinv(invGk, na*norb*nspin, na*norb*nspin, Gkmat)

        ! Extract diagonal elements to output
        do ia=1,na
            do ispin=1,nspin
                do iorb=1,norb
                    i = (ia-1)*norb*nspin+(ispin-1)*norb+iorb
                    gk(iorb,ispin,ia) = Gkmat(i,i)
                enddo
            enddo
        enddo
    end subroutine gkdiag
    
    subroutine green_bethe(iw,ispin,iorb,gf)
        integer, intent(in) :: iw, ispin, iorb
        double complex, intent(out) :: gf
        integer :: ie
        double precision :: e

        gf = cmplx(0.d0,0.d0)
        ! trapezoidal rule integration
        ! note that at E=-1,1, DOS(E) is exactly 0 so that
        ! we do not care about the boundary term
        do ie=2,nwmodel-1
            e = dosw(ie)

            gf = gf + dos(ie,iorb,ispin,1)/(cmplx(0.d0,omega(iw))+mu-e &
                                    -sigma(iw,iorb,ispin,1))
        enddo

        gf = gf*2.d0/(nwmodel-1)

    end subroutine green_bethe

    subroutine import_green_ftn(found)
        use utils, only: die
        logical :: found
        integer :: iw,na2,nspin2,norb2,nw2,&
                   ia2,ispin2,iorb2,iorb,ispin,ia
        double complex :: G_all(nw), G0_all(nw), Sigma_all(nw)
        double precision :: w, reg,img, reg0,img0, resig,imsig

        if (master) then
            inquire (file=FN_GR_SAVE, exist=found)
            if (.not.found) then
                write(*,*) "warning: Green's function save file not found."
                return
            endif
            open(unit=IO_GR_DATA, file=FN_GR_SAVE, status="old")
            read(IO_GR_DATA,*) na2,nspin2,norb2,nw2

            write(*,*) "Found Green's function save file."
            write(*,*) "Dimensions : "
            write(*,"(3x,a,I5)") "na    = ", na2
            write(*,"(3x,a,I5)") "nspin = ", nspin2
            write(*,"(3x,a,I5)") "norb  = ", norb2
            write(*,"(3x,a,I5)") "nw    = ", nw2

            ! dimension matching
            ! @TODO reject different beta ?
            if (na2/=na .or. nspin2/=nspin .or. norb2/=norb .or. &
                nw2/=nw) then
                call die("import_green_ftn", &
                    "Dimension of the saved Green's function does not match")
                return
            endif
        endif

        do ia=1,na
            do ispin=1,nspin
                do iorb=1,norb
                    if (master) then
                        read(IO_GR_DATA, *) ia2, ispin2, iorb2

                        do iw=1,nw
                            read(IO_GR_DATA,"(7E20.10)") &
                                w, reg, img, reg0, img0, resig, imsig

                            G_all(iw) = cmplx(reg,img)
                            G0_all(iw) = cmplx(reg0,img0)
                            Sigma_all(iw) = cmplx(resig,imsig)
                        enddo
                    endif

                    call mpi_scatterv(G_all, nw_procs, nw_offsets, &
                        mpi_double_complex, G(:,iorb,ispin,ia), nwloc, &
                        mpi_double_complex, 0, comm, mpierr)

                    call mpi_scatterv(G0_all, nw_procs, nw_offsets, &
                        mpi_double_complex, G0(:,iorb,ispin,ia), nwloc, &
                        mpi_double_complex, 0, comm, mpierr)

                    call mpi_scatterv(Sigma_all, nw_procs, nw_offsets, &
                        mpi_double_complex, Sigma(:,iorb,ispin,ia), nwloc, &
                        mpi_double_complex, 0, comm, mpierr)

                enddo
            enddo
        enddo

        if (master) then
            write(*,*) "successfully imported ", FN_GR_SAVE
            close(IO_GR_DATA)
        endif
    end subroutine import_green_ftn

    subroutine find_density

        integer :: iorb,ispin,ia

        double precision, allocatable :: density_loc(:,:,:)
        double precision :: nocc
        allocate(density_loc(norb,nspin,na))

        density_loc = 0.d0
        do ia=1,na
            do ispin=1,nspin
                do iorb=1,norb
                    density_loc(iorb,ispin,ia) = sum(real(G(:,iorb,ispin,ia)))
                enddo
            enddo
        enddo

        call mpi_allreduce(density_loc,density,norb*nspin*na,&
            mpi_double_precision,mpi_sum,comm,mpierr)

        density = 2.d0*density/beta + 0.5d0

        if (nspin==1) then
            density = density*2.d0
        endif

        if (master) then
            write(*,*) 
            write(*,*) "Particle occupancy"
            do ia=1,na
                do ispin=1,nspin
                    do iorb=1,norb
                        write(*,"(1x,A,I2,A,I2,A,I2,A,F6.3)") &
                            "n(",iorb,",",ispin,",",ia,") = ", &
                            density(iorb,ispin,ia)
                    enddo
                enddo
            enddo
            write(*,"(1x,A,F6.3)") "Total       = ", sum(density)
        endif

        deallocate(density_loc)
    end subroutine find_density
end module dmft_green
