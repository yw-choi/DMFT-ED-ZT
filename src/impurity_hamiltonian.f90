module impurity_hamiltonian

    use mpi
    use dmft_params, only: nsite, norb, nbath, nspin, U, Up, Jex, Jp, mu, &
                           em_input, ek_input, vmk_input, init_h_imp, &
                           em_present, ek_present, vmk_present
    use dmft_grid, only: nwloc, omega
    use ed_basis, only: basis_t, ed_basis_get, ed_basis_idx
    use ed_projection, only: fit_h_imp_params 
    use utils, only: die

    implicit none

contains

    subroutine initial_h_imp_params(em, ek, vmk)
        double precision, intent(out) :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)

        integer :: ispin, iorb, ibath, i, j

        select case (init_h_imp)
            case (1)
                if (.not.em_present.or..not.ek_present.or..not.vmk_present) then
                    call die("initial_h_imp_params", "input is not given")
                endif
                em = em_input
                ek = ek_input
                vmk = vmk_input
            case (2)
                call import_h_imp_params( em, ek, vmk )
            case (3)
                call init_random( em, ek, vmk )
            case (4)
                call init_fit( em, ek, vmk )
            case default
                call die("initial_h_imp_params", "not implemented")
        end select

        if (master) then
            write(*,*) 
            write(*,*) "Initial impurity Hamiltonian parameters"
            do ispin=1,nspin
                write(*,*) 
                write(*,*) "Spin ",ispin
                write(*,*) "Impurity/Bath Levels"
                do i=1,norb
                    write(*,"(1x,A,I2,F12.6)") "orb  ",i,em(i,ispin)
                enddo
                do i=1,nbath
                    write(*,"(1x,A,I2,F12.6)") "bath ",i,ek(i,ispin)
                enddo
                write(*,*)
                write(*,*) "Impurity/Bath Hybridization"
                write(*,"(7x)", advance="no")
                do i=1,norb
                    write(*,"(4x,A3,I1,4x)",advance="no") "orb",i
                enddo
                write(*,*)
                do i=1,nbath
                    write(*,"(1x,a4,I2)",advance="no") "bath",i
                    do j=1,norb
                        write(*,"(F12.6)",advance="no") vmk(j,i,ispin)
                    enddo
                    write(*,*)
                enddo
            enddo
            write(*,*)
        endif

        call mpi_barrier(comm, mpierr)

    end subroutine initial_h_imp_params

    subroutine hx( em, ek, vmk, basis, x, y, x_all, py )

        type(basis_t), intent(in) :: basis
        double precision, intent(in) :: &
            x(basis%nloc), &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)
                                        
        logical, intent(in) :: py
        double precision, intent(out) :: y(basis%nloc), x_all(basis%ntot)

        double precision :: rowsum
        integer :: bra, ket, a, b
        integer :: iorb, jorb, ispin, isite, jsite, j, jbath
        integer :: ni(2), nj(2), imin, imax
        logical :: jb

        call mpi_allgatherv(x, basis%nloc, mpi_double_precision, x_all,&
            basis%nlocals, basis%offsets, mpi_double_precision, comm, mpierr)

        bloop: do b = 1, basis%nloc
            bra = ed_basis_get(basis, b)
            rowsum = 0.0d0

            do iorb = 1, norb
                if (BTEST(bra,iorb-1)) then
                    ni(1) = 1
                else
                    ni(1) = 0
                endif

                if (BTEST(bra,nsite+iorb-1)) then
                    ni(2) = 1 
                else
                    ni(2) = 0
                endif

                ! ==================================
                ! 1. onsite energies up/dn, 
                ! 2. intra orbital density-density
                ! ==================================
                rowsum = rowsum + ((em(iorb,1)-mu)*ni(1) &
                            + (em(iorb,2)-mu)*ni(2) &
                            + U*ni(1)*ni(2))*x(b)

                ! ==================================
                ! 3. inter orbital density-density
                ! ==================================
                do jorb = iorb+1, norb
                    ! n_j,up
                    if (BTEST(bra,jorb-1)) then
                        nj(1) = 1
                    else
                        nj(1) = 0
                    endif
                    ! n_j,dn
                    if (BTEST(bra,nsite+jorb-1)) then
                        nj(2) = 1
                    else
                        nj(2) = 0
                    endif
                    rowsum = rowsum + ((Up-Jex)*(ni(1)*nj(1)+ni(2)*nj(2)) &
                                    + Up*(ni(1)*nj(2)+ni(2)*nj(1)))*x(b)
                enddo

                ! ==================================
                ! 4. hybridization 
                ! ==================================
                do ispin=1,2
                    do jbath=1,nbath
                        if (abs(vmk(iorb,jbath,ispin))<1.d-16) then
                            cycle
                        endif

                        isite = (ispin-1)*nsite+iorb
                        jsite = (ispin-1)*nsite+norb+jbath

                        jb = BTEST(bra, jsite-1)

                        if (ni(ispin).eq.0 .and. jb) then
                            ! c^+_{iorb,ispin} b_{jbath,ispin}
                            ket = IBSET(IBCLR(bra,jsite-1), isite-1)
                            a = ed_basis_idx(basis, ket)
                            rowsum = rowsum &
                               +sgn(bra,isite,jsite)*vmk(iorb,jbath,ispin)&
                               *x_all(a)

                        else if (ni(ispin).ne.0 .and. .not.jb) then
                            ! b^+_{ibath,ispin} c_{iorb,ispin}
                            ket = IBSET(IBCLR(bra,isite-1),jsite-1)
                            a = ed_basis_idx(basis, ket)
                            rowsum = rowsum &
                               -sgn(bra,isite,jsite)*vmk(iorb,jbath,ispin)&
                               *x_all(a)
                        endif
                    enddo
                enddo

                ! ==================================
                ! 5. spin-flip & pair-hopping
                ! ==================================
                do jorb=1,norb
                    if (iorb==jorb) cycle

                    ! n_j,up
                    if (BTEST(bra,jorb-1)) then
                        nj(1) = 1
                    else
                        nj(1) = 0
                    endif

                    ! n_j,dn
                    if (BTEST(bra,nsite+jorb-1)) then
                        nj(2) = 1
                    else
                        nj(2) = 0
                    endif

                    imin = min(iorb,jorb)
                    imax = max(iorb,jorb)

                    ! spin-flip
                    ! Jp c^+_{i,up} c_{j,up} c^+_{j,dn} c_{i,dn} 
                    if (ni(1) == 0 .and. nj(1) /= 0 .and. &
                        nj(2) == 0 .and. ni(2) /= 0) then

                        ket = IBCLR(bra,nsite+iorb-1)
                        ket = IBSET(ket,nsite+jorb-1)
                        ket = IBCLR(ket,jorb-1)
                        ket = IBSET(ket,iorb-1)
                        a = ed_basis_idx(basis, ket)
                        rowsum = rowsum &
                            -sgn2(bra,imin,imax)*Jp*x_all(a)
                    endif

                    ! pair-hopping
                    ! Jp c^+_{i,up} c_{j,up} c^+_{i,dn} c_{j,dn}
                    if (ni(1).eq.0.and.ni(2).eq.0.and.&
                        nj(1).ne.0.and.nj(2).ne.0) then

                        ket = IBCLR(bra,nsite+jorb-1)
                        ket = IBSET(ket,nsite+iorb-1)
                        ket = IBCLR(ket,jorb-1)
                        ket = IBSET(ket,iorb-1)
                        a = ed_basis_idx(basis, ket)
                        rowsum = rowsum &
                            +sgn2(bra,imin,imax)*Jp*x_all(a)
                    endif
                enddo 
            enddo
            
            ! ==================================
            ! 6. bath onsite
            ! ==================================
            do jbath=1,nbath
                if (BTEST(bra,norb+jbath-1)) then
                    rowsum = rowsum + ek(jbath,1)*x(b)
                endif
                if (BTEST(bra,nsite+norb+jbath-1)) then
                    rowsum = rowsum + ek(jbath,2)*x(b)
                endif
            enddo

            if (py) then
                y(b) = y(b) + rowsum
            else
                y(b) = rowsum
            endif
        enddo bloop
    end subroutine hx

    ! (-1)^{ sum_{k=i}^{j-1} n_k }
    ! it is assumed that i<j
    integer function sgn(state, i, j) result(s)
        integer :: i, j, k
        integer :: state
        s = 1
        do k=i,j-1 
            if (BTEST(state,k-1)) then
                s = -s
            endif
        enddo
    end function sgn

    ! (-1)^{ sum_{k=i}^{j-1} n_{k,up}+n_{k,dn} }
    ! it is assumed that i<j, i,j <= Nsite
    integer function sgn2(state, i, j) result(s)
        integer :: i, j, k
        integer :: state
        s = 1
        do k=i,j-1 
            if (BTEST(state,k-1)) then
                s = -s
            endif
            if (BTEST(state,nsite+k-1)) then
                s = -s
            endif
        enddo
    end function sgn2

    subroutine import_h_imp_params( em, ek, vmk )
        use io_units
        double precision, intent(out) :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)

        integer :: norb2, nbath2, nspin2, ispin, iorb, ibath

        logical :: found

        if (master) then
            inquire(file=FN_H_PARAMS, exist=found)
            if (.not.found) then
                call die("import_h_imp_params", "h_imp.dat not found")
            endif
            open(unit=IO_H_PARAMS,file=FN_H_PARAMS,status="old")
            read(IO_H_PARAMS,*) norb2, nbath2, nspin2
            if (norb2/=norb.or.nbath2/=nbath.or.nspin2/=nspin) then
                call die("import_h_imp_params", "dimension mismatch.")
            endif
            do ispin=1,nspin
                do iorb=1,norb
                    read(IO_H_PARAMS,*) em(iorb,ispin)
                enddo
                do ibath=1,nbath
                    read(IO_H_PARAMS,*) ek(ibath,ispin)
                enddo
                do ibath=1,nbath
                    do iorb=1,norb
                        read(IO_H_PARAMS,*) vmk(iorb,ibath,ispin)
                    enddo
                enddo
            enddo

            if (nspin==1) then
                em(:,2) = em(:,1)
                ek(:,2) = ek(:,1)
                vmk(:,:,2) = vmk(:,:,1)
            endif

            close(IO_H_PARAMS)
        endif
        call mpi_bcast(em, norb*2, mpi_double_precision, 0, comm, mpierr)
        call mpi_bcast(ek, nbath*2, mpi_double_precision, 0, comm, mpierr)
        call mpi_bcast(vmk, norb*nbath*2, mpi_double_precision, 0, comm, mpierr)
    end subroutine import_h_imp_params

    subroutine init_random( em, ek, vmk )
        double precision, intent(out) :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)
        integer :: ispin, iorb, ibath

        double precision :: r

        if (master) then
            call random_seed
            do ispin=1,2
                do iorb=1,norb
                    call random_number(r)
                    em(iorb,ispin) = r*2-1
                enddo
                do ibath=1,nbath
                    call random_number(r)
                    ek(ibath,ispin) = r*2-1
                enddo
                do ibath=1,nbath
                    do iorb=1,norb
                        call random_number(r)
                        vmk(iorb,ibath,ispin) = r*2-1
                    enddo
                enddo
            enddo
        endif
        
        call mpi_bcast(em, norb*2, mpi_double_precision, 0, comm, mpierr)
        call mpi_bcast(ek, nbath*2, mpi_double_precision, 0, comm, mpierr)
        call mpi_bcast(vmk, norb*nbath*2, mpi_double_precision, 0, comm, mpierr)
    end subroutine init_random
    
    subroutine init_fit( em, ek, vmk )
        double precision, intent(out) :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)

        integer :: iw
        double complex, allocatable :: &
            G_cl(:,:,:),      & ! G_cl(nwloc,norb,nspin) 
            G_loc(:,:,:),     & ! G_loc(nwloc,norb,nspin)
            D_cl(:,:,:)         ! D_cl(nwloc,norb,nspin)
        allocate(G_cl(nwloc,norb,nspin), G_loc(nwloc,norb,nspin))
        allocate(D_cl(nwloc,norb,nspin))

        do iw=1,nwloc
            G_cl(iw,:,:) = 1/cmplx(0.d0,omega(iw))
            G_loc(iw,:,:) = 1/cmplx(0.d0,omega(iw))
        enddo
        D_cl = 0.25d0*G_cl

        call init_random( em, ek, vmk )

        ! find new em, ek, vmk 
        call fit_h_imp_params( em, ek, vmk, D_cl, G_cl, G_loc )
    end subroutine init_fit
end module impurity_hamiltonian
