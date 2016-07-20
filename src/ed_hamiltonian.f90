module ed_hamiltonian

    use io_units
    use mpi
    use dmft_lattice, only: hk
    use dmft_params, only: norb, U, Jex, mu, Up, Jp, nspin, tbham, na
    use utils, only: die
    use ed_params, only: nbath, kind_basis, nsite, ek_in, vk_in, nbathperorb, &
                         read_h_imp_params
    use ed_basis, only: basis_t, ed_basis_get,  ed_basis_idx
    implicit none

    public :: &
        generate_hamiltonian, &
        ed_hamiltonian_init,  &
        dump_hamiltonian_params,  &
        multiply_H

    double precision, allocatable, public :: &
        ek(:,:,:),    &   ! ek(nsite,2,na)         impurity/bath onsite energies
        vk(:,:,:,:)       ! vk(norb,nbath,2,na)    impurity-bath hybridization

    character(len=100), parameter :: FN_H_IMP_PARAMS = "h_imp.params"
    private
contains

    subroutine ed_hamiltonian_init
        integer :: ispin,iorb,i,j,ia
        logical :: found

        allocate(ek(nsite,2,na),vk(norb,nbath,2,na))

        found = .false.
        if (read_h_imp_params) then
            call import_h_imp_params(found)
        endif

        if (.not.found) then
            ! input levels, hybridization
            do ia=1,na
                do ispin=1,nspin
                    ek(:,ispin,ia)   = ek_in
                    vk(:,:,ispin,ia) = vk_in
                enddo
                if (nspin.eq.1) then
                    ek(:,2,ia) = ek_in
                    vk(:,:,2,ia) = vk_in
                endif
            enddo
        endif

        if (master) then
            write(*,*) 
            write(*,*) "Initial impurity/bath levels"
            do ia=1,na
                do ispin=1,nspin
                    write(*,*) 
                    write(*,"(a,I2,a,I2,a)") "(ia,ispin) = (",ia,",",ispin,")"
                    do i=1,norb
                        write(*,"(1x,A,I2,F12.6)") "orb ",i,ek(i,ispin,ia)
                    enddo
                    do i=norb+1,norb+nbath
                        write(*,"(1x,A,I2,F12.6)") "bath ",(i-norb),ek(i,ispin,ia)
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
                            write(*,"(F12.6)",advance="no") vk(j,i,ispin,ia)
                        enddo
                        write(*,*)
                    enddo
                enddo
            enddo
            write(*,*)
        endif
        call mpi_barrier(comm, mpierr)
    end subroutine ed_hamiltonian_init

    ! generates the Hamiltonian for the sector in basis
    subroutine generate_hamiltonian(ia, basis, H)
        integer, intent(in) :: ia
        type(basis_t), intent(in) :: basis
        double precision, intent(out) :: H(basis%nloc,basis%ntot)

        ! local variables
        integer(kind=kind_basis) :: basis_i, basis_j
        integer :: i,j,iorb,jorb,ibath,ispin, isite, jsite
        integer :: ni(2),nj(2),nb,sgntot, itmp

        if (.not.allocated(ek) .or. .not. allocated(vk)) then
            write(*,*) "Cluster Hamiltonian parameters are net initialized."
            call die("generate_hamiltonian", "ek, vk not allocated.")
        endif

        H = 0.0d0

        ! calculates <i|H|j>. 
        ! 1. find nonvanishing coefficients of H|j> = SUM_k |k> <k|H|j>
        ! 2. set H(i,j) = <i|H|j>
        jloop: do j=1,basis%nloc
            basis_j = ed_basis_get(basis,j)
            iorbloop: do iorb=1,norb
                if (BTEST(basis_j,iorb-1)) then
                    ni(1) = 1
                else
                    ni(1) = 0
                endif

                if (BTEST(basis_j,nsite+iorb-1)) then
                    ni(2) = 1 
                else
                    ni(2) = 0
                endif
                
                ! ==================================
                ! 1. onsite energies up/dn, 
                ! 2. intra orbital density-density
                ! ==================================
                H(j,j) = H(j,j) + (ek(iorb,1,ia)-mu)*ni(1) &
                            + (ek(iorb,2,ia)-mu)*ni(2) &
                            + U*ni(1)*ni(2)
                ! ==================================

                ! ==================================
                ! 3. inter orbital density-density
                ! ==================================
                do jorb=iorb+1,norb
                    ! n_j,up
                    if (BTEST(basis_j,jorb-1)) then
                        nj(1) = 1
                    else
                        nj(1) = 0
                    endif
                    ! n_j,dn
                    if (BTEST(basis_j,nsite+jorb-1)) then
                        nj(2) = 1
                    else
                        nj(2) = 0
                    endif

                    ! interorbital density-density
                    H(j,j) = H(j,j) + (Up-Jex)*(ni(1)*nj(1)+ni(2)*nj(2)) &
                                    + Up*(ni(1)*nj(2)+ni(2)*nj(1))
                enddo
                ! ==================================

                ! ==================================
                ! 4. hybridization
                ! ==================================
                do ispin=1,2
                    do ibath=1,nbath
                        isite = (ispin-1)*nsite+iorb
                        jsite = (ispin-1)*nsite+norb+ibath

                        if (ni(ispin).eq.0 .and. BTEST(basis_j,jsite-1)) then
                            ! c^+_{iorb,ispin} b_{ibath,ispin}
                            sgntot = sgn(basis_j,1,jsite-1)
                            basis_i = IBCLR(basis_j,jsite-1)

                            sgntot = sgntot*sgn(basis_i,1,isite-1)
                            basis_i = IBSET(basis_i,isite-1)

                            i = ed_basis_idx(basis, basis_i)

                            H(i,j) = H(i,j) + sgntot*vk(iorb,ibath,ispin,ia) 

                        else if (ni(ispin).ne.0 .and. &
                                 .not.BTEST(basis_j,jsite-1)) then
                            ! b^+_{ibath,ispin} c_{iorb,ispin}
                            sgntot = sgn(basis_j,1,isite-1)
                            basis_i = IBCLR(basis_j,isite-1)

                            sgntot = sgntot*sgn(basis_i,1,jsite-1)
                            basis_i = IBSET(basis_i,jsite-1)

                            i = ed_basis_idx(basis, basis_i)

                            H(i,j) = H(i,j) + sgntot*vk(iorb,ibath,ispin,ia) 

                        endif
                    enddo
                enddo
                ! ==================================

                ! ==================================
                ! 5. spin-flip & pair-hopping
                ! ==================================
                do jorb=1,norb
                    ! n_j,up
                    if (BTEST(basis_j,jorb-1)) then
                        nj(1) = 1
                    else
                        nj(1) = 0
                    endif
                    ! n_j,dn
                    if (BTEST(basis_j,nsite+jorb-1)) then
                        nj(2) = 1
                    else
                        nj(2) = 0
                    endif

                    ! spin-flip
                    ! -Jp c^+_{j,up} c_{j,dn} c^+_{i,dn} c_{i,up}
                    if (ni(1).ne.0.and.ni(2).eq.0.and.&
                        nj(1).eq.0.and.nj(2).ne.0) then

                        sgntot = sgn(basis_j,1,iorb-1)
                        basis_i = IBCLR(basis_j,iorb-1)

                        sgntot = sgntot*sgn(basis_i,1,nsite+iorb-1)
                        basis_i = IBSET(basis_i,nsite+iorb-1)

                        sgntot = sgntot*sgn(basis_i,1,nsite+jorb-1)
                        basis_i = IBCLR(basis_i,nsite+jorb-1)

                        sgntot = sgntot*sgn(basis_i,1,jorb-1)
                        basis_i = IBSET(basis_i,jorb-1)

                        i = ed_basis_idx(basis, basis_i)

                        H(i,j) = H(i,j) - sgntot*Jp

                    endif

                    ! pair-hopping
                    ! -Jp c^+_{j,up} c^+_{j,dn} c_{i,up} c_{i,dn}
                    if (ni(1).eq.0.and.ni(2).eq.0.and.&
                        nj(1).ne.0.and.nj(2).ne.0) then

                        sgntot = sgn(basis_j,1,nsite+iorb-1)
                        basis_i = IBCLR(basis_j,nsite+iorb-1)

                        sgntot = sgntot*sgn(basis_i,1,iorb-1)
                        basis_i = IBCLR(basis_i,iorb-1)

                        sgntot = sgntot*sgn(basis_i,1,nsite+jorb-1)
                        basis_i = IBSET(basis_i,nsite+jorb-1)

                        sgntot = sgntot*sgn(basis_i,1,jorb-1)
                        basis_i = IBSET(basis_i,jorb-1)

                        i = ed_basis_idx(basis, basis_i)

                        H(i,j) = H(i,j) - sgntot*Jp
                    endif
                enddo 
            enddo iorbloop

            ! ==================================
            ! 6. bath onsite
            ! ==================================
            do ibath=1,nbath
                if (BTEST(basis_j,norb+ibath-1)) then
                    H(j,j)=H(j,j)+ek(norb+ibath,1,ia)
                endif
                if (BTEST(basis_j,nsite+norb+ibath-1)) then
                    H(j,j)=H(j,j)+ek(norb+ibath,nspin,ia)
                endif
            enddo
        enddo jloop

    end subroutine generate_hamiltonian
    
    ! y = H*x
    ! matrix element is on-the-fly generated
    subroutine multiply_H(ia, basis, x, y, x_all)
        integer, intent(in) :: ia
        type(basis_t), intent(in) :: basis
        double precision, intent(in) :: x(basis%nloc)
        double precision, intent(out) :: x_all(basis%ntot)
        double precision, intent(out) :: y(basis%nloc)

        double precision :: rowsum
        integer(kind=kind_basis) :: bra,ket,a,b
        integer :: iorb, jorb, ispin, isite, jsite, j, jbath
        integer :: ni(2),nj(2),imin,imax
        logical :: jb

        call mpi_allgatherv(x, basis%nloc, mpi_double_precision, x_all,&
            basis%nlocals, basis%offsets, mpi_double_precision, comm, mpierr)

        y = 0.0d0

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
                rowsum = rowsum + &
                        ((ek(iorb,1,ia)-mu)*ni(1) + (ek(iorb,2,ia)-mu)*ni(2) &
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
                    do j=1,nbathperorb
                        jbath = (j-1)*norb+iorb
                        isite = (ispin-1)*nsite+iorb
                        jsite = (ispin-1)*nsite+norb+jbath

                        jb = BTEST(bra, jsite-1)

                        if (ni(ispin).eq.0 .and. jb) then
                            ! c^+_{iorb,ispin} b_{jbath,ispin}
                            ket = IBSET(IBCLR(bra,jsite-1), isite-1)
                            a = ed_basis_idx(basis, ket)
                            rowsum = rowsum &
                               +sgn(bra,isite,jsite)*vk(iorb,jbath,ispin,ia)&
                               *x_all(a)

                        else if (ni(ispin).ne.0 .and. .not.jb) then
                            ! b^+_{ibath,ispin} c_{iorb,ispin}
                            ket = IBSET(IBCLR(bra,isite-1),jsite-1)
                            a = ed_basis_idx(basis, ket)
                            rowsum = rowsum &
                               -sgn(bra,isite,jsite)*vk(iorb,jbath,ispin,ia)&
                               *x_all(a)
                        endif
                    enddo
                enddo

                ! ==================================
                ! 5. spin-flip & pair-hopping
                ! ==================================
                do jorb=1,norb
                    if (iorb==jorb) continue
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
                    rowsum = rowsum + ek(norb+jbath,1,ia)*x(b)
                endif
                if (BTEST(bra,nsite+norb+jbath-1)) then
                    rowsum = rowsum + ek(norb+jbath,nspin,ia)*x(b)
                endif
            enddo

            y(b) = rowsum
        enddo bloop

    end subroutine multiply_H

    ! (-1)^{ sum_{k=i}^{j-1} n_k }
    ! it is assumed that i<j
    integer function sgn(state, i, j) result(s)
        integer :: i, j, k
        integer(kind=kind_basis) :: state
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
        integer(kind=kind_basis) :: state
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

    subroutine dump_hamiltonian_params
        integer :: ia,iorb,ibath,ispin,i

        if (master) then
            open(unit=IO_H_PARAMS,file=FN_H_IMP_PARAMS,status="replace")
            ! na nspin norb nbath
            write(IO_H_PARAMS,*) na,2,norb,nbath

            do ia=1,na
                do ispin=1,2
                    do i=1,nsite
                        write(IO_H_PARAMS,*) ek(i,ispin,ia)
                    enddo 

                    do iorb=1,norb
                        do ibath=1,nbath
                            write(IO_H_PARAMS,*) vk(iorb,ibath,ispin,ia)
                        enddo
                    enddo
                enddo
            enddo

            close(IO_H_PARAMS)
        endif

    end subroutine dump_hamiltonian_params

    subroutine import_h_imp_params(found)
        logical, intent(out) :: found

        integer :: na2, dummy, norb2, nbath2
        integer :: ia,ibath,iorb,ispin,i

        if (master) then
            inquire(file=FN_H_IMP_PARAMS, exist=found)
            if (.not.found) then
                write(*,*) "warning: h_imp.params not found."
                return
            endif

            write(*,*) "Reading impurity Hamiltonian parameters from file."

            open(unit=IO_H_PARAMS,file=FN_H_IMP_PARAMS,status="old")
            ! na nspin norb nbath
            read(IO_H_PARAMS,*) na2,dummy,norb2,nbath2

            if (na2/=na.or.dummy/=2.or.norb2/=norb.or.nbath2/=nbath) then
                call die("import_h_imp_params", "dimension mismatch.")
                return
            endif

            do ia=1,na
                do ispin=1,2
                    do i=1,nsite
                        read(IO_H_PARAMS,*) ek(i,ispin,ia)
                    enddo 

                    do iorb=1,norb
                        do ibath=1,nbath
                            read(IO_H_PARAMS,*) vk(iorb,ibath,ispin,ia)
                        enddo
                    enddo
                enddo
            enddo
            close(IO_H_PARAMS)
        endif

        call mpi_bcast(ek, nsite*2*na, mpi_double_precision, 0, comm, mpierr)
        call mpi_bcast(vk, norb*nbath*2*na, mpi_double_precision, 0, comm, mpierr)
    end subroutine import_h_imp_params
end module ed_hamiltonian
