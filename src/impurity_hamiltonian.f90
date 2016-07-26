module impurity_hamiltonian

    use mpi
    use dmft_params, only: nsite, norb, nbath, nspin, U, Up, Jex, Jp, mu, &
                           em_input, ek_input, vmk_input

    use ed_basis, only: basis_t, ed_basis_get, ed_basis_idx

    implicit none

contains

    subroutine initial_h_imp_params(em, ek, vmk)
        double precision, intent(out) :: &
            em(norb,2), &
            ek(nbath,2), &
            vmk(norb,nbath,2)

        integer :: ispin, iorb, ibath

        ! @TODO read h_imp params from the save file

        em = em_input
        ek = ek_input
        vmk = vmk_input

        if (master) then
            write(*,*) "Initial mpurity Hamiltonian parameters"
            do ispin=1,nspin
                write(*,*) "spin ", ispin

                write(*,*) "impurity levels"
                do iorb=1,norb
                    write(*,"(F6.3)") em(iorb,ispin)
                enddo

                write(*,*) "bath levels"
                do ibath=1,nbath
                    write(*,"(F6.3)") ek(ibath,ispin)
                enddo

                write(*,*) "hybridization"
                do ibath=1,nbath
                    do iorb=1,norb
                        write(*,"(F6.3)",advance="no") vmk(iorb,ibath,ispin)
                    enddo
                    write(*,*)
                enddo
                write(*,*)
            enddo
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
end module impurity_hamiltonian
