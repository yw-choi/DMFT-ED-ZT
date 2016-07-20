module ed_solver
    use mpi
    use utils, only: die
    use dmft_params, only: norb, nspin, mu, na
    use dmft_grid, only: nwloc, omega
    use ed_projection, only: project_to_impurity_model
    use ed_params, only: ed_read_params, nbath, nsite, nsector, &
                         nev, nstep, sectors, ek_in, vk_in, diag_method

    use ed_hamiltonian, only: ed_hamiltonian_init, ek,vk, dump_hamiltonian_params
    use ed_green, only: ed_green_init, cluster_green_ftn, G_cl, g_coeffs 
    use ed_diag_arpack, only: diag_arpack
    use timer, only: t1_green_loop, t2_green_loop, print_elapsed_time
    use ed_eigpair, only: allocate_eigpairs
    use io_units

    implicit none

contains

    subroutine ed_init
        call ed_read_params
        call ed_hamiltonian_init
        call ed_green_init

        call allocate_eigpairs(nev)
    end subroutine ed_init

    subroutine ed_solve(iloop,ia,G0,Sigma)
        integer, intent(in) :: iloop, ia
        double complex, intent(in) :: G0(nwloc,norb,nspin)
        double complex, intent(out) :: Sigma(nwloc,norb,nspin)
        integer :: iw,iorb,ispin

        ! Find ek,vk by fitting the cluster quantity(e.g. G0cl) to G0.
        if (iloop>1) then
            call project_to_impurity_model(G0, ek(:,:,ia), &
                                           vk(:,:,:,ia))
        endif

        call dump_hamiltonian_params

        ! Diagonalize the AIM Hamiltonian characterized by ek,vk and
        ! return the eigpairs.
        select case(diag_method)
            ! case ("full")
            !     call diag_full(ia)
            case ("arpack")
                call diag_arpack(ia)
            case default
                call die("ed_solve", "Diagonalization method not implemented.")
        end select

        t1_green_loop = mpi_wtime(mpierr)
        call cluster_green_ftn(ia)
        t2_green_loop = mpi_wtime(mpierr)

        if (master) then
            call print_elapsed_time(" Cluster Green's function calculation time",&
                t1_green_loop,t2_green_loop)
        endif

        ! the self-energy 
        call cluster_self_energy(ia,Sigma)

    end subroutine ed_solve

    subroutine cluster_self_energy(ia, Sigma)
        integer, intent(in) :: ia
        double complex, intent(out) :: Sigma(nwloc,norb,nspin)
        double complex :: sig
        integer :: iw,iorb,ispin

        do ispin=1,nspin
            do iorb=1,norb
                do iw=1,nwloc
                    Sigma(iw,iorb,ispin) = cmplx(0.0d0,omega(iw))+mu &
                        -ek(iorb,ispin,ia)-delta_cl(ia,iw,iorb,ispin)&
                        -1.d0/G_cl(iw,iorb,ispin)
                enddo
            enddo
        enddo
    end subroutine cluster_self_energy

    double complex function delta_cl(ia,iw,iorb,ispin)
        integer :: iw,iorb,ibath,ispin,ia

        delta_cl = cmplx(0.0d0,0.0d0)

        do ibath=1,nbath
            delta_cl = delta_cl + vk(iorb,ibath,ispin,ia)*vk(iorb,ibath,ispin,ia)&
                        /(cmplx(0.0d0,omega(iw))-ek(norb+ibath,ispin,ia))
        enddo
    end function delta_cl

    subroutine ed_post_processing

        integer :: ia,ispin,iorb,istep,iev

        character(len=100) :: fn

        ! dump green function coefficients
        if (master) then
            write(*,*) "Writing cluster Green's function coefficients..."
            open(unit=IO_G_COEFFS, file="g_coeffs.dat", status="replace")
            write(IO_G_COEFFS,*) na,2,norb
            do ia=1,na
                do ispin=1,2
                    do iorb=1,norb
                        write(IO_G_COEFFS,"(2I5)") &
                            g_coeffs(iorb,ispin,ia)%nstep,&
                            g_coeffs(iorb,ispin,ia)%nev
                        do iev=1,g_coeffs(iorb,ispin,ia)%nev
                            write(IO_G_COEFFS,"(L,1x,4I3,2E)") &
                                g_coeffs(iorb,ispin,ia)%even,&
                                ia,ispin,iorb,iev, &
                                g_coeffs(iorb,ispin,ia)%eigval(iev), &
                                g_coeffs(iorb,ispin,ia)%prob(iev)

                            write(IO_G_COEFFS,"(4E)") &
                                (g_coeffs(iorb,ispin,ia)%ap(istep,iev), &
                                 g_coeffs(iorb,ispin,ia)%bp(istep,iev), &
                                 g_coeffs(iorb,ispin,ia)%an(istep,iev), &
                                 g_coeffs(iorb,ispin,ia)%bn(istep,iev), &
                                 istep=1,g_coeffs(iorb,ispin,ia)%nstep)
                        enddo
                    enddo
                enddo
            enddo
            close(IO_G_COEFFS)
        endif
        call mpi_barrier(comm,mpierr)

    end subroutine ed_post_processing
end module ed_solver
