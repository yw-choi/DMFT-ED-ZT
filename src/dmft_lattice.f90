module dmft_lattice
    use mpi
    use constants
    use dmft_grid
    use dmft_params
    use io_units 
    use utils, only: die

    implicit none

    public :: &
        dmft_lattice_init

    double complex, allocatable, public :: &
        Hk(:,:,:,:)    ! Hk(nk,norb,norb,nspin) tight-binding hamiltonian

    double precision, allocatable, public :: &
        dos(:,:,:),  & ! dos(ne,norb,nspin) tight-binding dos
        egrid(:),    & ! egrid(ne) dos integration energy grid
        occ0(:,:)      ! occ0(norb,nspin) tight-binding occupancy

    integer, parameter, public :: ne = 2000
    double precision, parameter, public :: &
        emin = -1.d0, &
        emax =  1.d0

    private
contains

    subroutine dmft_lattice_init
        integer :: ie
        double precision :: de

        if (master) then
            write(*,*) "Setting up the tight-binding Hamiltonian..."
        endif

        ! dos energy grid (used only if dos integration is needed)
        allocate(dos(ne,norb,nspin))
        allocate(occ0(norb,nspin))
        allocate(egrid(ne))

        egrid = 0.d0
        de = (emax-emin)/(ne-1)
        do ie=1,ne
            egrid(ie) = emin + de*(ie-1)
        enddo

        select case (tbham)
            case (0)
                call read_tb_hamiltonian
            case (1)
                call tbh_model_1
            case (2)
                call tbh_model_2
            case default
                call die("dmft_lattice_init", &
                    "invalid tight-binding hamiltonian type")
        end select

        call dump_dos
    end subroutine dmft_lattice_init

    subroutine read_tb_hamiltonian
        ! @TODO
        call die("read_tb_hamiltonian", "not implemented yet")

    end subroutine read_tb_hamiltonian

    subroutine tbh_model_1
        ! @TODO
        call die("read_tb_hamiltonian", "not implemented yet")
    end subroutine tbh_model_1

    subroutine tbh_model_2
        ! all degenerate Bethe lattice, infinite coordination
        ! half bandwidth = 1

        integer :: iorb,ispin,ie

        double precision :: e, de

        dos = 0.d0
        do ispin=1,nspin
            do iorb=1,norb
                do ie=1,ne
                    e = egrid(ie)
                    dos(ie,iorb,ispin) = 2.d0/pi*sqrt(max(1.d0-e*e, 0.d0))
                enddo 
            enddo
        enddo

        do ispin=1,nspin
            do iorb=1,norb
                ! integration by trapezoidal rule
                ! boundary values are exactly 0 in this case.
                occ0(iorb,ispin) = sum(dos(:,iorb,ispin))
                occ0(iorb,ispin) = occ0(iorb,ispin)*(emax-emin)/(ne-1)
            enddo
        enddo
        occ0 = occ0/2

    end subroutine tbh_model_2

    subroutine dump_dos
        integer :: ispin,iorb,ie
        double precision :: de, nup, ndown

        if (master) then
            de = (emax-emin)/(ne-1)
            open(unit=IO_TB_DOS, file=FN_TB_DOS, status="replace")
            write(IO_TB_DOS,"(a1,I5,2F8.3,2I5)") "#", ne, emin, emax, norb, nspin
            do ie=1,ne
                write(IO_TB_DOS,"(F16.8,1x)",advance="no") egrid(ie)
                do ispin=1,nspin
                    do iorb=1,norb
                        write(IO_TB_DOS,"(F16.8,1x)",advance="no") &
                            dos(ie,iorb,ispin)
                    enddo
                enddo
                write(IO_TB_DOS,*)
            enddo
            close(IO_TB_DOS)

            ! print tight-binding occupancies
            write(*,*) "Tight-binding occupancies"
            write(*,*) "iorb    nup       ndown"

            do iorb=1,norb
                if (nspin==1) then
                    nup = occ0(iorb,1)
                    ndown = occ0(iorb,1)
                else
                    nup = occ0(iorb,1)
                    ndown = occ0(iorb,2)
                endif
                write(*,"(I5,4x,F6.4,4x,F6.4)") iorb,nup,ndown
            enddo
        endif

        call mpi_barrier(comm, mpierr)
    end subroutine dump_dos

end module dmft_lattice
