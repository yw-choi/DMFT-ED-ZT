module dmft_green
    use mpi
    use io_units
    use dmft_params, only: norb, nbath, nspin, read_gf_from_file
    use dmft_grid, only: nwloc, omega

    implicit none

    public :: &
        initial_green_ftn, &
        calc_delta

    private 
contains

    subroutine initial_green_ftn(G, G0, Sigma)
        double complex, intent(inout) :: &
            G(nwloc,norb,nspin), &
            G0(nwloc,norb,nspin), &
            Sigma(nwloc,norb,nspin)

        logical :: found

        found = .false.
        if (read_gf_from_file) then
            ! call import_green_ftn(found, G, G0, Sigma)
        endif

        if (found) then
            return
        endif

        ! call local_green_ftn(G, Sigma)

    end subroutine initial_green_ftn

    ! calculates hybridization function from ek, vmk 
    subroutine calc_delta(ek, vmk, delta)
        double precision, intent(in) :: &
            ek(nbath,2), vmk(norb,nbath,2)
        double complex, intent(out) :: &
            Delta(nwloc,norb,nspin)

        integer :: iw,iorb,ispin,ibath
        double complex :: d, iwn

        do ispin=1,nspin
            do iorb=1,norb
                do iw=1,nwloc
                    Delta(iw,iorb,ispin) = 0.d0
                    do ibath=1,nbath
                        iwn = cmplx(0.d0,omega(iw)) 
                        d = vmk(iorb,ibath,ispin)**2/(iwn-ek(ibath,ispin))
                        Delta(iw,iorb,ispin) = Delta(iw,iorb,ispin) + d
                    enddo
                enddo
            enddo
        enddo
    end subroutine calc_Delta 

end module dmft_green
