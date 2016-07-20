module ed_basis

    use mpi
    use dmft_params, only: norb
    use ed_params, only: nbath, KIND_BASIS, nsite
    use numeric_utils, only: icom
    use alloc, only: re_alloc, de_alloc
    
    implicit none

    public :: &
        generate_basis, &
        ed_basis_get,   &
        ed_basis_idx,   &
        get_bitidx,     &
        dealloc_basis

    type, public :: basis_t
        integer(kind=KIND_BASIS) :: nloc

        integer(kind=KIND_BASIS) :: ntot
        integer(kind=KIND_BASIS) :: nup
        integer(kind=KIND_BASIS) :: ndown

        integer :: ne_up
        integer :: ne_down

        integer(kind=KIND_BASIS), pointer :: up(:)
        integer(kind=KIND_BASIS), pointer :: down(:)

        integer(kind=KIND_BASIS), pointer :: idx_up(:)
        integer(kind=KIND_BASIS), pointer :: idx_down(:)

        integer(kind=KIND_BASIS), allocatable :: nlocals(:)
        integer(kind=KIND_BASIS), allocatable :: offsets(:)
    end type basis_t

    private
contains

    subroutine generate_basis( ne_up, ne_down, basis )
        integer, intent(in) :: ne_up, ne_down
        type(basis_t), intent(out) :: basis

        ! local variables
        integer(kind=KIND_BASIS) :: ispin, i, j, nam
        integer(kind=KIND_BASIS) :: minrange, maxrange, counts, nbit
        integer(kind=KIND_BASIS) :: nud(2)

        basis%ne_up = ne_up
        basis%ne_down = ne_down

        basis%nup   = icom(nsite,ne_up)
        basis%ndown = icom(nsite,ne_down)
        basis%ntot  = basis%nup * basis%ndown

        basis%nloc  = basis%ntot/nprocs
        nam = mod(basis%ntot,nprocs)
        if (taskid.lt.nam) basis%nloc = basis%nloc + 1

        allocate(basis%nlocals(0:nprocs-1),basis%offsets(0:nprocs-1))
        call mpi_allgather(basis%nloc,1,mpi_integer,basis%nlocals(0),1,mpi_integer,comm,mpierr)

        basis%offsets(0) = 0 
        do i = 1, nprocs-1
            basis%offsets(i) = basis%offsets(i-1) + basis%nlocals(i-1)
        enddo

        nullify(basis%up,basis%down,basis%idx_down,basis%idx_up)

        call re_alloc(basis%up, 1, basis%nup, 'ed_basis', 'basis_up')
        call re_alloc(basis%down, 1, basis%ndown, 'ed_basis', 'basis_down')

        nud(1) = ne_up
        nud(2) = ne_down

        do ispin=1,2
            ! ref : arXiv:1307.7542, Appendix A
            minrange = 0
            maxrange = 0

            do i=1,nud(ispin)
                minrange = minrange + 2**(i-1)
                maxrange = maxrange + 2**(nsite-i)
            enddo
            
            if (ispin.eq.1) then
                call re_alloc(basis%idx_up, minrange, maxrange, &
                    'ed_basis', 'basis_idx_up')
            else
                call re_alloc(basis%idx_down, minrange, maxrange, &
                    'ed_basis', 'basis_idx_down')
            endif
            
            counts = 0

            do i=minrange,maxrange
                nbit = 0
                do j=0,nsite-1
                    if (BTEST(i,j)) then
                        nbit = nbit + 1
                    endif
                enddo

                if (nbit.eq.nud(ispin)) then
                    counts = counts + 1
                    if (ispin.eq.1) then
                        basis%up(counts) = i
                        basis%idx_up(i) = counts
                    else
                        basis%down(counts) = i
                        basis%idx_down(i) = counts
                    endif
                endif
            enddo
        enddo
    end subroutine generate_basis

    subroutine dealloc_basis(basis)
        type(basis_t), intent(inout) :: basis
        call de_alloc(basis%up, 'ed_basis', 'basis_up')
        call de_alloc(basis%down, 'ed_basis', 'basis_down')
        call de_alloc(basis%idx_up, 'ed_basis', 'basis_idx_up')
        call de_alloc(basis%idx_down, 'ed_basis', 'basis_idx_down')
    end subroutine dealloc_basis

    ! ref : arXiv:1307.7542 eq (6)
    integer(kind=kind_basis) function ed_basis_get(basis,idx_loc) 
        type(basis_t), intent(in) :: basis
        integer(kind=KIND_BASIS), intent(in) :: idx_loc

        ! local variables
        integer(kind=KIND_BASIS) :: iup, idown, idx

        ! local idx to global idx
        idx = idx_loc + basis%offsets(taskid)

        iup = mod(idx-1,basis%nup)+1
        idown = (idx-1)/basis%nup+1

        ed_basis_get = basis%up(iup)+2**(nsite)*basis%down(idown)
    end function ed_basis_get

    ! ref : arXiv:1307.7542 eq (8)
    integer(kind=KIND_BASIS) function ed_basis_idx(basis, basis_i)
        type(basis_t), intent(in) :: basis
        integer(kind=kind_basis) :: basis_i

        ! local variables
        integer(kind=KIND_BASIS) :: basis_i_up, basis_i_down
        
        basis_i_up = mod(basis_i,2**(nsite))
        basis_i_down = basis_i/(2**(nsite))

        ed_basis_idx = (basis%idx_down(basis_i_down)-1)*basis%nup + &
                             basis%idx_up(basis_i_up)
    end function ed_basis_idx

    integer(kind=KIND_BASIS) function get_bitidx(isite,ispin)
        integer(kind=KIND_BASIS) :: isite, ispin
        get_bitidx = (ispin-1)*Nsite + isite-1
    end function get_bitidx
end module ed_basis
