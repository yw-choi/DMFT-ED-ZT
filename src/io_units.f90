module io_units

    implicit none

    integer, parameter :: IO_TB_DOS = 11
    character(len=100), parameter :: FN_TB_DOS = "tbdos.dat"

    integer, parameter :: IO_GR_DATA     = 100
    character(len=100), parameter :: FN_GR_SAVE = "green.save"

    integer, parameter :: IO_H_PARAMS    = 101
    integer, parameter :: IO_MEM_REPORT  = 102
    integer, parameter :: IO_G_COEFFS    = 103
end module io_units
