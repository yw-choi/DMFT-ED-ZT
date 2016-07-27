module io_units

    implicit none

    integer, parameter :: IO_TB_DOS = 11
    character(len=100), parameter :: FN_TB_DOS = "tbdos.dat"

    integer, parameter :: IO_GR     = 100
    character(len=100), parameter :: FN_GR = "green.save"

    integer, parameter :: IO_H_PARAMS    = 101
    character(len=100), parameter :: FN_H_PARAMS = "h_imp.save"

    integer, parameter :: IO_MEM_REPORT  = 102
    integer, parameter :: IO_G_COEFFS    = 103
end module io_units
