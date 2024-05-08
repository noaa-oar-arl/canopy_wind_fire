MODULE canopy_files_mod

!-------------------------------------------------------------------------------
! Name:     Files
! Purpose:  Contains FORTRAN units and file names.
! Revised:  15 Jul 2022  Original version.  (P. C. Campbell)
! Revised:  30 Nov 2023  Added supplementary canopy profile, file_canvars.  (P.C. Campbell)
!-------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER                       :: cdfid_m
    INTEGER,            PARAMETER :: max_mm       = 10000
    INTEGER,            PARAMETER :: iutnml       =  8

    CHARACTER(LEN=256)            :: file_vars    ( max_mm )
    CHARACTER(LEN=256)            :: file_canvars ( max_mm )
    CHARACTER(LEN=256)            :: file_out     ( 1 )
    CHARACTER(LEN=*), PARAMETER   :: file_nml     = 'input/namelist.canopy'


END MODULE canopy_files_mod
