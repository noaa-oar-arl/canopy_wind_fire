MODULE canopy_canmet_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Met/Sfc Input Variable Descriptions
! Purpose:  Contains canopy met and sfc input  variable descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!-------------------------------------------------------------------------------
    use canopy_const_mod, ONLY: rk

    IMPLICIT NONE
!! .... defines canopy options (read from user namelist)

    ! Generic 2D met/sfc input variables that should be passed to canopy calculations
    TYPE :: variable_type
        real(rk)   :: lat          !latitude of cell/point
        real(rk)   :: lon          !longitude of cell/point
        real(rk)   :: ch           !canopy height (m)
        real(rk)   :: ugrd10m      !u wind speed at reference height above canopy (m/s)
        real(rk)   :: vgrd10m      !v wind speed at reference height above canopy (m/s)
        real(rk)   :: clu          !clumping index
        real(rk)   :: lai          !leaf area index
        integer    :: vtype        !vegetation type
        real(rk)   :: canfrac      !canopy fraction
        real(rk)   :: fricv        !friction velocity (u*) (m/s)
        real(rk)   :: csz          !cosine of solar zenith angle
        real(rk)   :: sfcr         !surface roughness length (m)
        real(rk)   :: mol          !Monin-Obukhov length (m)
        real(rk)   :: frp          !fire radiative power (MW)
        real(rk)   :: href         !reference height above the canopy (m)
        integer    :: sotyp        !soil type
        real(rk)   :: pressfc      !surface pressure (hPa)
        real(rk)   :: dswrf        !instantaneous downward shortwave radiation (W/m2)
        real(rk)   :: shtfl        !instantaneous surface sensible heat net flux (W/m2)
        real(rk)   :: tmpsfc       !surface temperature (K)
        real(rk)   :: tmp2m        !2-meter temperature (K)
        real(rk)   :: spfh2m       !2-meter specific humidity (kg/kg)
        real(rk)   :: hpbl         !height of planetary boundary layer (m)
        real(rk)   :: prate_ave    !mass precipitation rate (kg/m2 s)
        real(rk)   :: soilw1       !volumetric soil moisture level 1 (m3/m3)
        real(rk)   :: soilw2       !volumetric soil moisture level 2 (m3/m3)
        real(rk)   :: soilw3       !volumetric soil moisture level 3 (m3/m3)
        real(rk)   :: soilw4       !volumetric soil moisture level 4 (m3/m3)
        real(rk)   :: wilt         !wilting point (proportion)
        real(rk)   :: ozone_w126   !ozone W126 values (ppm-hours)

    end TYPE variable_type

    type(variable_type), allocatable :: variables( : ), variables_2d( : , :)

    ! Generic 3D input variables
    TYPE :: variable_type_1d
        real(rk)   :: lev         !Input mid-level heights associated with 3D input option below
    end TYPE variable_type_1d

    type(variable_type_1d), allocatable :: variables_1d( : )

    TYPE :: variable_type_3d
        real(rk)   :: pavd         !Plant Area Volume Density (PAVD) profile (m2/m3)
    end TYPE variable_type_3d

    type(variable_type_3d), allocatable :: variables_3d( : , : , :)

    ! Generic set of observed canopy profile input variable levels (14) from point text file
    TYPE :: variable_type_can
        real(rk)   :: lat    !latitude of cell/point
        real(rk)   :: lon    !longitude of cell/point
        real(rk)   :: lev01  !Input canopy profile levels
        real(rk)   :: pavd01 !!Input canopy PAVD profile
        real(rk)   :: lev02
        real(rk)   :: pavd02
        real(rk)   :: lev03
        real(rk)   :: pavd03
        real(rk)   :: lev04
        real(rk)   :: pavd04
        real(rk)   :: lev05
        real(rk)   :: pavd05
        real(rk)   :: lev06
        real(rk)   :: pavd06
        real(rk)   :: lev07
        real(rk)   :: pavd07
        real(rk)   :: lev08
        real(rk)   :: pavd08
        real(rk)   :: lev09
        real(rk)   :: pavd09
        real(rk)   :: lev10
        real(rk)   :: pavd10
        real(rk)   :: lev11
        real(rk)   :: pavd11
        real(rk)   :: lev12
        real(rk)   :: pavd12
        real(rk)   :: lev13
        real(rk)   :: pavd13
        real(rk)   :: lev14
        real(rk)   :: pavd14
    end TYPE variable_type_can

    type(variable_type_can), allocatable :: variables_can( : )


    ! Met/Sfc variable reassignment names  above reference conditions from the model
    real(rk)       ::    latref          !latitude of cell/point
    real(rk)       ::    lonref          !longitude of cell/point
    real(rk)       ::    hcmref          !Input Canopy Height
    real(rk)       ::    uref            !Input above canopy/reference 10-m U wind speed
    real(rk)       ::    vref            !Input above canopy/reference 10-m V wind speed
    real(rk)       ::    ubzref          !Input above canopy/reference 10-m model wind speed
    real(rk)       ::    cluref          !Input canopy clumping index
    real(rk)       ::    lairef          !Input leaf area index
    integer        ::    vtyperef        !Input vegetation type (VIIRS)
    real(rk)       ::    canfracref      !Input canopy fraction of grid cell
    real(rk)       ::    ustref          !Input friction velocity
    real(rk)       ::    cszref          !Input cosine of zenith angle
    real(rk)       ::    z0ref           !Input total/surface roughness length
    real(rk)       ::    molref          !Input Monin-Obukhov Length
    real(rk)       ::    frpref          !Input fire radiative power
    real(rk)       ::    hgtref          !reference height above the canopy (m)
    integer        ::    sotypref        !soil type
    real(rk)       ::    pressfcref      !surface pressure (hPa)
    real(rk)       ::    dswrfref        !instantaneous downward shortwave radiation (W/m2)
    real(rk)       ::    shtflref        !instantaneous surface sensible heat net flux (W/m2)
    real(rk)       ::    tmpsfcref       !surface temperature (K)
    real(rk)       ::    tmp2mref        !2-meter temperature (K)
    real(rk)       ::    spfh2mref       !2-meter specific humidity (kg/kg)
    real(rk)       ::    hpblref         !height of planetary boundary layer (m)
    real(rk)       ::    prate_averef    !mass precipitation rate (kg/m2 s)
    real(rk)       ::    soilw1ref       !volumetric soil moisture (m3/m3) Layer 1
    real(rk)       ::    soilw2ref       !volumetric soil moisture (m3/m3) Layer 2
    real(rk)       ::    soilw3ref       !volumetric soil moisture (m3/m3) Layer 3
    real(rk)       ::    soilw4ref       !volumetric soil moisture (m3/m3) Layer 4
    real(rk)       ::    wiltref         !wilting point (proportion)
    real(rk)       ::    ozone_w126ref   !ozone W126 values (ppm-hours)

!    real(rk)       ::    lev01ref, lev02ref, lev03ref, lev04ref, lev05ref, & !Input canopy profile levels
!                         lev06ref, lev07ref, lev08ref, lev09ref, lev10ref, &
!                         lev11ref, lev12ref, lev13ref, lev14ref
!    real(rk)       ::    pavd01ref, pavd02ref, pavd03ref, pavd04ref, pavd05ref, & !Input canopy PAVD profile
!                         pavd06ref, pavd07ref, pavd08ref, pavd09ref, pavd10ref, &
!                         pavd11ref, pavd12ref, pavd13ref, pavd14ref
    real(rk), allocatable   :: pavdref ( : ), pavd_arr ( : ) !plant area volume density (m2/m3)
    real(rk), allocatable   :: levref ( : ), lev_arr  ( : ) !reference vertical levels with 3d input data

END MODULE canopy_canmet_mod
