MODULE canopy_canvars_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Variables Descriptions
! Purpose:  Contains canopy variables descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!-------------------------------------------------------------------------------
    use canopy_const_mod, ONLY: rk

    IMPLICIT NONE

!-------------------------------------------------------------------------------
! Canopy scalars in the model
!-------------------------------------------------------------------------------

    integer        ::    firetype      !1 = Above Canopy Fire; 0 = Below Canopy Fire
    integer        ::    cansublays    !Number of sub-canopy layers
    integer        ::    midflamepoint !Indice of the mid-flame point
    real(rk)       ::    cdrag         !Drag coefficient (nondimensional)
    real(rk)       ::    pai           !Plant/foliage area index (nondimensional)
    real(rk)       ::    zcanmax       !Height of maximum foliage area density (z/h) (nondimensional)
    real(rk)       ::    sigmau        !Standard deviation of shape function above zcanmax (z/h)
    real(rk)       ::    sigma1        !Standard deviation of shape function below zcanmax (z/h)

!-------------------------------------------------------------------------------
! Allocatable canopy variable arrays
!-------------------------------------------------------------------------------

    real(rk), allocatable :: zk                  ( : )          ! in-canopy heights (m)
    real(rk), allocatable :: zhc                 ( : )          ! z/h
    real(rk), allocatable :: fainc               ( : )          ! incremental foliage shape function
    real(rk), allocatable :: fafracz             ( : )          ! incremental fractional foliage shape function
    real(rk), allocatable :: fafraczInt          ( : )          ! integral of incremental fractional foliage shape function
    real(rk), allocatable :: fsun                ( : )          ! Sunlit/Shaded fraction from photolysis correction factor
    real(rk), allocatable :: tleaf_sun           ( : )          ! Leaf temp for sunlit leaves (K)
    real(rk), allocatable :: tleaf_shade         ( : )          ! Leaf temp for shaded leaves (K)
    real(rk), allocatable :: tleaf_ave           ( : )          ! Average Leaf temp for sunlit and shaded leaves (K)
    real(rk), allocatable :: ppfd_sun            ( : )          ! PPFD for sunlit leaves (umol phot/m2 s)
    real(rk), allocatable :: ppfd_shade          ( : )          ! PPFD for shaded leaves (umol phot/m2 s)
    real(rk), allocatable :: ppfd_ave            ( : )          ! Average PPFD for sunlit and shaded leaves (umol phot/m2 s)
    real(rk), allocatable :: tleaf_sun24_tmp     ( : , :, : )          ! Leaf temp for sunlit leaves (K)
    real(rk), allocatable :: tleaf_shade24_tmp   ( : , :, : )          ! Leaf temp for shaded leaves (K)
    real(rk), allocatable :: tleaf_ave24_tmp     ( : , :, : )          ! Average Leaf temp for sunlit and shaded leaves (K)
    real(rk), allocatable :: ppfd_sun24_tmp      ( : , :, : )          ! PPFD for sunlit leaves (umol phot/m2 s)
    real(rk), allocatable :: ppfd_shade24_tmp    ( : , :, : )          ! PPFD for shaded leaves (umol phot/m2 s)
    real(rk), allocatable :: tleaf_sun240_tmp    ( : , :, : )          ! Leaf temp for sunlit leaves (K)
    real(rk), allocatable :: tleaf_shade240_tmp  ( : , :, : )          ! Leaf temp for shaded leaves (K)
    real(rk), allocatable :: tleaf_ave240_tmp    ( : , :, : )          ! Average Leaf temp for sunlit and shaded leaves (K)
    real(rk), allocatable :: ppfd_sun240_tmp     ( : , :, : )          ! PPFD for sunlit leaves (umol phot/m2 s)
    real(rk), allocatable :: ppfd_shade240_tmp   ( : , :, : )          ! PPFD for shaded leaves (umol phot/m2 s)
    real(rk), allocatable :: tmp2mref_tmp        ( : , : )          ! 2-meter (AGL) input reference air temperature (K)
    real(rk), allocatable :: ubzref_tmp          ( : , : )          ! 10-meter (AGL) input reference wind speed (m/s)
    real(rk), allocatable :: tleaf_sun24         ( : , : )          ! Leaf temp for sunlit leaves (K)
    real(rk), allocatable :: tleaf_shade24       ( : , : )          ! Leaf temp for shaded leaves (K)
    real(rk), allocatable :: tleaf_ave24         ( : , : )          ! Average Leaf temp for sunlit and shaded leaves (K)
    real(rk), allocatable :: ppfd_sun24          ( : , : )          ! PPFD for sunlit leaves (umol phot/m2 s)
    real(rk), allocatable :: ppfd_shade24        ( : , : )          ! PPFD for shaded leaves (umol phot/m2 s)
    real(rk), allocatable :: tleaf_sun240        ( : , : )          ! Leaf temp for sunlit leaves (K)
    real(rk), allocatable :: tleaf_shade240      ( : , : )          ! Leaf temp for shaded leaves (K)
    real(rk), allocatable :: tleaf_ave240        ( : , : )          ! Average Leaf temp for sunlit and shaded leaves (K)
    real(rk), allocatable :: ppfd_sun240         ( : , : )          ! PPFD for sunlit leaves (umol phot/m2 s)
    real(rk), allocatable :: ppfd_shade240       ( : , : )          ! PPFD for shaded leaves (umol phot/m2 s)
    real(rk), allocatable :: daily_maxt2m        ( : )              ! Daily maximum 2-meter (AGL) input reference air temperature (K)
    real(rk), allocatable :: daily_mint2m        ( : )              ! Daily minimum 2-meter (AGL) input reference air temperature (K)
    real(rk), allocatable :: daily_maxws10m        ( : )            ! Daily maximum 10-meter (AGL) input wind speed (m/s)
    real(rk), allocatable :: tleaf_sun24_tmp_3d     ( : , : , : , : )          ! Leaf temp for sunlit leaves (K)
    real(rk), allocatable :: tleaf_shade24_tmp_3d   ( : , : , : , : )          ! Leaf temp for shaded leaves (K)
    real(rk), allocatable :: tleaf_ave24_tmp_3d     ( : , : , : , : )          ! Average Leaf temp for sunlit and shaded leaves (K)
    real(rk), allocatable :: ppfd_sun24_tmp_3d      ( : , : , : , : )          ! PPFD for sunlit leaves (umol phot/m2 s)
    real(rk), allocatable :: ppfd_shade24_tmp_3d    ( : , : , : , : )          ! PPFD for shaded leaves (umol phot/m2 s)
    real(rk), allocatable :: tleaf_sun240_tmp_3d    ( : , : , : , : )          ! Leaf temp for sunlit leaves (K)
    real(rk), allocatable :: tleaf_shade240_tmp_3d  ( : , : , : , : )          ! Leaf temp for shaded leaves (K)
    real(rk), allocatable :: tleaf_ave240_tmp_3d    ( : , : , : , : )          ! Average Leaf temp for sunlit and shaded leaves (K)
    real(rk), allocatable :: ppfd_sun240_tmp_3d     ( : , : , : , : )          ! PPFD for sunlit leaves (umol phot/m2 s)
    real(rk), allocatable :: ppfd_shade240_tmp_3d   ( : , : , : , : )          ! PPFD for shaded leaves (umol phot/m2 s)
    real(rk), allocatable :: tmp2mref_tmp_3d        ( : , : , : )              ! 2-meter (AGL) input reference air temperature (K)
    real(rk), allocatable :: ubzref_tmp_3d          ( : , : , : )              ! 10-meter (AGL) input reference wind speed (m/s)
    real(rk), allocatable :: tleaf_sun24_3d         ( : , : , : )          ! Leaf temp for sunlit leaves (K)
    real(rk), allocatable :: tleaf_shade24_3d       ( : , : , : )          ! Leaf temp for shaded leaves (K)
    real(rk), allocatable :: tleaf_ave24_3d         ( : , : , : )          ! Average Leaf temp for sunlit and shaded leaves (K)
    real(rk), allocatable :: ppfd_sun24_3d          ( : , : , : )          ! PPFD for sunlit leaves (umol phot/m2 s)
    real(rk), allocatable :: ppfd_shade24_3d        ( : , : , : )          ! PPFD for shaded leaves (umol phot/m2 s)
    real(rk), allocatable :: tleaf_sun240_3d        ( : , : , : )          ! Leaf temp for sunlit leaves (K)
    real(rk), allocatable :: tleaf_shade240_3d      ( : , : , : )          ! Leaf temp for shaded leaves (K)
    real(rk), allocatable :: tleaf_ave240_3d        ( : , : , : )          ! Average Leaf temp for sunlit and shaded leaves (K)
    real(rk), allocatable :: ppfd_sun240_3d         ( : , : , : )          ! PPFD for sunlit leaves (umol phot/m2 s)
    real(rk), allocatable :: ppfd_shade240_3d       ( : , : , : )          ! PPFD for shaded leaves (umol phot/m2 s)
    real(rk), allocatable :: daily_maxt2m_2d        ( : , : )                  ! Daily maximum 2-meter (AGL) input reference air temperature (K)
    real(rk), allocatable :: daily_mint2m_2d        ( : , : )                  ! Daily minimum 2-meter (AGL) input reference air temperature (K)
    real(rk), allocatable :: daily_maxws10m_2d      ( : , : )                  ! Daily maximum 10-meter (AGL) input wind speed (m/s)

    real(rk), allocatable :: canBOT                 ( : )                  ! Canopy bottom wind reduction factors
    real(rk), allocatable :: canTOP                 ( : )                  ! Canopy top wind reduction factors
    real(rk), allocatable :: canWIND                ( : , : )              ! canopy wind speeds (m/s)
    real(rk), allocatable :: canWIND_3d             ( : , : , : )          ! canopy wind speeds -- 3D (m/s)
    real(rk), allocatable :: lad                    ( : , : )              ! Leaf Area Density calculated from foliage shape function (m2/m3)
    real(rk), allocatable :: lad_3d                 ( : , : , : )          ! Leaf Area Density calculated from foliage shape function (m2/m3)
    real(rk), allocatable :: tka                    ( : , : )              ! Ambient temperature within/above canopy (K)
    real(rk), allocatable :: tka_3d                 ( : , : , : )          ! Ambient temperature within/above canopy (K)
    real(rk), allocatable :: pressa                 ( : , : )              ! Ambient pressure within/above canopy (mb)
    real(rk), allocatable :: pressa_3d              ( : , : , : )          ! Ambient pressure within/above canopy (mb)
    real(rk), allocatable :: relhuma                ( : , : )              ! Ambient relative humidity within/above canopy (%)
    real(rk), allocatable :: relhuma_3d             ( : , : , : )          ! Ambient relative humidity within/above canopy (%)
    real(rk), allocatable :: spechuma               ( : , : )              ! Ambient specific humidity within/above canopy (g/kg)
    real(rk), allocatable :: spechuma_3d            ( : , : , : )          ! Ambient specific humidity within/above canopy (g/kg)

    real(rk), allocatable :: dx                  ( : )          ! Model grid cell distance/resolution (m)
    real(rk), allocatable :: dx_2d               ( : , : )      ! Model grid cell distance/resolution -- 2D (m)
    real(rk), allocatable :: waf                 ( : )          ! Calculated Wind Adjustment Factor
    real(rk), allocatable :: waf_2d              ( : , : )      ! Calculated Wind Adjustment Factor -- 2D
    real(rk), allocatable :: d_h                 ( : )          ! Zero plane displacement heights (z/h)
    real(rk), allocatable :: d_h_2d              ( : , : )      ! Zero plane displacement heights (z/h) -- 2D
    real(rk), allocatable :: zo_h                ( : )          ! Surface (soil+veg) roughness lengths (z/h)
    real(rk), allocatable :: zo_h_2d             ( : , : )      ! Surface (soil+veg) roughness lengths (z/h) -- 2D
    real(rk), allocatable :: Kz                  ( :, : )       ! Eddy Diffusivities (m2/s)
    real(rk), allocatable :: Kz_3d               ( : , : , : )  ! Eddy Diffusivities -- 3D (m2/s)
    real(rk), allocatable :: rjcf                ( :, : )       ! Photolysis Attenuation Correction Factors
    real(rk), allocatable :: rjcf_3d             ( : , : , : )  ! Photolysis Attenuation Correction Factors -- 3D
    real(rk), allocatable :: flameh              ( : )          ! Flame Height (m)
    real(rk), allocatable :: flameh_2d           ( : , : )      ! Flame Height -- 2D (m)
    real(rk), allocatable :: emi_isop            ( :, : )       ! Isoprene biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_isop_3d         ( : , : , : )  ! Isoprene biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_myrc            ( :, : )       ! Myrcene biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_myrc_3d         ( : , : , : )  ! Myrcene biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_sabi            ( :, : )       ! Sabinene biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_sabi_3d         ( : , : , : )  ! Sabinene biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_limo            ( :, : )       ! Limonene biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_limo_3d         ( : , : , : )  ! Limonene biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_care            ( :, : )       ! Carene biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_care_3d         ( : , : , : )  ! Carene biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_ocim            ( :, : )       ! Ocimene biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_ocim_3d         ( : , : , : )  ! Ocimene biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_bpin            ( :, : )       ! Beta-Pinene biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_bpin_3d         ( : , : , : )  ! Beta-Pinene biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_apin            ( :, : )       ! Alpha-Pinene biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_apin_3d         ( : , : , : )  ! Alpha-Pinene biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_mono            ( :, : )       ! Other Monoterpenes biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_mono_3d         ( : , : , : )  ! Other Mononterpenes biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_farn            ( :, : )       ! Farnesene biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_farn_3d         ( : , : , : )  ! Farnesene biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_cary            ( :, : )       ! Caryophyllene biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_cary_3d         ( : , : , : )  ! Caryophyllene biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_sesq            ( :, : )       ! Other Sesquiterpenes biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_sesq_3d         ( : , : , : )  ! Other Sesquiterpenes biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_mbol            ( :, : )       ! 232-MBO biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_mbol_3d         ( : , : , : )  ! 232-MBO biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_meth            ( :, : )       ! Methanol biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_meth_3d         ( : , : , : )  ! Methanol biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_acet            ( :, : )       ! Acetone biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_acet_3d         ( : , : , : )  ! Acetone biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_co              ( :, : )       ! CO biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_co_3d           ( : , : , : )  ! CO biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_bvoc            ( :, : )       ! BIDI VOC biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_bvoc_3d         ( : , : , : )  ! BIDI VOC biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_svoc            ( :, : )       ! Stress VOC biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_svoc_3d         ( : , : , : )  ! Stress VOC biogenic emissions (kg/m2 s) -- 3D
    real(rk), allocatable :: emi_ovoc            ( :, : )       ! Other VOC biogenic emissions (kg/m2 s)
    real(rk), allocatable :: emi_ovoc_3d         ( : , : , : )  ! Other VOC biogenic emissions (kg/m2 s) -- 3D

    real(rk), allocatable :: ddep_o3             ( :, : )       ! Dry Deposition velocity for O3 (m/s)
    real(rk), allocatable :: ddep_o3_3d          ( : , : , : )  ! Dry Deposition velocity for O3 (m/s) -- 3D
!-------------------------------------------------------------------------------
! Canopy-App Program and version descriptors.
!-------------------------------------------------------------------------------

    CHARACTER(LEN=16),  PARAMETER     :: progname   = 'Canopy-App'
    CHARACTER(LEN=10),  PARAMETER     :: vdate      = '02/13/2022'
    CHARACTER(LEN=8),   PARAMETER     :: ver        = 'V1.0.0'

!-------------------------------------------------------------------------------
! Define output NETCDF data structures.
!-------------------------------------------------------------------------------

    TYPE fld1dtdata
        REAL(rk),        POINTER   :: fld        ( : )
        CHARACTER(LEN=16)          :: fldname
        CHARACTER(LEN=80)          :: long_name
        CHARACTER(LEN=80)          :: units
        CHARACTER(LEN=80)          :: cartesian_axis = "T"
        CHARACTER(LEN=80)          :: calendar_type = "JULIAN"
        CHARACTER(LEN=80)          :: calendar = "JULIAN"
        REAL                       :: fillvalue
        CHARACTER(LEN=16)          :: dimnames   ( 4 )
        INTEGER                    :: istart     ( 4 )
        INTEGER                    :: iend       ( 4 )
    END TYPE fld1dtdata

    TYPE fld1ddata
        REAL(rk),        POINTER   :: fld        ( : )
        CHARACTER(LEN=16)          :: fldname
        CHARACTER(LEN=80)          :: long_name
        CHARACTER(LEN=80)          :: units
        REAL                       :: fillvalue
        CHARACTER(LEN=16)          :: dimnames   ( 4 )
        INTEGER                    :: istart     ( 4 )
        INTEGER                    :: iend       ( 4 )
    END TYPE fld1ddata

    TYPE fld2ddata
        REAL(rk),        POINTER   :: fld        ( : , : )
        CHARACTER(LEN=16)          :: fldname
        CHARACTER(LEN=80)          :: long_name
        CHARACTER(LEN=16)          :: units
        REAL                       :: fillvalue
        CHARACTER(LEN=16)          :: dimnames   ( 4 )
        INTEGER                    :: istart     ( 4 )
        INTEGER                    :: iend       ( 4 )
    END TYPE fld2ddata

    TYPE fld3ddata
        REAL(rk),        POINTER   :: fld        ( : , : , : )
        CHARACTER(LEN=16)          :: fldname
        CHARACTER(LEN=80)          :: long_name
        CHARACTER(LEN=16)          :: units
        REAL                       :: fillvalue
        CHARACTER(LEN=16)          :: dimnames   ( 4 )
        INTEGER                    :: istart     ( 4 )
        INTEGER                    :: iend       ( 4 )
    END TYPE fld3ddata

!-------------------------------------------------------------------------------
! Assign number of time independent and varying 2D/3D fields at cell centers.
!-------------------------------------------------------------------------------

    INTEGER           :: nfld1dt        ! time field
    INTEGER           :: nfld1dz        ! time-independent 1d cell centers
    INTEGER           :: nfld2dxy       ! time-independent 2d cell centers
    INTEGER           :: nfld2dxyt      ! time-varying 2d cell centers
    INTEGER           :: nfld3dxyzt     ! time-varying 3d cell centers


!-------------------------------------------------------------------------------
! Time field.
!-------------------------------------------------------------------------------

    TYPE(fld1dtdata), ALLOCATABLE, TARGET :: fld1dt ( : )
    TYPE(fld1dtdata), POINTER     :: g_time

!-------------------------------------------------------------------------------
! Time-independent 1d fields at cell centers.
!-------------------------------------------------------------------------------

    TYPE(fld1ddata), ALLOCATABLE, TARGET :: fld1dz ( : )
    TYPE(fld1ddata), POINTER     :: g_level

!-------------------------------------------------------------------------------
! Time-independent 2d fields at cell centers.
!-------------------------------------------------------------------------------

    TYPE(fld2ddata), ALLOCATABLE, TARGET :: fld2dxy ( : )
    TYPE(fld2ddata), POINTER     :: g_lat
    TYPE(fld2ddata), POINTER     :: g_lon

!-------------------------------------------------------------------------------
! Time-varying 2d fields at cell centers for output NETCDF
!-------------------------------------------------------------------------------

    TYPE(fld2ddata), ALLOCATABLE, TARGET :: fld2dxyt ( : )
    TYPE(fld2ddata), POINTER     :: c_waf
    TYPE(fld2ddata), POINTER     :: c_flameh
    TYPE(fld2ddata), POINTER     :: c_canheight
    TYPE(fld2ddata), POINTER     :: c_dh
    TYPE(fld2ddata), POINTER     :: c_zoh

!-------------------------------------------------------------------------------
! Time-varying 3d fields at cell centers for output NETCDF
!-------------------------------------------------------------------------------

    TYPE(fld3ddata), ALLOCATABLE, TARGET :: fld3dxyzt ( : )
    TYPE(fld3ddata), POINTER     :: c_lad
    TYPE(fld3ddata), POINTER     :: c_canwind
    TYPE(fld3ddata), POINTER     :: c_Kz
    TYPE(fld3ddata), POINTER     :: c_rjcf
    TYPE(fld3ddata), POINTER     :: c_emi_isop
    TYPE(fld3ddata), POINTER     :: c_emi_myrc
    TYPE(fld3ddata), POINTER     :: c_emi_sabi
    TYPE(fld3ddata), POINTER     :: c_emi_limo
    TYPE(fld3ddata), POINTER     :: c_emi_care
    TYPE(fld3ddata), POINTER     :: c_emi_ocim
    TYPE(fld3ddata), POINTER     :: c_emi_bpin
    TYPE(fld3ddata), POINTER     :: c_emi_apin
    TYPE(fld3ddata), POINTER     :: c_emi_mono
    TYPE(fld3ddata), POINTER     :: c_emi_farn
    TYPE(fld3ddata), POINTER     :: c_emi_cary
    TYPE(fld3ddata), POINTER     :: c_emi_sesq
    TYPE(fld3ddata), POINTER     :: c_emi_mbol
    TYPE(fld3ddata), POINTER     :: c_emi_meth
    TYPE(fld3ddata), POINTER     :: c_emi_acet
    TYPE(fld3ddata), POINTER     :: c_emi_co
    TYPE(fld3ddata), POINTER     :: c_emi_bvoc
    TYPE(fld3ddata), POINTER     :: c_emi_svoc
    TYPE(fld3ddata), POINTER     :: c_emi_ovoc

END MODULE canopy_canvars_mod
