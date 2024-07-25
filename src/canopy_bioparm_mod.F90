module canopy_bioparm_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_BIOP( EMI_IND, LU_OPT, VTYPE, &
        EF, CT1, CEO, ANEW, AGRO, AMAT, AOLD, ROOTA, ROOTB, &
        CAQ, TAQ, DTAQ, CHT, THT, DTHT, CLT, TLT, DTLT, &
        CHW, THW, DTHW)

!-----------------------------------------------------------------------

! Description:
!     gets biogenic emissions factors and parameters from MEGAN2.1

! Preconditions:
!     emissions index, lu option, and vegtype

! Subroutines and Functions Called:

! Revision History:
!     Prototype 02/23 by PCC, based on Guenther et al. (2012)
! Citation:
!Guenther, A. B., et al.,: The Model of Emissions of Gases and Aerosols from
!Nature version 2.1 (MEGAN2.1): an extended and updated framework for
!modeling biogenic emissions, Geosci. Model Dev., 5, 1471–1492,
!https://doi.org/10.5194/gmd-5-1471-2012, 2012.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Feb 2023 P.C. Campbell: Initial version
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

        use canopy_const_mod, ONLY: rk

! Arguments:
!     IN/OUT
        INTEGER,     INTENT( IN )       :: EMI_IND         ! Input biogenic emissions index
        INTEGER,     INTENT( IN )       :: LU_OPT          ! integer for LU type from model mapped to Massman et al. (default = 0/VIIRS)
        INTEGER,     INTENT( IN )       :: VTYPE           ! Grid cell dominant vegetation type
        REAL(RK),    INTENT( OUT )      :: EF              ! Out Mapped EF ((ug/m2 hr)
        REAL(RK),    INTENT( OUT )      :: CT1             ! Out Activation energy (kJ/mol)
        REAL(RK),    INTENT( OUT )      :: CEO             ! Out Empirical coefficient
        REAL(RK),    INTENT( OUT )      :: ANEW, AGRO, AMAT, AOLD   !Empirical factors or coefficients for: growing, mature, and old/senescing foliage, as per Table 4 of Guenther et al., 2012
        REAL(RK),    INTENT( OUT )      :: ROOTA, ROOTB    ! Coefficients A and B used for PFT dependent cumulative root depth fraction [m-1]
        REAL(RK),    INTENT( OUT )      :: CAQ             ! coefficient for poor Air Quality stress
        REAL(RK),    INTENT( OUT )      :: TAQ             ! threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),    INTENT( OUT )      :: DTAQ            ! delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),    INTENT( OUT )      :: CHT             ! coefficient for high temperature stress
        REAL(RK),    INTENT( OUT )      :: THT             ! threshold for high temperature stress (K)
        REAL(RK),    INTENT( OUT )      :: DTHT            ! delta threshold high temperature stress (K)
        REAL(RK),    INTENT( OUT )      :: CLT             ! coefficient for low temperature stress
        REAL(RK),    INTENT( OUT )      :: TLT             ! threshold for low temperature stress (K)
        REAL(RK),    INTENT( OUT )      :: DTLT            ! delta threshold low temperature stress (K)
        REAL(RK),    INTENT( OUT )      :: CHW             ! coefficient for high wind stress
        REAL(RK),    INTENT( OUT )      :: THW             ! threshold for high wind stress (m/s)
        REAL(RK),    INTENT( OUT )      :: DTHW            ! delta threshold high wind stress (m/s)

!  LOCAL
        REAL(RK) :: EF1,EF2,EF3,EF4,EF5,EF6,EF7    ! Plant Emission factors (EF) (ug/m2 hr)
        REAL(RK) :: EF8,EF9,EF10,EF11,EF12,EF13    !
        REAL(RK) :: EF14,EF15                      !


! Plant-Dependent emissions capacity/factors (EFs) for Isoprene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_ISOP    =  600.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_ISOP    =  3000.0_rk    ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_ISOP    =  1.0_rk       ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_ISOP    =  7000.0_rk    ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_ISOP    =  10000.0_rk   ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_ISOP    =  7000.0_rk    ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_ISOP    =  10000.0_rk   ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_ISOP    =  11000.0_rk   ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_ISOP    =  2000.0_rk    ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_ISOP   =  4000.0_rk    ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_ISOP   =  4000.0_rk    ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_ISOP   =  1600.0_rk    ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_ISOP   =  800.0_rk     ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_ISOP   =  200.0_rk     ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_ISOP   =  1.0_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Isoprene as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_ISOP  = 0.05_rk
        REAL(RK),          PARAMETER     :: AGRO_ISOP  = 0.6_rk
        REAL(RK),          PARAMETER     :: AMAT_ISOP  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_ISOP  = 0.9_rk

! Plant-Dependent emissions capacity/factors (EFs) for Myrcene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_MYRC    =  70.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_MYRC    =  70.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_MYRC    =  60.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_MYRC    =  80.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_MYRC    =  30.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_MYRC    =  80.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_MYRC    =  30.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_MYRC    =  30.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_MYRC    =  30.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_MYRC   =  50.0_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_MYRC   =  30.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_MYRC   =  0.3_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_MYRC   =  0.3_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_MYRC   =  0.3_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_MYRC   =  0.3_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Myrcene as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_MYRC  = 2.0_rk
        REAL(RK),          PARAMETER     :: AGRO_MYRC  = 1.8_rk
        REAL(RK),          PARAMETER     :: AMAT_MYRC  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_MYRC  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for Sabinene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_SABI    =  70.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_SABI    =  70.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_SABI    =  40.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_SABI    =  80.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_SABI    =  50.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_SABI    =  80.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_SABI    =  50.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_SABI    =  50.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_SABI    =  50.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_SABI   =  70.0_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_SABI   =  50.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_SABI   =  0.7_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_SABI   =  0.7_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_SABI   =  0.7_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_SABI   =  0.7_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Sabinene as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_SABI  = 2.0_rk
        REAL(RK),          PARAMETER     :: AGRO_SABI  = 1.8_rk
        REAL(RK),          PARAMETER     :: AMAT_SABI  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_SABI  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for Limonene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_LIMO    =  100.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_LIMO    =  100.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_LIMO    =  130.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_LIMO    =  80.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_LIMO    =  80.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_LIMO    =  80.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_LIMO    =  80.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_LIMO    =  80.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_LIMO    =  60.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_LIMO   =  100.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_LIMO   =  60.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_LIMO   =  0.7_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_LIMO   =  0.7_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_LIMO   =  0.7_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_LIMO   =  0.7_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Limonene as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_LIMO  = 2.0_rk
        REAL(RK),          PARAMETER     :: AGRO_LIMO  = 1.8_rk
        REAL(RK),          PARAMETER     :: AMAT_LIMO  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_LIMO  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for 3-Carene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_CARE    =  160.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_CARE    =  160.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_CARE    =  80.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_CARE    =  40.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_CARE    =  30.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_CARE    =  40.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_CARE    =  30.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_CARE    =  30.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_CARE    =  30.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_CARE   =  100.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_CARE   =  30.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_CARE   =  0.3_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_CARE   =  0.3_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_CARE   =  0.3_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_CARE   =  0.3_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for 3-Carene as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_CARE  = 2.0_rk
        REAL(RK),          PARAMETER     :: AGRO_CARE  = 1.8_rk
        REAL(RK),          PARAMETER     :: AMAT_CARE  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_CARE  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for t-Beta-Ocimene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_OCIM    =  70.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_OCIM    =  70.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_OCIM    =  60.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_OCIM    =  150.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_OCIM    =  120.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_OCIM    =  150.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_OCIM    =  120.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_OCIM    =  120.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_OCIM    =  90.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_OCIM   =  150.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_OCIM   =  90.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_OCIM   =  2.0_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_OCIM   =  2.0_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_OCIM   =  2.0_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_OCIM   =  2.0_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for t-Beta-Ocimene as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_OCIM  = 2.0_rk
        REAL(RK),          PARAMETER     :: AGRO_OCIM  = 1.8_rk
        REAL(RK),          PARAMETER     :: AMAT_OCIM  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_OCIM  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for Beta-Pinene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_BPIN    =  300.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_BPIN    =  300.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_BPIN    =  200.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_BPIN    =  120.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_BPIN    =  130.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_BPIN    =  120.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_BPIN    =  130.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_BPIN    =  130.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_BPIN    =  100.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_BPIN   =  150.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_BPIN   =  100.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_BPIN   =  1.5_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_BPIN   =  1.5_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_BPIN   =  1.5_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_BPIN   =  1.5_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Beta-Pinene as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_BPIN  = 2.0_rk
        REAL(RK),          PARAMETER     :: AGRO_BPIN  = 1.8_rk
        REAL(RK),          PARAMETER     :: AMAT_BPIN  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_BPIN  = 1.05_rk


! Plant-Dependent emissions capacity/factors (EFs) for Alpha-Pinene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_APIN    =  500.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_APIN    =  500.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_APIN    =  510.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_APIN    =  600.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_APIN    =  400.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_APIN    =  600.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_APIN    =  400.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_APIN    =  400.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_APIN    =  200.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_APIN   =  300.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_APIN   =  200.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_APIN   =  2.0_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_APIN   =  2.0_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_APIN   =  2.0_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_APIN   =  2.0_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Alpha-Pinene as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_APIN  = 2.0_rk
        REAL(RK),          PARAMETER     :: AGRO_APIN  = 1.8_rk
        REAL(RK),          PARAMETER     :: AMAT_APIN  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_APIN  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for Other Monoterpenes (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! ! Other Monoterpenes category (34 compounds):  See Table 1 of Guenther et al. (2012)
        REAL(RK),          PARAMETER     :: EF1_MONO    =  180.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_MONO    =  180.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_MONO    =  170.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_MONO    =  150.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_MONO    =  150.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_MONO    =  150.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_MONO    =  150.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_MONO    =  150.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_MONO    =  110.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_MONO   =  200.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_MONO   =  110.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_MONO   =  5.0_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_MONO   =  5.0_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_MONO   =  5.0_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_MONO   =  5.0_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Other Monoterpenes as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_MONO  = 2.0_rk
        REAL(RK),          PARAMETER     :: AGRO_MONO  = 1.8_rk
        REAL(RK),          PARAMETER     :: AMAT_MONO  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_MONO  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for Alpha-Farnesene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_FARN    =  40.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_FARN    =  40.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_FARN    =  40.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_FARN    =  60.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_FARN    =  40.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_FARN    =  60.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_FARN    =  40.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_FARN    =  40.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_FARN    =  40.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_FARN   =  40.0_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_FARN   =  40.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_FARN   =  3.0_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_FARN   =  3.0_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_FARN   =  3.0_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_FARN   =  4.0_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Alpha-Farnesene as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_FARN  = 0.4_rk
        REAL(RK),          PARAMETER     :: AGRO_FARN  = 0.6_rk
        REAL(RK),          PARAMETER     :: AMAT_FARN  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_FARN  = 0.95_rk

! Plant-Dependent emissions capacity/factors (EFs) for Beta-Caryophyllene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_CARY    =  80.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_CARY    =  80.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_CARY    =  80.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_CARY    =  60.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_CARY    =  40.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_CARY    =  60.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_CARY    =  40.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_CARY    =  40.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_CARY    =  50.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_CARY   =  50.0_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_CARY   =  50.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_CARY   =  1.0_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_CARY   =  1.0_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_CARY   =  1.0_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_CARY   =  4.0_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Beta-Caryophyllene as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_CARY  = 0.4_rk
        REAL(RK),          PARAMETER     :: AGRO_CARY  = 0.6_rk
        REAL(RK),          PARAMETER     :: AMAT_CARY  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_CARY  = 0.95_rk

! Plant-Dependent emissions capacity/factors (EFs) for Other Sesquieterpenes (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! Other Sesquiterpenes category (30 compounds):  See Table 1 of Guenther et al. (2012)
        REAL(RK),          PARAMETER     :: EF1_SESQ    =  120.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_SESQ    =  120.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_SESQ    =  120.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_SESQ    =  120.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_SESQ    =  100.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_SESQ    =  120.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_SESQ    =  100.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_SESQ    =  100.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_SESQ    =  100.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_SESQ   =  100.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_SESQ   =  100.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_SESQ   =  1.0_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_SESQ   =  1.0_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_SESQ   =  1.0_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_SESQ   =  1.0_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Other Sesquieterpenes as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_SESQ  = 0.4_rk
        REAL(RK),          PARAMETER     :: AGRO_SESQ  = 0.6_rk
        REAL(RK),          PARAMETER     :: AMAT_SESQ  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_SESQ  = 0.95_rk

! Plant-Dependent emissions capacity/factors (EFs) for 232-MBO (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_MBOL    =  700.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_MBOL    =  60.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_MBOL    =  0.01_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_MBOL    =  0.01_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_MBOL    =  0.01_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_MBOL    =  0.01_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_MBOL    =  0.01_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_MBOL    =  2.0_rk       ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_MBOL    =  0.01_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_MBOL   =  0.01_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_MBOL   =  0.01_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_MBOL   =  0.01_rk      ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_MBOL   =  0.01_rk      ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_MBOL   =  0.01_rk      ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_MBOL   =  0.01_rk      ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for 232-MBO as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_MBOL  = 0.05_rk
        REAL(RK),          PARAMETER     :: AGRO_MBOL  = 0.6_rk
        REAL(RK),          PARAMETER     :: AMAT_MBOL  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_MBOL  = 0.9_rk

! Plant-Dependent emissions capacity/factors (EFs) for Methanol (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_METH    =  900.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_METH    =  900.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_METH    =  900.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_METH    =  500.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_METH    =  900.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_METH    =  500.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_METH    =  900.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_METH    =  900.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_METH    =  900.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_METH   =  900.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_METH   =  900.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_METH   =  500.0_rk     ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_METH   =  500.0_rk     ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_METH   =  500.0_rk     ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_METH   =  900.0_rk     ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Methanol as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_METH  = 3.5_rk
        REAL(RK),          PARAMETER     :: AGRO_METH  = 3.0_rk
        REAL(RK),          PARAMETER     :: AMAT_METH  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_METH  = 1.2_rk

! Plant-Dependent emissions capacity/factors (EFs) for Acetone (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_ACET    =  240.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_ACET    =  240.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_ACET    =  240.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_ACET    =  240.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_ACET    =  240.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_ACET    =  240.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_ACET    =  240.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_ACET    =  240.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_ACET    =  240.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_ACET   =  240.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_ACET   =  240.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_ACET   =  80.0_rk      ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_ACET   =  80.0_rk      ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_ACET   =  80.0_rk      ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_ACET   =  80.0_rk      ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Acetone as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_ACET  = 1.0_rk
        REAL(RK),          PARAMETER     :: AGRO_ACET  = 1.0_rk
        REAL(RK),          PARAMETER     :: AMAT_ACET  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_ACET  = 1.0_rk

! Plant-Dependent emissions capacity/factors (EFs) for Carbon Monoxide (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_CO      =  600.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_CO      =  600.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_CO      =  600.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_CO      =  600.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_CO      =  600.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_CO      =  600.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_CO      =  600.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_CO      =  600.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_CO      =  600.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_CO     =  600.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_CO     =  600.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_CO     =  600.0_rk     ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_CO     =  600.0_rk     ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_CO     =  600.0_rk     ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_CO     =  600.0_rk     ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Carbon Monoxide as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_CO  = 1.0_rk
        REAL(RK),          PARAMETER     :: AGRO_CO  = 1.0_rk
        REAL(RK),          PARAMETER     :: AMAT_CO  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_CO  = 1.0_rk

! Plant-Dependent emissions capacity/factors (EFs) for BIDI VOC species (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! Bidirectional VOC (5 compounds): See Table 1 of Guenther et al. (2012)
        REAL(RK),          PARAMETER     :: EF1_BVOC    =  500.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_BVOC    =  500.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_BVOC    =  500.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_BVOC    =  500.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_BVOC    =  500.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_BVOC    =  500.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_BVOC    =  500.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_BVOC    =  500.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_BVOC    =  500.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_BVOC   =  500.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_BVOC   =  500.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_BVOC   =  80.0_rk      ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_BVOC   =  80.0_rk      ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_BVOC   =  80.0_rk      ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_BVOC   =  80.0_rk      ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for BIDI VOC species as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_BVOC  = 1.0_rk
        REAL(RK),          PARAMETER     :: AGRO_BVOC  = 1.0_rk
        REAL(RK),          PARAMETER     :: AMAT_BVOC  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_BVOC  = 1.0_rk

! Plant-Dependent emissions capacity/factors (EFs) for Stress VOCs (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! Stress VOC (15 compounds): See Table 1 of Guenther et al. (2012)
        REAL(RK),          PARAMETER     :: EF1_SVOC    =  300.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_SVOC    =  300.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_SVOC    =  300.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_SVOC    =  300.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_SVOC    =  300.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_SVOC    =  300.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_SVOC    =  300.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_SVOC    =  300.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_SVOC    =  300.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_SVOC   =  300.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_SVOC   =  300.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_SVOC   =  300.0_rk     ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_SVOC   =  300.0_rk     ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_SVOC   =  300.0_rk     ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_SVOC   =  300.0_rk     ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Stress VOCs as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_SVOC  = 1.0_rk
        REAL(RK),          PARAMETER     :: AGRO_SVOC  = 1.0_rk
        REAL(RK),          PARAMETER     :: AMAT_SVOC  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_SVOC  = 1.0_rk

! Plant-Dependent emissions capacity/factors (EFs) for Other VOCs (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! Other VOC (49 compounds): See Table 1 of Guenther et al. (2012)
        REAL(RK),          PARAMETER     :: EF1_OVOC    =  140.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_OVOC    =  140.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_OVOC    =  140.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_OVOC    =  140.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_OVOC    =  140.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_OVOC    =  140.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_OVOC    =  140.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_OVOC    =  140.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_OVOC    =  140.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_OVOC   =  140.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_OVOC   =  140.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_OVOC   =  140.0_rk     ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_OVOC   =  140.0_rk     ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_OVOC   =  140.0_rk     ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_OVOC   =  140.0_rk     ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Other VOCs as per Table 4 of Guenther et al., 2012
        REAL(RK),          PARAMETER     :: ANEW_OVOC  = 1.0_rk
        REAL(RK),          PARAMETER     :: AGRO_OVOC  = 1.0_rk
        REAL(RK),          PARAMETER     :: AMAT_OVOC  = 1.0_rk
        REAL(RK),          PARAMETER     :: AOLD_OVOC  = 1.0_rk

! Species-Dependent Parameterized Canopy Model Parameters (Table 4 of Guenther et al., 2012)
        REAL(RK),          PARAMETER     :: CT1_ISOP         =  95.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_ISOP         =  2.0_rk     !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_ISOP         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_ISOP         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_ISOP        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_ISOP         =  1.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_ISOP         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_ISOP        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_ISOP         =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_ISOP         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_ISOP        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_ISOP         =  1.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_ISOP         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_ISOP        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_MYRC         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_MYRC         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_MYRC         =  5.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_MYRC         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_MYRC        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_MYRC         =  5.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_MYRC         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_MYRC        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_MYRC         =  5.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_MYRC         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_MYRC        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_MYRC         =  5.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_MYRC         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_MYRC        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_SABI         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_SABI         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_SABI         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_SABI         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_SABI        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_SABI         =  1.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_SABI         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_SABI        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_SABI         =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_SABI         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_SABI        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_SABI         =  1.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_SABI         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_SABI        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_LIMO         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_LIMO         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_LIMO         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_LIMO         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_LIMO        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_LIMO         =  1.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_LIMO         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_LIMO        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_LIMO         =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_LIMO         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_LIMO        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_LIMO         =  5.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_LIMO         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_LIMO        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_CARE         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_CARE         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_CARE         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_CARE         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_CARE        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_CARE         =  1.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_CARE         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_CARE        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_CARE         =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_CARE         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_CARE        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_CARE         =  5.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_CARE         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_CARE        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_OCIM         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_OCIM         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_OCIM         =  5.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_OCIM         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_OCIM        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_OCIM         =  5.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_OCIM         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_OCIM        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_OCIM         =  5.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_OCIM         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_OCIM        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_OCIM         =  5.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_OCIM         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_OCIM        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_BPIN         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_BPIN         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_BPIN         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_BPIN         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_BPIN        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_BPIN         =  1.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_BPIN         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_BPIN        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_BPIN         =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_BPIN         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_BPIN        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_BPIN         =  5.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_BPIN         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_BPIN        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_APIN         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_APIN         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_APIN         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_APIN         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_APIN        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_APIN         =  1.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_APIN         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_APIN        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_APIN         =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_APIN         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_APIN        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_APIN         =  5.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_APIN         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_APIN        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_MONO         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_MONO         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_MONO         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_MONO         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_MONO        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_MONO         =  1.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_MONO         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_MONO        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_MONO         =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_MONO         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_MONO        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_MONO         =  5.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_MONO         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_MONO        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_FARN         =  130.0_rk   !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_FARN         =  2.37_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_FARN         =  5.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_FARN         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_FARN        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_FARN         =  5.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_FARN         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_FARN        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_FARN         =  5.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_FARN         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_FARN        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_FARN         =  5.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_FARN         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_FARN        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_CARY         =  130.0_rk   !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_CARY         =  2.37_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_CARY         =  5.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_CARY         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_CARY        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_CARY         =  5.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_CARY         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_CARY        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_CARY         =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_CARY         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_CARY        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_CARY         =  5.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_CARY         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_CARY        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_SESQ         =  130.0_rk   !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_SESQ         =  2.37_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_SESQ         =  5.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_SESQ         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_SESQ        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_SESQ         =  5.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_SESQ         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_SESQ        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_SESQ         =  5.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_SESQ         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_SESQ        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_SESQ         =  5.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_SESQ         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_SESQ        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_MBOL         =  95.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_MBOL         =  2.0_rk     !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_MBOL         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_MBOL         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_MBOL        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_MBOL         =  1.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_MBOL         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_MBOL        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_MBOL         =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_MBOL         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_MBOL        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_MBOL         =  1.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_MBOL         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_MBOL        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_METH         =  60.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_METH         =  1.6_rk     !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_METH         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_METH         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_METH        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_METH         =  1.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_METH         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_METH        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_METH         =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_METH         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_METH        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_METH         =  1.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_METH         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_METH        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_ACET         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_ACET         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_ACET         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_ACET         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_ACET        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_ACET         =  1.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_ACET         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_ACET        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_ACET         =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_ACET         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_ACET        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_ACET         =  1.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_ACET         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_ACET        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_CO           =  60.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_CO           =  1.6_rk     !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_CO           =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_CO           =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_CO          =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_CO           =  1.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_CO           =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_CO          =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_CO           =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_CO           =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_CO          =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_CO           =  1.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_CO           =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_CO          =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_BVOC         =  95.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_BVOC         =  2.0_rk     !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_BVOC         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_BVOC         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_BVOC        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_BVOC         =  1.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_BVOC         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_BVOC        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_BVOC         =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_BVOC         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_BVOC        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_BVOC         =  1.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_BVOC         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_BVOC        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_SVOC         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_SVOC         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_SVOC         =  5.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_SVOC         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_SVOC        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_SVOC         =  5.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_SVOC         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_SVOC        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_SVOC         =  5.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_SVOC         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_SVOC        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_SVOC         =  5.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_SVOC         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_SVOC        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(RK),          PARAMETER     :: CT1_OVOC         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_OVOC         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CAQ_OVOC         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(RK),          PARAMETER     :: TAQ_OVOC         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: DTAQ_OVOC        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(RK),          PARAMETER     :: CHT_OVOC         =  1.0_rk     !coefficient for high temperature stress
        REAL(RK),          PARAMETER     :: THT_OVOC         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: DTHT_OVOC        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(RK),          PARAMETER     :: CLT_OVOC         =  1.0_rk     !coefficient for low temperature stress
        REAL(RK),          PARAMETER     :: TLT_OVOC         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: DTLT_OVOC        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(RK),          PARAMETER     :: CHW_OVOC         =  1.0_rk     !coefficient for high wind stress
        REAL(RK),          PARAMETER     :: THW_OVOC         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(RK),          PARAMETER     :: DTHW_OVOC        =  8.0_rk     !delta threshold for high wind stress (m/s)

! Set tree and species dependent coefficients
        if (EMI_IND .eq. 1 ) then
            CT1  = CT1_ISOP
            CEO  = CEO_ISOP
            EF1  = EF1_ISOP
            EF2  = EF2_ISOP
            EF3  = EF3_ISOP
            EF4  = EF4_ISOP
            EF5  = EF5_ISOP
            EF6  = EF6_ISOP
            EF7  = EF7_ISOP
            EF8  = EF8_ISOP
            EF9  = EF9_ISOP
            EF10 = EF10_ISOP
            EF11 = EF11_ISOP
            EF12 = EF12_ISOP
            EF13 = EF13_ISOP
            EF14 = EF14_ISOP
            EF15 = EF15_ISOP
            ANEW  = ANEW_ISOP
            AGRO  = AGRO_ISOP
            AMAT  = AMAT_ISOP
            AOLD  = AOLD_ISOP
            CAQ   = CAQ_ISOP
            TAQ   = TAQ_ISOP
            DTAQ  = DTAQ_ISOP
            CHT   = CHT_ISOP
            THT   = THT_ISOP
            DTHT  = DTHT_ISOP
            CLT   = CLT_ISOP
            TLT   = TLT_ISOP
            DTLT  = DTLT_ISOP
            CHW   = CHW_ISOP
            THW   = THW_ISOP
            DTHW  = DTHW_ISOP
        else if (EMI_IND .eq. 2 ) then
            CT1 = CT1_MYRC
            CEO = CEO_MYRC
            EF1  = EF1_MYRC
            EF2  = EF2_MYRC
            EF3  = EF3_MYRC
            EF4  = EF4_MYRC
            EF5  = EF5_MYRC
            EF6  = EF6_MYRC
            EF7  = EF7_MYRC
            EF8  = EF8_MYRC
            EF9  = EF9_MYRC
            EF10 = EF10_MYRC
            EF11 = EF11_MYRC
            EF12 = EF12_MYRC
            EF13 = EF13_MYRC
            EF14 = EF14_MYRC
            EF15 = EF15_MYRC
            ANEW  = ANEW_MYRC
            AGRO  = AGRO_MYRC
            AMAT  = AMAT_MYRC
            AOLD  = AOLD_MYRC
            CAQ   = CAQ_MYRC
            TAQ   = TAQ_MYRC
            DTAQ  = DTAQ_MYRC
            CHT   = CHT_MYRC
            THT   = THT_MYRC
            DTHT  = DTHT_MYRC
            CLT   = CLT_MYRC
            TLT   = TLT_MYRC
            DTLT  = DTLT_MYRC
            CHW   = CHW_MYRC
            THW   = THW_MYRC
            DTHW  = DTHW_MYRC
        else if (EMI_IND .eq. 3 ) then
            CT1 = CT1_SABI
            CEO = CEO_SABI
            EF1  = EF1_SABI
            EF2  = EF2_SABI
            EF3  = EF3_SABI
            EF4  = EF4_SABI
            EF5  = EF5_SABI
            EF6  = EF6_SABI
            EF7  = EF7_SABI
            EF8  = EF8_SABI
            EF9  = EF9_SABI
            EF10 = EF10_SABI
            EF11 = EF11_SABI
            EF12 = EF12_SABI
            EF13 = EF13_SABI
            EF14 = EF14_SABI
            EF15 = EF15_SABI
            ANEW  = ANEW_SABI
            AGRO  = AGRO_SABI
            AMAT  = AMAT_SABI
            AOLD  = AOLD_SABI
            CAQ   = CAQ_SABI
            TAQ   = TAQ_SABI
            DTAQ  = DTAQ_SABI
            CHT   = CHT_SABI
            THT   = THT_SABI
            DTHT  = DTHT_SABI
            CLT   = CLT_SABI
            TLT   = TLT_SABI
            DTLT  = DTLT_SABI
            CHW   = CHW_SABI
            THW   = THW_SABI
            DTHW  = DTHW_SABI
        else if (EMI_IND .eq. 4 ) then
            CT1 = CT1_LIMO
            CEO = CEO_LIMO
            EF1  = EF1_LIMO
            EF2  = EF2_LIMO
            EF3  = EF3_LIMO
            EF4  = EF4_LIMO
            EF5  = EF5_LIMO
            EF6  = EF6_LIMO
            EF7  = EF7_LIMO
            EF8  = EF8_LIMO
            EF9  = EF9_LIMO
            EF10 = EF10_LIMO
            EF11 = EF11_LIMO
            EF12 = EF12_LIMO
            EF13 = EF13_LIMO
            EF14 = EF14_LIMO
            EF15 = EF15_LIMO
            ANEW  = ANEW_LIMO
            AGRO  = AGRO_LIMO
            AMAT  = AMAT_LIMO
            AOLD  = AOLD_LIMO
            CAQ   = CAQ_LIMO
            TAQ   = TAQ_LIMO
            DTAQ  = DTAQ_LIMO
            CHT   = CHT_LIMO
            THT   = THT_LIMO
            DTHT  = DTHT_LIMO
            CLT   = CLT_LIMO
            TLT   = TLT_LIMO
            DTLT  = DTLT_LIMO
            CHW   = CHW_LIMO
            THW   = THW_LIMO
            DTHW  = DTHW_LIMO
        else if (EMI_IND .eq. 5 ) then
            CT1 = CT1_CARE
            CEO = CEO_CARE
            EF1  = EF1_CARE
            EF2  = EF2_CARE
            EF3  = EF3_CARE
            EF4  = EF4_CARE
            EF5  = EF5_CARE
            EF6  = EF6_CARE
            EF7  = EF7_CARE
            EF8  = EF8_CARE
            EF9  = EF9_CARE
            EF10 = EF10_CARE
            EF11 = EF11_CARE
            EF12 = EF12_CARE
            EF13 = EF13_CARE
            EF14 = EF14_CARE
            EF15 = EF15_CARE
            ANEW  = ANEW_CARE
            AGRO  = AGRO_CARE
            AMAT  = AMAT_CARE
            AOLD  = AOLD_CARE
            CAQ   = CAQ_CARE
            TAQ   = TAQ_CARE
            DTAQ  = DTAQ_CARE
            CHT   = CHT_CARE
            THT   = THT_CARE
            DTHT  = DTHT_CARE
            CLT   = CLT_CARE
            TLT   = TLT_CARE
            DTLT  = DTLT_CARE
            CHW   = CHW_CARE
            THW   = THW_CARE
            DTHW  = DTHW_CARE
        else if (EMI_IND .eq. 6 ) then
            CT1 = CT1_OCIM
            CEO = CEO_OCIM
            EF1  = EF1_OCIM
            EF2  = EF2_OCIM
            EF3  = EF3_OCIM
            EF4  = EF4_OCIM
            EF5  = EF5_OCIM
            EF6  = EF6_OCIM
            EF7  = EF7_OCIM
            EF8  = EF8_OCIM
            EF9  = EF9_OCIM
            EF10 = EF10_OCIM
            EF11 = EF11_OCIM
            EF12 = EF12_OCIM
            EF13 = EF13_OCIM
            EF14 = EF14_OCIM
            EF15 = EF15_OCIM
            ANEW  = ANEW_OCIM
            AGRO  = AGRO_OCIM
            AMAT  = AMAT_OCIM
            AOLD  = AOLD_OCIM
            CAQ   = CAQ_OCIM
            TAQ   = TAQ_OCIM
            DTAQ  = DTAQ_OCIM
            CHT   = CHT_OCIM
            THT   = THT_OCIM
            DTHT  = DTHT_OCIM
            CLT   = CLT_OCIM
            TLT   = TLT_OCIM
            DTLT  = DTLT_OCIM
            CHW   = CHW_OCIM
            THW   = THW_OCIM
            DTHW  = DTHW_OCIM
        else if (EMI_IND .eq. 7 ) then
            CT1 = CT1_BPIN
            CEO = CEO_BPIN
            EF1  = EF1_BPIN
            EF2  = EF2_BPIN
            EF3  = EF3_BPIN
            EF4  = EF4_BPIN
            EF5  = EF5_BPIN
            EF6  = EF6_BPIN
            EF7  = EF7_BPIN
            EF8  = EF8_BPIN
            EF9  = EF9_BPIN
            EF10 = EF10_BPIN
            EF11 = EF11_BPIN
            EF12 = EF12_BPIN
            EF13 = EF13_BPIN
            EF14 = EF14_BPIN
            EF15 = EF15_BPIN
            ANEW  = ANEW_BPIN
            AGRO  = AGRO_BPIN
            AMAT  = AMAT_BPIN
            AOLD  = AOLD_BPIN
            CAQ   = CAQ_BPIN
            TAQ   = TAQ_BPIN
            DTAQ  = DTAQ_BPIN
            CHT   = CHT_BPIN
            THT   = THT_BPIN
            DTHT  = DTHT_BPIN
            CLT   = CLT_BPIN
            TLT   = TLT_BPIN
            DTLT  = DTLT_BPIN
            CHW   = CHW_BPIN
            THW   = THW_BPIN
            DTHW  = DTHW_BPIN
        else if (EMI_IND .eq. 8 ) then
            CT1 = CT1_APIN
            CEO = CEO_APIN
            EF1  = EF1_APIN
            EF2  = EF2_APIN
            EF3  = EF3_APIN
            EF4  = EF4_APIN
            EF5  = EF5_APIN
            EF6  = EF6_APIN
            EF7  = EF7_APIN
            EF8  = EF8_APIN
            EF9  = EF9_APIN
            EF10 = EF10_APIN
            EF11 = EF11_APIN
            EF12 = EF12_APIN
            EF13 = EF13_APIN
            EF14 = EF14_APIN
            EF15 = EF15_APIN
            ANEW  = ANEW_APIN
            AGRO  = AGRO_APIN
            AMAT  = AMAT_APIN
            AOLD  = AOLD_APIN
            CAQ   = CAQ_APIN
            TAQ   = TAQ_APIN
            DTAQ  = DTAQ_APIN
            CHT   = CHT_APIN
            THT   = THT_APIN
            DTHT  = DTHT_APIN
            CLT   = CLT_APIN
            TLT   = TLT_APIN
            DTLT  = DTLT_APIN
            CHW   = CHW_APIN
            THW   = THW_APIN
            DTHW  = DTHW_APIN
        else if (EMI_IND .eq. 9 ) then
            CT1 = CT1_MONO
            CEO = CEO_MONO
            EF1  = EF1_MONO
            EF2  = EF2_MONO
            EF3  = EF3_MONO
            EF4  = EF4_MONO
            EF5  = EF5_MONO
            EF6  = EF6_MONO
            EF7  = EF7_MONO
            EF8  = EF8_MONO
            EF9  = EF9_MONO
            EF10 = EF10_MONO
            EF11 = EF11_MONO
            EF12 = EF12_MONO
            EF13 = EF13_MONO
            EF14 = EF14_MONO
            EF15 = EF15_MONO
            ANEW  = ANEW_MONO
            AGRO  = AGRO_MONO
            AMAT  = AMAT_MONO
            AOLD  = AOLD_MONO
            CAQ   = CAQ_MONO
            TAQ   = TAQ_MONO
            DTAQ  = DTAQ_MONO
            CHT   = CHT_MONO
            THT   = THT_MONO
            DTHT  = DTHT_MONO
            CLT   = CLT_MONO
            TLT   = TLT_MONO
            DTLT  = DTLT_MONO
            CHW   = CHW_MONO
            THW   = THW_MONO
            DTHW  = DTHW_MONO
        else if (EMI_IND .eq. 10 ) then
            CT1 = CT1_FARN
            CEO = CEO_FARN
            EF1  = EF1_FARN
            EF2  = EF2_FARN
            EF3  = EF3_FARN
            EF4  = EF4_FARN
            EF5  = EF5_FARN
            EF6  = EF6_FARN
            EF7  = EF7_FARN
            EF8  = EF8_FARN
            EF9  = EF9_FARN
            EF10 = EF10_FARN
            EF11 = EF11_FARN
            EF12 = EF12_FARN
            EF13 = EF13_FARN
            EF14 = EF14_FARN
            EF15 = EF15_FARN
            ANEW  = ANEW_FARN
            AGRO  = AGRO_FARN
            AMAT  = AMAT_FARN
            AOLD  = AOLD_FARN
            CAQ   = CAQ_FARN
            TAQ   = TAQ_FARN
            DTAQ  = DTAQ_FARN
            CHT   = CHT_FARN
            THT   = THT_FARN
            DTHT  = DTHT_FARN
            CLT   = CLT_FARN
            TLT   = TLT_FARN
            DTLT  = DTLT_FARN
            CHW   = CHW_FARN
            THW   = THW_FARN
            DTHW  = DTHW_FARN
        else if (EMI_IND .eq. 11 ) then
            CT1 = CT1_CARY
            CEO = CEO_CARY
            EF1  = EF1_CARY
            EF2  = EF2_CARY
            EF3  = EF3_CARY
            EF4  = EF4_CARY
            EF5  = EF5_CARY
            EF6  = EF6_CARY
            EF7  = EF7_CARY
            EF8  = EF8_CARY
            EF9  = EF9_CARY
            EF10 = EF10_CARY
            EF11 = EF11_CARY
            EF12 = EF12_CARY
            EF13 = EF13_CARY
            EF14 = EF14_CARY
            EF15 = EF15_CARY
            ANEW  = ANEW_CARY
            AGRO  = AGRO_CARY
            AMAT  = AMAT_CARY
            AOLD  = AOLD_CARY
            CAQ   = CAQ_CARY
            TAQ   = TAQ_CARY
            DTAQ  = DTAQ_CARY
            CHT   = CHT_CARY
            THT   = THT_CARY
            DTHT  = DTHT_CARY
            CLT   = CLT_CARY
            TLT   = TLT_CARY
            DTLT  = DTLT_CARY
            CHW   = CHW_CARY
            THW   = THW_CARY
            DTHW  = DTHW_CARY
        else if (EMI_IND .eq. 12 ) then
            CT1 = CT1_SESQ
            CEO = CEO_SESQ
            EF1  = EF1_SESQ
            EF2  = EF2_SESQ
            EF3  = EF3_SESQ
            EF4  = EF4_SESQ
            EF5  = EF5_SESQ
            EF6  = EF6_SESQ
            EF7  = EF7_SESQ
            EF8  = EF8_SESQ
            EF9  = EF9_SESQ
            EF10 = EF10_SESQ
            EF11 = EF11_SESQ
            EF12 = EF12_SESQ
            EF13 = EF13_SESQ
            EF14 = EF14_SESQ
            EF15 = EF15_SESQ
            ANEW  = ANEW_SESQ
            AGRO  = AGRO_SESQ
            AMAT  = AMAT_SESQ
            AOLD  = AOLD_SESQ
            CAQ   = CAQ_SESQ
            TAQ   = TAQ_SESQ
            DTAQ  = DTAQ_SESQ
            CHT   = CHT_SESQ
            THT   = THT_SESQ
            DTHT  = DTHT_SESQ
            CLT   = CLT_SESQ
            TLT   = TLT_SESQ
            DTLT  = DTLT_SESQ
            CHW   = CHW_SESQ
            THW   = THW_SESQ
            DTHW  = DTHW_SESQ
        else if (EMI_IND .eq. 13 ) then
            CT1 = CT1_MBOL
            CEO = CEO_MBOL
            EF1  = EF1_MBOL
            EF2  = EF2_MBOL
            EF3  = EF3_MBOL
            EF4  = EF4_MBOL
            EF5  = EF5_MBOL
            EF6  = EF6_MBOL
            EF7  = EF7_MBOL
            EF8  = EF8_MBOL
            EF9  = EF9_MBOL
            EF10 = EF10_MBOL
            EF11 = EF11_MBOL
            EF12 = EF12_MBOL
            EF13 = EF13_MBOL
            EF14 = EF14_MBOL
            EF15 = EF15_MBOL
            ANEW  = ANEW_MBOL
            AGRO  = AGRO_MBOL
            AMAT  = AMAT_MBOL
            AOLD  = AOLD_MBOL
            CAQ   = CAQ_MBOL
            TAQ   = TAQ_MBOL
            DTAQ  = DTAQ_MBOL
            CHT   = CHT_MBOL
            THT   = THT_MBOL
            DTHT  = DTHT_MBOL
            CLT   = CLT_MBOL
            TLT   = TLT_MBOL
            DTLT  = DTLT_MBOL
            CHW   = CHW_MBOL
            THW   = THW_MBOL
            DTHW  = DTHW_MBOL
        else if (EMI_IND .eq. 14 ) then
            CT1 = CT1_METH
            CEO = CEO_METH
            EF1  = EF1_METH
            EF2  = EF2_METH
            EF3  = EF3_METH
            EF4  = EF4_METH
            EF5  = EF5_METH
            EF6  = EF6_METH
            EF7  = EF7_METH
            EF8  = EF8_METH
            EF9  = EF9_METH
            EF10 = EF10_METH
            EF11 = EF11_METH
            EF12 = EF12_METH
            EF13 = EF13_METH
            EF14 = EF14_METH
            EF15 = EF15_METH
            ANEW  = ANEW_METH
            AGRO  = AGRO_METH
            AMAT  = AMAT_METH
            AOLD  = AOLD_METH
            CAQ   = CAQ_METH
            TAQ   = TAQ_METH
            DTAQ  = DTAQ_METH
            CHT   = CHT_METH
            THT   = THT_METH
            DTHT  = DTHT_METH
            CLT   = CLT_METH
            TLT   = TLT_METH
            DTLT  = DTLT_METH
            CHW   = CHW_METH
            THW   = THW_METH
            DTHW  = DTHW_METH
        else if (EMI_IND .eq. 15 ) then
            CT1 = CT1_ACET
            CEO = CEO_ACET
            EF1  = EF1_ACET
            EF2  = EF2_ACET
            EF3  = EF3_ACET
            EF4  = EF4_ACET
            EF5  = EF5_ACET
            EF6  = EF6_ACET
            EF7  = EF7_ACET
            EF8  = EF8_ACET
            EF9  = EF9_ACET
            EF10 = EF10_ACET
            EF11 = EF11_ACET
            EF12 = EF12_ACET
            EF13 = EF13_ACET
            EF14 = EF14_ACET
            EF15 = EF15_ACET
            ANEW  = ANEW_ACET
            AGRO  = AGRO_ACET
            AMAT  = AMAT_ACET
            AOLD  = AOLD_ACET
            CAQ   = CAQ_ACET
            TAQ   = TAQ_ACET
            DTAQ  = DTAQ_ACET
            CHT   = CHT_ACET
            THT   = THT_ACET
            DTHT  = DTHT_ACET
            CLT   = CLT_ACET
            TLT   = TLT_ACET
            DTLT  = DTLT_ACET
            CHW   = CHW_ACET
            THW   = THW_ACET
            DTHW  = DTHW_ACET
        else if (EMI_IND .eq. 16 ) then
            CT1 = CT1_CO
            CEO = CEO_CO
            EF1  = EF1_CO
            EF2  = EF2_CO
            EF3  = EF3_CO
            EF4  = EF4_CO
            EF5  = EF5_CO
            EF6  = EF6_CO
            EF7  = EF7_CO
            EF8  = EF8_CO
            EF9  = EF9_CO
            EF10 = EF10_CO
            EF11 = EF11_CO
            EF12 = EF12_CO
            EF13 = EF13_CO
            EF14 = EF14_CO
            EF15 = EF15_CO
            ANEW  = ANEW_CO
            AGRO  = AGRO_CO
            AMAT  = AMAT_CO
            AOLD  = AOLD_CO
            CAQ   = CAQ_CO
            TAQ   = TAQ_CO
            DTAQ  = DTAQ_CO
            CHT   = CHT_CO
            THT   = THT_CO
            DTHT  = DTHT_CO
            CLT   = CLT_CO
            TLT   = TLT_CO
            DTLT  = DTLT_CO
            CHW   = CHW_CO
            THW   = THW_CO
            DTHW  = DTHW_CO
        else if (EMI_IND .eq. 17 ) then
            CT1 = CT1_BVOC
            CEO = CEO_BVOC
            EF1  = EF1_BVOC
            EF2  = EF2_BVOC
            EF3  = EF3_BVOC
            EF4  = EF4_BVOC
            EF5  = EF5_BVOC
            EF6  = EF6_BVOC
            EF7  = EF7_BVOC
            EF8  = EF8_BVOC
            EF9  = EF9_BVOC
            EF10 = EF10_BVOC
            EF11 = EF11_BVOC
            EF12 = EF12_BVOC
            EF13 = EF13_BVOC
            EF14 = EF14_BVOC
            EF15 = EF15_BVOC
            ANEW  = ANEW_BVOC
            AGRO  = AGRO_BVOC
            AMAT  = AMAT_BVOC
            AOLD  = AOLD_BVOC
            CAQ   = CAQ_BVOC
            TAQ   = TAQ_BVOC
            DTAQ  = DTAQ_BVOC
            CHT   = CHT_BVOC
            THT   = THT_BVOC
            DTHT  = DTHT_BVOC
            CLT   = CLT_BVOC
            TLT   = TLT_BVOC
            DTLT  = DTLT_BVOC
            CHW   = CHW_BVOC
            THW   = THW_BVOC
            DTHW  = DTHW_BVOC
        else if (EMI_IND .eq. 18 ) then
            CT1 = CT1_SVOC
            CEO = CEO_SVOC
            EF1  = EF1_SVOC
            EF2  = EF2_SVOC
            EF3  = EF3_SVOC
            EF4  = EF4_SVOC
            EF5  = EF5_SVOC
            EF6  = EF6_SVOC
            EF7  = EF7_SVOC
            EF8  = EF8_SVOC
            EF9  = EF9_SVOC
            EF10 = EF10_SVOC
            EF11 = EF11_SVOC
            EF12 = EF12_SVOC
            EF13 = EF13_SVOC
            EF14 = EF14_SVOC
            EF15 = EF15_SVOC
            ANEW  = ANEW_SVOC
            AGRO  = AGRO_SVOC
            AMAT  = AMAT_SVOC
            AOLD  = AOLD_SVOC
            CAQ   = CAQ_SVOC
            TAQ   = TAQ_SVOC
            DTAQ  = DTAQ_SVOC
            CHT   = CHT_SVOC
            THT   = THT_SVOC
            DTHT  = DTHT_SVOC
            CLT   = CLT_SVOC
            TLT   = TLT_SVOC
            DTLT  = DTLT_SVOC
            CHW   = CHW_SVOC
            THW   = THW_SVOC
            DTHW  = DTHW_SVOC
        else   ! EMI_IND = 19
            CT1 = CT1_OVOC
            CEO = CEO_OVOC
            EF1  = EF1_OVOC
            EF2  = EF2_OVOC
            EF3  = EF3_OVOC
            EF4  = EF4_OVOC
            EF5  = EF5_OVOC
            EF6  = EF6_OVOC
            EF7  = EF7_OVOC
            EF8  = EF8_OVOC
            EF9  = EF9_OVOC
            EF10 = EF10_OVOC
            EF11 = EF11_OVOC
            EF12 = EF12_OVOC
            EF13 = EF13_OVOC
            EF14 = EF14_OVOC
            EF15 = EF15_OVOC
            ANEW  = ANEW_OVOC
            AGRO  = AGRO_OVOC
            AMAT  = AMAT_OVOC
            AOLD  = AOLD_OVOC
            CAQ   = CAQ_OVOC
            TAQ   = TAQ_OVOC
            DTAQ  = DTAQ_OVOC
            CHT   = CHT_OVOC
            THT   = THT_OVOC
            DTHT  = DTHT_OVOC
            CLT   = CLT_OVOC
            TLT   = TLT_OVOC
            DTLT  = DTLT_OVOC
            CHW   = CHW_OVOC
            THW   = THW_OVOC
            DTHW  = DTHW_OVOC
        end if

        if (LU_OPT .eq. 0 .or. LU_OPT .eq. 1) then !VIIRS or MODIS  LU types

! Simple MEGAN (Table 3 in Guenther et al., 2012) PFT to VIIRS/MODIS VTYPE mapping
            if (VTYPE .eq. 1) then !VIIRS Cat 1 Evergreen Needleleaf
                !--> Average Needleleaf Evergreen Temperate Tree and Needleleaf Evergreen Boreal Tree

                EF = (EF1+EF2)/2.0_rk
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                ROOTA = 6.706_rk
                ROOTB = 2.175_rk

            else if (VTYPE .eq. 2) then !VIIRS/MODIS Cat 2 Evergreen Broadleaf
                !--> Average Broadleaf Evergreen Tropical Tree and Broadleaf Evergreen Temperate Tree

                EF = (EF4+EF5)/2.0_rk
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                ROOTA = 7.344_rk
                ROOTB = 1.303_rk

            else if (VTYPE .eq. 3) then !VIIRS/MODIS Cat 3 Deciduous Needleaf
                !--> Average Needleleaf Deciduous Boreal Tree

                EF = EF3
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                ROOTA = 7.066_rk
                ROOTB = 1.953_rk

            else if (VTYPE .eq. 4) then !VIIRS/MODIS Cat 4 Deciduous Broadleaf
                !--> Average Broadleaf Deciduous Tropical Tree, Broadleaf Deciduous Temperate Tree,
                ! and Broadleaf Deciduous Boreal Tree

                EF = (EF6+EF7+EF8)/3.0_rk
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                ROOTA = 5.990_rk
                ROOTB = 1.955_rk

            else if (VTYPE .eq. 5) then !VIIRS/MODIS Cat 5 Mixed forests
                !--> Avearge of all above EF1-EF8 PFTs.

                EF = (EF1+EF2+EF3+EF4+EF5+EF6+EF7+EF8)/8.0_rk
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                ROOTA = 4.453_rk
                ROOTB = 1.631_rk

            else if (VTYPE .ge. 6 .and. VTYPE .le. 7) then !VIIRS/MODIS Cat 6-7 Closed/Open Shrublands
                !--> Avearge Broadleaf Evergreen Temperate Shrub, Broadleaf Deciduous Temperate Shrub,
                ! and Broadleaf Deciduous Boreal Shrub

                EF = (EF9+EF10+EF11)/3.0_rk
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                ROOTA = (6.326_rk + 7.718_rk)/2.0_rk
                ROOTB = (1.567_rk + 1.262_rk)/2.0_rk

            else if (VTYPE .ge. 8 .and. VTYPE .le. 10) then !VIIRS/MODIS Cat 8-10 Savannas and Grasslands
                !--> Avearge Arctic C3 Grass, Cool C3 Grass, Warm C4 Grass)

                EF = (EF12+EF13+EF14)/3.0_rk
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                ROOTA = (7.604_rk + 8.235_rk + 10.74_rk)/3.0_rk
                ROOTB = (2.300_rk + 1.627_rk + 2.608_rk)/3.0_rk

            else if (VTYPE .eq. 12 ) then !VIIRS/MODIS Cat 12 Croplands
                !--> Crop1

                EF = EF15
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                ROOTA = 5.558_rk
                ROOTB = 2.614_rk

            else if (VTYPE .eq. 14 ) then !VIIRS/MODIS Cat 14 Croplands/Natural Mosaic
                !--> Crop1

                EF = EF15
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                ROOTA = 5.558_rk
                ROOTB = 2.614_rk

            else if (VTYPE .ge. 18 .and. VTYPE .le. 19) then !VIIRS/MODIS Cat 18-19 Wooded/Mixed Tundra
                !--> Crop1

                EF = (EF9+EF10+EF11)/3.0_rk  !Assume EF for wooded/mixed tundra = shrublands above
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                ROOTA = (6.326_rk + 7.718_rk)/2.0_rk
                ROOTB = (1.567_rk + 1.262_rk)/2.0_rk

            else

                EF = 0.0_rk  !set EF=0 for non-valid VTYPEs

            end if

        else
            write(*,*)  'Wrong LU_OPT choice of ', LU_OPT, 'in namelist, only VIIRS (0) or MODIS (1) available right now...exiting'
            call exit(2)
        end if

    END SUBROUTINE CANOPY_BIOP

end module canopy_bioparm_mod
