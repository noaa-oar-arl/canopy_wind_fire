module canopy_bioemi_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_BIO( ZK, FCLAI, FCH, LAI, FSUN, &
        PPFD_SUN, PPFD_SHADE, TLEAF_SUN, TLEAF_SHADE, PPFD24_SUN, &
        PPFD24_SHADE, TLEAF24_AVE,  &
        PPFD240_SUN, PPFD240_SHADE, &
        TLEAF240_AVE, TEMP2, LU_OPT, &
        VTYPE, MODRES, CCE, VERT, CO2OPT, CO2SET, &
        LEAFAGEOPT, PASTLAI, CURRENTLAI, TSTEPLAI, &
        LOSSOPT, LOSSSET, LOSSIND, LIFETIME, USTAR, &
        SOIMOPT, SOIM1, SOIM2, SOIM3, SOIM4, SOID1, SOID2, SOID3, &
        SOID4, WILT, AQOPT, W126_SET, W126_REF,  &
        MODLAYS, EMI_IND, EMI_OUT)

!-----------------------------------------------------------------------

! Description:
!     computes parameterized canopy biogenic emissions

! Preconditions:
!     in-canopy FCLAI, model LAI, etc.

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Clifton et al. (2022) algorithms
! Citation:
! Clifton, O. E. et al. (2022). Large eddy simulation for investigating
! coupled forest canopy and turbulence influences on atmospheric chemistry.
! Journal of Advances in Modeling Earth Systems, 14, e2022MS003078.
! https://doi.org/10.1029/2022MS003078
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Jan 2023 P.C. Campbell: Initial canopy isoprene only version
!     Feb 2023 P.C. Campbell: Modified for multiple biogenic species
!     Jul 2023 P.C. Campbell: Restructured to use FSUN, TLEAF, and PPFD
!                             as inputs
!     Sept 2023: QZ Rasool : Modifications for  LeafAge Response for multiple
!     BVOCs
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod,  ONLY: rk,rgasuniv   !constants for canopy models
        use canopy_utils_mod,  ONLY: interp_linear1_internal, &
            GET_GAMMA_CO2,GET_GAMMA_LEAFAGE, &
            GET_GAMMA_SOIM, GET_GAMMA_AQ, GET_CANLOSS_BIO
        use canopy_bioparm_mod
        use canopy_tleaf_mod

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )       :: ZK(:)              ! Model heights (m)
        REAL(RK),    INTENT( IN )       :: FCLAI(:)           ! Fractional (z) shapes of the
        ! plant surface distribution (nondimensional), i.e., a Fractional Culmulative LAI
        REAL(RK),    INTENT( IN )       :: FCH                ! Canopy height (m)
        REAL(RK),    INTENT( IN )       :: LAI                ! Total Leaf Area Index
        REAL(RK),    INTENT( IN )       :: FSUN(:)            ! Sunlit/Shaded fraction from photolysis correction factor
        REAL(RK),    INTENT( IN )       :: TLEAF_SUN(:)       ! Leaf temp for sunlit leaves (K)
        REAL(RK),    INTENT( IN )       :: TLEAF_SHADE(:)     ! Leaf temp for shaded leaves (K)
        REAL(RK),    INTENT( IN )       :: PPFD_SUN(:)        ! PPFD for sunlit leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: PPFD_SHADE(:)      ! PPFD for shaded leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: PPFD24_SUN(:)      ! PPFD for sunlit leaves (umol phot/m2 s) --24 hr ave
        REAL(RK),    INTENT( IN )       :: PPFD24_SHADE(:)    ! PPFD for shaded leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: TLEAF24_AVE(:)     ! Average Leaf temp (K)
        REAL(RK),    INTENT( IN )       :: PPFD240_SUN(:)     ! PPFD for sunlit leaves (umol phot/m2 s) -- 240 hr ave
        REAL(RK),    INTENT( IN )       :: PPFD240_SHADE(:)   ! PPFD for shaded leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: TLEAF240_AVE(:)    ! Average Leaf temp (K)

        REAL(RK),    INTENT( IN )       :: TEMP2           ! Model input 2-m Temperature (K)
        INTEGER,     INTENT( IN )       :: LU_OPT          ! integer for LU type from model mapped to Massman et al. (default = 0/VIIRS)
        INTEGER,     INTENT( IN )       :: VTYPE           ! Grid cell dominant vegetation type
        REAL(RK),    INTENT( IN )       :: MODRES          ! Canopy model input vertical resolution (m)
        REAL(RK),    INTENT( IN )       :: CCE             ! MEGAN Canopy environment coefficient.
        INTEGER,     INTENT( IN )       :: VERT            ! MEGAN vertical integration option (default = 0/no integration)
        INTEGER,     INTENT( IN )       :: CO2OPT          ! Option for co2 inhibition calculation
        REAL(RK),    INTENT( IN )       :: CO2SET          ! User set atmospheric CO2 conc [ppmv]
        INTEGER,     INTENT( IN )       :: SOIMOPT         ! Option for soil moisture factor
        REAL(RK),    INTENT( IN )       :: SOIM1           ! Volumetric soil moisture layer 1 [m3/m3]
        REAL(RK),    INTENT( IN )       :: SOIM2           ! Volumetric soil moisture layer 2 [m3/m3]
        REAL(RK),    INTENT( IN )       :: SOIM3           ! Volumetric soil moisture layer 3 [m3/m3]
        REAL(RK),    INTENT( IN )       :: SOIM4           ! Volumetric soil moisture layer 4 [m3/m3]
        REAL(RK),    INTENT( IN )       :: SOID1           ! Soil depth layer 1 [cm]
        REAL(RK),    INTENT( IN )       :: SOID2           ! Soil depth layer 2 [cm]
        REAL(RK),    INTENT( IN )       :: SOID3           ! Soil depth layer 3 [cm]
        REAL(RK),    INTENT( IN )       :: SOID4           ! Soil depth layer 4 [cm]
        REAL(RK),    INTENT( IN )       :: WILT            ! Wilting point [proportion]
        INTEGER,     INTENT( IN )       :: AQOPT           ! Option for aq stress calculation
        REAL(RK),    INTENT( IN )       :: W126_SET        ! User set ozone W126 [ppm-hours]
        REAL(RK),    INTENT( IN )       :: W126_REF        ! GFS calculated, ozone W126 [ppm-hours]

        INTEGER,    INTENT( IN )        :: LEAFAGEOPT     ! leafage_opt (0= ON, 1= off i.e. GAMMALEAFAGE =1, in canopy_readnml.F90)
        REAL(RK),    INTENT( IN )       :: PASTLAI        ! Past LAI [cm2/cm2]
        REAL(RK),    INTENT( IN )       :: CURRENTLAI     ! Current LAI [cm2/cm2]
        REAL(RK),    INTENT( IN )       :: TSTEPLAI       !Number of days between the past and current LAI

        INTEGER,     INTENT( IN )       :: MODLAYS         ! Input total model layers
        INTEGER,     INTENT( IN )       :: LOSSOPT         ! Option for canopy loss factor when summing top of canopy emissions
        REAL(RK),    INTENT( IN )       :: LIFETIME        ! Above canopy chemical lifetime of VOC (s)
        REAL(RK),    INTENT( IN )       :: LOSSSET         ! Input value for constant canopy loss factor applied used with loss_opt=2 (Default = 0.96)
        INTEGER,     INTENT( IN )       :: LOSSIND         ! Input integer for applying canopy loss factor to all species (=0) or only specific biogenics specie indices (> 0)

        REAL(RK),    INTENT( IN )       :: USTAR           ! Above canopy friction velocity (m/s)
        INTEGER,     INTENT( IN )       :: EMI_IND         ! Input biogenic emissions index
        REAL(RK),    INTENT( OUT )      :: EMI_OUT(:)      ! Output canopy layer volume emissions (kg m-3 s-1)

! Local Variables
        REAL(RK) :: GammaTLEAF_SUN_NUM(SIZE(ZK))   ! Numerator in Tleaf sun activity factor
        REAL(RK) :: GammaTLEAF_SHADE_NUM(SIZE(ZK)) ! Numerator in Tleaf shade activity factor
        REAL(RK) :: GammaTLEAF_SUN_DEN(SIZE(ZK))   ! Denominator in Tleaf sun activity factor
        REAL(RK) :: GammaTLEAF_SHADE_DEN(SIZE(ZK)) ! Denominator in Tleaf sun activity factor
        REAL(RK) :: GammaTLEAF_SUN(SIZE(ZK))       ! Tleaf sun activity factor
        REAL(RK) :: GammaTLEAF_SHADE(SIZE(ZK))     ! Tleaf shade activity factor
        REAL(RK) :: GammaTLEAF_AVE(SIZE(ZK))       ! Average Tleaf activity factor
        REAL(RK) :: CP_SUN(SIZE(ZK))               ! Normalized emission capacity sun at PPFD = 1000 umol phot/m2 s
        REAL(RK) :: CP_SHADE(SIZE(ZK))             ! Normalized emission capacity shade at PPFD = 1000 umol phot/m2 s
        REAL(RK) :: ALPHA_P_SUN(SIZE(ZK))          ! Quantum yield of isoprene sunlit (mol/mol)
        REAL(RK) :: ALPHA_P_SHADE(SIZE(ZK))        ! Quantum yield of isoprene shade (mol/mol)
        REAL(RK) :: GammaPPFD_SUN(SIZE(ZK))        ! PPFD activity factor sun (unitless)
        REAL(RK) :: GammaPPFD_SHADE(SIZE(ZK))      ! PPFD activity factor shade (unitless)
        REAL(RK) :: GammaPPFD_AVE(SIZE(ZK))        ! PPFD activity factor ave sun and shade
        REAL(RK) :: E_OPT(SIZE(ZK))                ! maximum normalized emission capacity
        REAL(RK) :: TLEAF_OPT(SIZE(ZK))            ! Tleaf at which E_OPT occurs (K)
        REAL(RK) :: FLAI(SIZE(ZK))                 ! Fractional LAI in layer
        REAL(RK) :: VPGWT(SIZE(ZK))                ! MEGANv3-like in-canopy weighting factor
        REAL(RK) :: GAUSS(SIZE(ZK))                ! MEGANv3-like in-canopy gaussian
        REAL(RK) :: CT1                            ! Activation energy (kJ/mol)
        REAL(RK) :: CEO                            ! Empirical coefficient
        REAL(RK) :: EF                             ! Final Mapped Emission factor (EF) (ug/m2 hr)


        REAL(RK) :: TABOVECANOPY !(SIZE(ZK))  ! Above Canopy Temp assigned = TEMP2 i.e., Model input 2-m Temperature (K for now)
        ! Empirical coeff.'s for Leaf Age factor calculations (see
        ! canopy_bioparm_mod or call canopy_biop)
        REAL(RK) :: ANEW
        REAL(RK) :: AGRO
        REAL(RK) :: AMAT
        REAL(RK) :: AOLD

        !Coeff.'s and threshold/delta threshold for air quality stress factors from canopy_biop
        REAL(RK) :: CAQ
        REAL(RK) :: TAQ  ![ppm-hours]
        REAL(RK) :: DTAQ ![ppm-hours]

        ! Coefficients A and B used for PFT dependent cumulative root depth fraction
        REAL(RK) :: ROOTA ! [m-1]
        REAL(RK) :: ROOTB ! [m-1]
        REAL(RK) :: GAMMASOIM ! Soil moisture factor

        REAL(RK) :: GAMMAAQ                        !Air quality stress factor

        REAL(RK) :: GAMMACO2                       ! CO2 inhibition factor (isoprene only)

        REAL(RK) :: GAMMALEAFAGE !(SIZE(ZK))                 ! LEAF AGE factor

        REAL(RK) :: CANLOSS_FAC                    !Canopy loss factor for summing option

        integer i, LAYERS

! Constant Canopy Parameters
        REAL(RK),          PARAMETER     :: PPFD0_SUN       =  200.0      !Constant PPFDo sunlit (umol/m2 s) (Guenther et al.,2012)
        REAL(RK),          PARAMETER     :: PPFD0_SHADE     =  50.0       !Constant PPFDo shaded (umol/m2 s) (Guenther et al.,2012)
        REAL(RK),          PARAMETER     :: CT2             =  230.0_rk   !Deactivation energy (kJ/mol) (Guenther et al., 2012)

! Calculate maximum normalized emission capacity (E_OPT) and Tleaf at E_OPT
        TLEAF_OPT = 313.0_rk + (0.6_rk * (TLEAF240_AVE-297.0_rk)) !Guenther et al. (2012)

! Calculate emission species/plant-dependent mapped emission factors and other important coefficients for gamma terms
        call canopy_biop(EMI_IND, LU_OPT, VTYPE, EF, CT1, CEO, ANEW, AGRO, AMAT, AOLD, ROOTA, ROOTB, CAQ, TAQ, DTAQ)

        E_OPT = CEO * EXP(0.05_rk * (TLEAF24_AVE-297.0_rk)) * EXP(0.05_rk * (TLEAF240_AVE-297.0_rk))

! Calculate gamma (activity) values for average Tleaf (Clifton et al., 2022; based on Guenther et al. 2012)
        GammaTLEAF_SUN_NUM = CT2*exp((CT1/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SUN)))
        GammaTLEAF_SUN_DEN = (CT2-CT1)*(1.0-exp((CT2/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SUN))))
        GammaTLEAF_SUN     = E_OPT*(GammaTLEAF_SUN_NUM/GammaTLEAF_SUN_DEN)

        GammaTLEAF_SHADE_NUM = CT2*exp((CT1/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SHADE)))
        GammaTLEAF_SHADE_DEN = (CT2-CT1)*(1.0-exp((CT2/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SHADE))))
        GammaTLEAF_SHADE     = E_OPT*(GammaTLEAF_SHADE_NUM/GammaTLEAF_SHADE_DEN)

        GammaTLEAF_AVE = (GammaTLEAF_SUN*FSUN) + (GammaTLEAF_SHADE*(1.0-FSUN)) ! average = sum sun and shade weighted by sunlit fraction

        GammaTLEAF_AVE = MAX( GammaTLEAF_AVE, 0.0_rk )

! Calculate gamma (activity) values for average PPFD (Clifton et al., 2022; based on Guenther et al. 2012)
        ALPHA_P_SUN = 0.004 - 0.0005*log(PPFD240_SUN)
        ALPHA_P_SHADE = 0.004 - 0.0005*log(PPFD240_SHADE)
        CP_SUN = 0.0468*(PPFD240_SUN**(0.6))*exp(0.005*(PPFD24_SUN-PPFD0_SUN))
        CP_SHADE = 0.0468*(PPFD240_SHADE**(0.6))*exp(0.005*(PPFD24_SHADE-PPFD0_SHADE))
        GammaPPFD_SUN   = CP_SUN*((ALPHA_P_SUN*PPFD_SUN)/SQRT(1.0 + (ALPHA_P_SUN**2.0) * (PPFD_SUN**2.0)))
        GammaPPFD_SHADE = CP_SHADE*((ALPHA_P_SHADE*PPFD_SHADE)/SQRT(1.0 + (ALPHA_P_SHADE**2.0) * (PPFD_SHADE**2.0)))

        GammaPPFD_AVE = (GammaPPFD_SUN*FSUN) + (GammaPPFD_SHADE*(1.0-FSUN)) ! average = sum sun and shade weighted by sunlit fraction

        GammaPPFD_AVE = MAX( GammaPPFD_AVE, 0.0_rk )

! Get CO2 inhibition factor for isoprene only

        if (EMI_IND .eq. 1) then  !Isoprene
            GAMMACO2 = GET_GAMMA_CO2(CO2OPT,CO2SET)
        else
            GAMMACO2 = 1.0_rk
        end if

! Get Soil Moisture Factor
        GAMMASOIM =  GET_GAMMA_SOIM(SOIMOPT,SOIM1,SOIM2,SOIM3,SOIM4,SOID1,SOID2,SOID3,SOID4,WILT, &
            ROOTA,ROOTB)

! Get LEAF AGE factor

        TABOVECANOPY  = TEMP2   !TEMP2 (above air temp) for TABOVECANOPY
        !do i=1, SIZE(ZK)
        GAMMALEAFAGE = GET_GAMMA_LEAFAGE(LEAFAGEOPT, PASTLAI, CURRENTLAI, TSTEPLAI, TABOVECANOPY, ANEW, AGRO, AMAT, AOLD)
        !end do

! Get AQ Stress Factor
!        print*, 'AQOPT=',AQOPT,'W126_REF=',W126_REF,'W126_SET=',W126_SET
!        print*, 'CAQ=',CAQ,'TAQ=',TAQ,'DTAQ=',DTAQ
        GAMMAAQ = GET_GAMMA_AQ(AQOPT, W126_REF, W126_SET, CAQ, TAQ, DTAQ)
!        print*,'GAMMAAQ=',GAMMAAQ

! Get canopy loss factor (only used in vertical summing options and empirical formulation and parameters based on isoprene)
! Note:  Allowed for other BVOCs but use caution when applying to compare with above canopy flux observations

        CANLOSS_FAC = 1.0_rk  !Initialize
        ! All species
        if (LOSSIND .eq. 0) then
            if (LOSSOPT .eq. 2) then !User set value from NL
                CANLOSS_FAC = LOSSSET
            else                !Try and calculate if turned on
                CANLOSS_FAC = GET_CANLOSS_BIO(LOSSOPT,LIFETIME,USTAR,FCH)
            end if
        end if

        !Only for a specific biogenic species/indice
        if (LOSSIND .eq. EMI_IND) then
            if (LOSSOPT .eq. 2) then !User set value from NL
                CANLOSS_FAC = LOSSSET
            else                !Try and calculate if turned on
                CANLOSS_FAC = GET_CANLOSS_BIO(LOSSOPT,LIFETIME,USTAR,FCH)
            end if
        end if

! Calculate emissions profile in the canopy
        EMI_OUT = 0.0_rk  ! set initial emissions profile to zero
        FLAI = 0.0_rk  ! set initial fractional FLAI (LAD) profile to zero

        if (VERT .eq. 0) then         !Full 3D leaf-level biogenic emissions (no averaging, summing, or integration)
            do i=1, SIZE(ZK)
                if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then           ! above ground level and at/below canopy top
                    if (i .lt. MODLAYS)  then
                        FLAI(i) = ((FCLAI(i+1) - FCLAI(i)) * LAI)/MODRES    !fractional LAI in each layer converted to LAD (m2 m-3)
                    else
                        FLAI(i) = FLAI(MODLAYS-1)
                    end if
                    EMI_OUT(i) = FLAI(i) * EF * GammaTLEAF_AVE(i) * GammaPPFD_AVE(i) * GAMMACO2 * CCE * &
                        GAMMALEAFAGE * GAMMASOIM  ! (ug m-3 hr-1)
                    EMI_OUT(i) = EMI_OUT(i) * 2.7777777777778E-13_rk    !convert emissions output to (kg m-3 s-1)
                end if
            end do
        else if (VERT .eq. 1) then       !"MEGANv3-like": Use weighting factors normalized to plant distribution shape (FCLAI)
            !across canopy layers
            LAYERS = min((floor(FCH/MODRES) + 1),MODLAYS)
            do i=1,  SIZE(ZK)
                if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then
                    if (i .lt. MODLAYS)  then
                        FLAI(i) = ((FCLAI(i+1) - FCLAI(i)) * LAI)/MODRES    !fractional LAI in each layer converted to LAD (m2 m-3)
                    else
                        FLAI(i) = FLAI(MODLAYS-1)
                    end if
                end if
            end do
            do i=1,  SIZE(ZK)
                if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then
                    VPGWT(i) = (FLAI(i))/sum(FLAI(1:LAYERS))
                else
                    VPGWT(i) = 0.0_rk !above canopy
                end if
            end do
            EMI_OUT(SIZE(ZK)) = LAI * EF * SUM(GammaTLEAF_AVE(1:LAYERS) * GammaPPFD_AVE(1:LAYERS) * &
                VPGWT(1:LAYERS)) * GAMMACO2 * CCE * GAMMALEAFAGE * GAMMASOIM * CANLOSS_FAC !put into top model layer (ug m-2 hr-1)
            EMI_OUT = EMI_OUT * 2.7777777777778E-13_rk    !convert emissions output to    (kg m-2 s-1)
        else if (VERT .eq. 2) then       !"MEGANv3-like": Add weighted sum of activity coefficients using normal distribution
            !across canopy layers using 5 layer numbers directly from MEGANv3
            !--warning: weights are not consistent with FCLAI distribution
            !used for biomass distribution used for sunlit/shaded in Gamma TLEAF and GammaPPFD.
            LAYERS = min((floor(FCH/MODRES) + 1),MODLAYS)
            do i=1, SIZE(ZK)
                if (ZK(i) .gt. FCH) then
                    GAUSS(i) = 0.0
                else if (ZK(i) .le. FCH .and. ZK(i) .gt. FCH*(4.0_rk/5.0_rk)) then  !Level 1 - 2
                    GAUSS(i)   = interp_linear1_internal((/ FCH*(4.0_rk/5.0_rk),FCH /), &
                        (/ 0.118464_rk,0.0_rk /),ZK(i))
                else if (ZK(i) .le. FCH*(4.0_rk/5.0_rk) .and. ZK(i) .gt. FCH*(3.0_rk/5.0_rk)) then  !Level 2 - 3
                    GAUSS(i)   = interp_linear1_internal((/ FCH*(3.0_rk/5.0_rk),FCH*(4.0_rk/5.0_rk) /), &
                        (/ 0.239314_rk,0.118464_rk /),ZK(i))
                else if (ZK(i) .le. FCH*(3.0_rk/5.0_rk) .and. ZK(i) .gt. FCH*(2.0_rk/5.0_rk)) then  !Level 3 - 4
                    GAUSS(i)   = interp_linear1_internal((/ FCH*(2.0_rk/5.0_rk),FCH*(3.0_rk/5.0_rk) /), &
                        (/ 0.284444_rk,0.239314_rk /),ZK(i))
                else if (ZK(i) .le. FCH*(2.0_rk/5.0_rk) .and. ZK(i) .gt. FCH*(1.0_rk/5.0_rk) ) then  !Level 4 - Bottom
                    GAUSS(i)   = interp_linear1_internal((/ FCH*(1.0_rk/5.0_rk),FCH*(2.0_rk/5.0_rk) /), &
                        (/ 0.239314_rk,0.284444_rk /),ZK(i))
                else if (ZK(i) .le. FCH*(1.0_rk/5.0_rk) ) then  !Level 4 - Bottom
                    GAUSS(i)   = interp_linear1_internal((/ ZK(1),FCH*(1.0_rk/5.0_rk) /), &
                        (/ 0.118464_rk,0.239314_rk /),ZK(i))
                end if
            end do

            do i=1, SIZE(ZK)
                VPGWT(i) = GAUSS(i)/sum(GAUSS(1:LAYERS))
            end do
            EMI_OUT(SIZE(ZK)) = LAI * EF * SUM(GammaTLEAF_AVE(1:LAYERS) * GammaPPFD_AVE(1:LAYERS) * &
                VPGWT(1:LAYERS)) * GAMMACO2 * CCE * GAMMALEAFAGE * GAMMASOIM * CANLOSS_FAC    !put into top model layer (ug m-2 hr-1)
            EMI_OUT = EMI_OUT * 2.7777777777778E-13_rk    !convert emissions output to    (kg m-2 s-1)
        else if (VERT .eq. 3) then       !"MEGANv3-like": Add weighted sum of activity coefficients equally
            !across canopy layers
            !--warning: weights are not consistent with FCLAI distribution
            !used for biomass distribution used for sunlit/shaded in Gamma TLEAF and GammaPPFD.
            LAYERS = min((floor(FCH/MODRES) + 1),MODLAYS)
            do i=1,  SIZE(ZK)
                VPGWT(i) = 1.0_rk/LAYERS
            end do
            EMI_OUT(SIZE(ZK)) = LAI * EF * SUM(GammaTLEAF_AVE(1:LAYERS) * GammaPPFD_AVE(1:LAYERS) * &
                VPGWT(1:LAYERS)) * GAMMACO2 * CCE * GAMMALEAFAGE * GAMMASOIM * CANLOSS_FAC    !put into top model layer (ug m-2 hr-1)
            EMI_OUT = EMI_OUT * 2.7777777777778E-13_rk    !convert emissions output to    (kg m-2 s-1)
        else
            write(*,*)  'Wrong BIOVERT_OPT choice of ', VERT, ' in namelist...exiting'
            call exit(2)
        end if

    END SUBROUTINE CANOPY_BIO

end module canopy_bioemi_mod
