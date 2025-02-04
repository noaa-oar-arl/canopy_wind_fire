module canopy_drydep_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_GAS_DRYDEP_ZHANG( TEMPA, PRESSA, &
        RELHUMA, FSUN, PPFD_SUN, PPFD_SHADE, SRAD, DEP_IND, DEP_OUT)

!-----------------------------------------------------------------------

! Description:
!     computes parameterized canopy gas dry deposition based on Zhang et al. (2003)

! Preconditions:
!     in-canopy FCLAI, model LAI, etc.

! Subroutines and Functions Called:

! Revision History:
!     Prototype 02/25 by PCC, based on Zhang et al. (2003) algorithms
! Citation:
! Zhang, L., Brook, J. R., and Vet, R.: A revised parameterization for gaseous dry
! deposition in air-quality models, Atmos. Chem. Phys., 3, 2067–2082,
! https://doi.org/10.5194/acp-3-2067-2003, 2003.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Feb 2025 P.C. Campbell: Initial Zhang et al. gas dry deposition
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod,  ONLY: rk                       !constants for canopy models
        use canopy_utils_mod,  ONLY: MolecDiff,rs_zhang_gas

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )       :: FSUN(:)            ! Sunlit/Shaded fraction from photolysis correction factor
        REAL(RK),    INTENT( IN )       :: PPFD_SUN(:)        ! PPFD for sunlit leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: PPFD_SHADE(:)      ! PPFD for shaded leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: TEMPA(:)           ! Ambient Temperature profile in canopy [K]
        REAL(RK),    INTENT( IN )       :: PRESSA(:)          ! Ambient Pressure profile in canopy [mb]
        REAL(RK),    INTENT( IN )       :: RELHUMA(:)         ! Ambient Relative Humidity profile in canopy [%]
        REAL(RK),    INTENT( IN )       :: SRAD               ! Incoming solar irradiation top of canopy (W/m^2)
        INTEGER,     INTENT( IN )       :: DEP_IND            ! Gas deposition species index (depends on gas mech, set in constants)
        REAL(RK),    INTENT( OUT )      :: DEP_OUT(:)         ! Output canopy layer gas dry deposition rate for each DEP_IND

! Local Variables
        REAL(RK) ::  PPFD(SIZE(FSUN))                         ! PPFD ave sun and shade (umol/m2 s)
        REAL(RK) ::  mdiffl(SIZE(FSUN))                       ! Molecular diffusivity for species l based on DEP_IND [cm2/s]
        REAL(RK) ::  rs(SIZE(FSUN))                           ! Stomatal resistance for species l based on DEP_IND [s/cm]

!        REAL(RK) :: GammaTLEAF_SUN_NUM(SIZE(ZK))   ! Numerator in Tleaf sun activity factor
!        REAL(RK) :: GammaTLEAF_SHADE_NUM(SIZE(ZK)) ! Numerator in Tleaf shade activity factor
!        REAL(RK) :: GammaTLEAF_SUN_DEN(SIZE(ZK))   ! Denominator in Tleaf sun activity factor
!        REAL(RK) :: GammaTLEAF_SHADE_DEN(SIZE(ZK)) ! Denominator in Tleaf sun activity factor
!        REAL(RK) :: GammaTLEAF_SUN(SIZE(ZK))       ! Tleaf sun activity factor
!        REAL(RK) :: GammaTLEAF_SHADE(SIZE(ZK))     ! Tleaf shade activity factor
!        REAL(RK) :: GammaTLEAF_AVE(SIZE(ZK))       ! Average Tleaf activity factor
!        REAL(RK) :: CP_SUN(SIZE(ZK))               ! Normalized emission capacity sun at PPFD = 1000 umol phot/m2 s
!        REAL(RK) :: CP_SHADE(SIZE(ZK))             ! Normalized emission capacity shade at PPFD = 1000 umol phot/m2 s
!        REAL(RK) :: ALPHA_P_SUN(SIZE(ZK))          ! Quantum yield of isoprene sunlit (mol/mol)
!        REAL(RK) :: ALPHA_P_SHADE(SIZE(ZK))        ! Quantum yield of isoprene shade (mol/mol)
!        REAL(RK) :: GammaPPFD_SUN(SIZE(ZK))        ! PPFD activity factor sun (unitless)
!        REAL(RK) :: GammaPPFD_SHADE(SIZE(ZK))      ! PPFD activity factor shade (unitless)
!        REAL(RK) :: GammaPPFD_AVE(SIZE(ZK))        ! PPFD activity factor ave sun and shade
!        REAL(RK) :: E_OPT(SIZE(ZK))                ! maximum normalized emission capacity
!        REAL(RK) :: TLEAF_OPT(SIZE(ZK))            ! Tleaf at which E_OPT occurs (K)
!        REAL(RK) :: FLAI(SIZE(ZK))                 ! Fractional LAI in layer
!        REAL(RK) :: VPGWT(SIZE(ZK))                ! MEGANv3-like in-canopy weighting factor
!        REAL(RK) :: GAUSS(SIZE(ZK))                ! MEGANv3-like in-canopy gaussian
!        REAL(RK) :: CT1                            ! Activation energy (kJ/mol)
!        REAL(RK) :: CEO                            ! Empirical coefficient
!        REAL(RK) :: EF                             ! Final Mapped Emission factor (EF) (ug/m2 hr)
!
!
!        REAL(RK) :: TABOVECANOPY !(SIZE(ZK))  ! Above Canopy Temp assigned = TEMP2 i.e., Model input 2-m Temperature (K for now)
!        ! Empirical coeff.'s for Leaf Age factor calculations (see
!        ! canopy_bioparm_mod or call canopy_biop)
!        REAL(RK) :: ANEW
!        REAL(RK) :: AGRO
!        REAL(RK) :: AMAT
!        REAL(RK) :: AOLD
!
!        !Coeff.'s and threshold/delta threshold for air quality stress factors from canopy_biop
!        REAL(RK) :: CAQ
!        REAL(RK) :: TAQ  ![ppm-hours]
!        REAL(RK) :: DTAQ ![ppm-hours]
!
!        !Coeff.'s and threshold/delta threshold for high temperature stress factors from canopy_biop
!        REAL(RK) :: CHT
!        REAL(RK) :: THT  ![K]
!        REAL(RK) :: DTHT ![K]
!
!        !Coeff.'s and threshold/delta threshold for low temperature stress factors from canopy_biop
!        REAL(RK) :: CLT
!        REAL(RK) :: TLT  ![K]
!        REAL(RK) :: DTLT ![K]
!
!        !Coeff.'s and threshold/delta threshold for high wind stress factors from canopy_biop
!        REAL(RK) :: CHW
!        REAL(RK) :: THW  ![m/s]
!        REAL(RK) :: DTHW ![m/s]
!
!
!        ! Coefficients A and B used for PFT dependent cumulative root depth fraction
!        REAL(RK) :: ROOTA ! [m-1]
!        REAL(RK) :: ROOTB ! [m-1]
!        REAL(RK) :: GAMMASOIM ! Soil moisture factor
!
!        REAL(RK) :: GAMMAAQ                        !Air quality stress factor
!        REAL(RK) :: GAMMAHT                        !High temperature stress factor
!        REAL(RK) :: GAMMALT                        !Low temperature stress factor
!        REAL(RK) :: GAMMAHW                        !High wind speed stress factor
!
!        REAL(RK) :: GAMMACO2                       ! CO2 inhibition factor (isoprene only)
!
!        REAL(RK) :: GAMMALEAFAGE !(SIZE(ZK))                 ! LEAF AGE factor
!
!        REAL(RK) :: CANLOSS_FAC                    !Canopy loss factor for summing option
!
        integer i
!
!! Constant Canopy Parameters
!        REAL(RK),          PARAMETER     :: PPFD0_SUN       =  200.0      !Constant PPFDo sunlit (umol/m2 s) (Guenther et al.,2012)
!        REAL(RK),          PARAMETER     :: PPFD0_SHADE     =  50.0       !Constant PPFDo shaded (umol/m2 s) (Guenther et al.,2012)
!        REAL(RK),          PARAMETER     :: CT2             =  230.0_rk   !Deactivation energy (kJ/mol) (Guenther et al., 2012)


        PPFD = (PPFD_SUN*FSUN) + (PPFD_SHADE*(1.0-FSUN)) ! average = sum sun and shade weighted by sunlit fraction

!! Calculate molecular diffusivity (cm^2/s) and stomatal resistance (cm/s) of species l using from DEP_IND
        do i=1, SIZE(TEMPA)
            mdiffl(i)  = MolecDiff(DEP_IND,TEMPA(i),PRESSA(i))
            rs(i)    = rs_zhang_gas(mdiffl(i),TEMPA(i),PRESSA(i),PPFD(i),SRAD,RELHUMA(i))
        end do

        dep_out = rs

        !
!! Calculate maximum normalized emission capacity (E_OPT) and Tleaf at E_OPT
!        TLEAF_OPT = 313.0_rk + (0.6_rk * (TLEAF240_AVE-297.0_rk)) !Guenther et al. (2012)
!
!! Calculate emission species/plant-dependent mapped emission factors and other important coefficients for gamma terms
!        call canopy_biop(EMI_IND, LU_OPT, VTYPE, EF, CT1, CEO, ANEW, AGRO, AMAT, AOLD, ROOTA, ROOTB, CAQ, TAQ, DTAQ, &
!            CHT, THT, DTHT, CLT, TLT, DTLT, CHW, THW, DTHW)
!
!!        print*,'CHT=',CHT,'THT=',THT,'DTHT=',DTHT
!!        print*,'CLT=',CLT,'TLT=',TLT,'DTLT=',DTLT
!!        print*,'CHW=',CHW,'THW=',THW,'DTHW=',DTHW
!        E_OPT = CEO * EXP(0.05_rk * (TLEAF24_AVE-297.0_rk)) * EXP(0.05_rk * (TLEAF240_AVE-297.0_rk))
!
!! Calculate gamma (activity) values for average Tleaf (Clifton et al., 2022; based on Guenther et al. 2012)
!        GammaTLEAF_SUN_NUM = CT2*exp((CT1/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SUN)))
!        GammaTLEAF_SUN_DEN = (CT2-CT1)*(1.0-exp((CT2/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SUN))))
!        GammaTLEAF_SUN     = E_OPT*(GammaTLEAF_SUN_NUM/GammaTLEAF_SUN_DEN)
!
!        GammaTLEAF_SHADE_NUM = CT2*exp((CT1/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SHADE)))
!        GammaTLEAF_SHADE_DEN = (CT2-CT1)*(1.0-exp((CT2/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SHADE))))
!        GammaTLEAF_SHADE     = E_OPT*(GammaTLEAF_SHADE_NUM/GammaTLEAF_SHADE_DEN)
!
!        GammaTLEAF_AVE = (GammaTLEAF_SUN*FSUN) + (GammaTLEAF_SHADE*(1.0-FSUN)) ! average = sum sun and shade weighted by sunlit fraction
!
!        GammaTLEAF_AVE = MAX( GammaTLEAF_AVE, 0.0_rk )
!
!! Calculate gamma (activity) values for average PPFD (Clifton et al., 2022; based on Guenther et al. 2012)
!        ALPHA_P_SUN = 0.004 - 0.0005*log(PPFD240_SUN)
!        ALPHA_P_SHADE = 0.004 - 0.0005*log(PPFD240_SHADE)
!        CP_SUN = 0.0468*(PPFD240_SUN**(0.6))*exp(0.005*(PPFD24_SUN-PPFD0_SUN))
!        CP_SHADE = 0.0468*(PPFD240_SHADE**(0.6))*exp(0.005*(PPFD24_SHADE-PPFD0_SHADE))
!        GammaPPFD_SUN   = CP_SUN*((ALPHA_P_SUN*PPFD_SUN)/SQRT(1.0 + (ALPHA_P_SUN**2.0) * (PPFD_SUN**2.0)))
!        GammaPPFD_SHADE = CP_SHADE*((ALPHA_P_SHADE*PPFD_SHADE)/SQRT(1.0 + (ALPHA_P_SHADE**2.0) * (PPFD_SHADE**2.0)))
!
!        GammaPPFD_AVE = (GammaPPFD_SUN*FSUN) + (GammaPPFD_SHADE*(1.0-FSUN)) ! average = sum sun and shade weighted by sunlit fraction
!
!        GammaPPFD_AVE = MAX( GammaPPFD_AVE, 0.0_rk )
!
!! Get CO2 inhibition factor for isoprene only
!
!        if (EMI_IND .eq. 1) then  !Isoprene
!            GAMMACO2 = GET_GAMMA_CO2(CO2OPT,CO2SET)
!        else
!            GAMMACO2 = 1.0_rk
!        end if
!
!! Get Soil Moisture Factor
!        GAMMASOIM =  GET_GAMMA_SOIM(SOIMOPT,SOIM1,SOIM2,SOIM3,SOIM4,SOID1,SOID2,SOID3,SOID4,WILT, &
!            ROOTA,ROOTB)
!
!! Get LEAF AGE factor
!
!        TABOVECANOPY  = TEMP2   !TEMP2 (above air temp) for TABOVECANOPY
!        !do i=1, SIZE(ZK)
!        GAMMALEAFAGE = GET_GAMMA_LEAFAGE(LEAFAGEOPT, PASTLAI, CURRENTLAI, TSTEPLAI, TABOVECANOPY, ANEW, AGRO, AMAT, AOLD)
!        !end do
!
!! Get AQ Stress Factor
!        GAMMAAQ = GET_GAMMA_AQ(AQOPT, W126_REF, W126_SET, CAQ, TAQ, DTAQ)
!
!! Get HT Stress Factor
!        GAMMAHT = GET_GAMMA_HT(HTOPT, DAILY_MAXT2, CHT, THT, DTHT)
!
!! Get LT Stress Factor
!        GAMMALT = GET_GAMMA_LT(LTOPT, DAILY_MINT2, CLT, TLT, DTLT)
!
!! Get HW Stress Factor
!        GAMMAHW = GET_GAMMA_HW(HWOPT, DAILY_MAXWS10, CHW, THW, DTHW)
!
!! Get canopy loss factor (only used in vertical summing options and empirical formulation and parameters based on isoprene)
!! Note:  Allowed for other BVOCs but use caution when applying to compare with above canopy flux observations
!
!        CANLOSS_FAC = 1.0_rk  !Initialize
!        ! All species
!        if (LOSSIND .eq. 0) then
!            if (LOSSOPT .eq. 2) then !User set value from NL
!                CANLOSS_FAC = LOSSSET
!            else                !Try and calculate if turned on
!                CANLOSS_FAC = GET_CANLOSS_BIO(LOSSOPT,LIFETIME,USTAR,FCH)
!            end if
!        end if
!
!        !Only for a specific biogenic species/indice
!        if (LOSSIND .eq. EMI_IND) then
!            if (LOSSOPT .eq. 2) then !User set value from NL
!                CANLOSS_FAC = LOSSSET
!            else                !Try and calculate if turned on
!                CANLOSS_FAC = GET_CANLOSS_BIO(LOSSOPT,LIFETIME,USTAR,FCH)
!            end if
!        end if
!
!! Calculate emissions profile in the canopy
!        EMI_OUT = 0.0_rk  ! set initial emissions profile to zero
!        FLAI = 0.0_rk  ! set initial fractional FLAI (LAD) profile to zero
!
!        if (VERT .eq. 0) then         !Full 3D leaf-level biogenic emissions (no averaging, summing, or integration)
!            do i=1, SIZE(ZK)
!                if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then           ! above ground level and at/below canopy top
!                    if (i .lt. MODLAYS)  then
!                        FLAI(i) = ((FCLAI(i+1) - FCLAI(i)) * LAI)/MODRES    !fractional LAI in each layer converted to LAD (m2 m-3)
!                    else
!                        FLAI(i) = FLAI(MODLAYS-1)
!                    end if
!                    EMI_OUT(i) = FLAI(i) * EF * GammaTLEAF_AVE(i) * GammaPPFD_AVE(i) * GAMMACO2 * CCE * &
!                        GAMMALEAFAGE * GAMMASOIM * GAMMAAQ * GAMMAHT * GAMMALT * GAMMAHW ! (ug m-3 hr-1)
!                    EMI_OUT(i) = EMI_OUT(i) * 2.7777777777778E-13_rk    !convert emissions output to (kg m-3 s-1)
!                end if
!            end do
!        else if (VERT .eq. 1) then       !"MEGANv3-like": Use weighting factors normalized to plant distribution shape (FCLAI)
!            !across canopy layers
!            LAYERS = min((floor(FCH/MODRES) + 1),MODLAYS)
!            do i=1,  SIZE(ZK)
!                if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then
!                    if (i .lt. MODLAYS)  then
!                        FLAI(i) = ((FCLAI(i+1) - FCLAI(i)) * LAI)/MODRES    !fractional LAI in each layer converted to LAD (m2 m-3)
!                    else
!                        FLAI(i) = FLAI(MODLAYS-1)
!                    end if
!                end if
!            end do
!            do i=1,  SIZE(ZK)
!                if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then
!                    VPGWT(i) = (FLAI(i))/sum(FLAI(1:LAYERS))
!                else
!                    VPGWT(i) = 0.0_rk !above canopy
!                end if
!            end do
!            EMI_OUT(SIZE(ZK)) = LAI * EF * SUM(GammaTLEAF_AVE(1:LAYERS) * GammaPPFD_AVE(1:LAYERS) * &
!                VPGWT(1:LAYERS)) * GAMMACO2 * CCE * GAMMALEAFAGE * GAMMASOIM * &
!                GAMMAAQ * GAMMAHT * GAMMALT * GAMMAHW * CANLOSS_FAC !put into top model layer (ug m-2 hr-1)
!            EMI_OUT = EMI_OUT * 2.7777777777778E-13_rk    !convert emissions output to    (kg m-2 s-1)
!        else if (VERT .eq. 2) then       !"MEGANv3-like": Add weighted sum of activity coefficients using normal distribution
!            !across canopy layers using 5 layer numbers directly from MEGANv3
!            !--warning: weights are not consistent with FCLAI distribution
!            !used for biomass distribution used for sunlit/shaded in Gamma TLEAF and GammaPPFD.
!            LAYERS = min((floor(FCH/MODRES) + 1),MODLAYS)
!            do i=1, SIZE(ZK)
!                if (ZK(i) .gt. FCH) then
!                    GAUSS(i) = 0.0
!                else if (ZK(i) .le. FCH .and. ZK(i) .gt. FCH*(4.0_rk/5.0_rk)) then  !Level 1 - 2
!                    GAUSS(i)   = interp_linear1_internal((/ FCH*(4.0_rk/5.0_rk),FCH /), &
!                        (/ 0.118464_rk,0.0_rk /),ZK(i))
!                else if (ZK(i) .le. FCH*(4.0_rk/5.0_rk) .and. ZK(i) .gt. FCH*(3.0_rk/5.0_rk)) then  !Level 2 - 3
!                    GAUSS(i)   = interp_linear1_internal((/ FCH*(3.0_rk/5.0_rk),FCH*(4.0_rk/5.0_rk) /), &
!                        (/ 0.239314_rk,0.118464_rk /),ZK(i))
!                else if (ZK(i) .le. FCH*(3.0_rk/5.0_rk) .and. ZK(i) .gt. FCH*(2.0_rk/5.0_rk)) then  !Level 3 - 4
!                    GAUSS(i)   = interp_linear1_internal((/ FCH*(2.0_rk/5.0_rk),FCH*(3.0_rk/5.0_rk) /), &
!                        (/ 0.284444_rk,0.239314_rk /),ZK(i))
!                else if (ZK(i) .le. FCH*(2.0_rk/5.0_rk) .and. ZK(i) .gt. FCH*(1.0_rk/5.0_rk) ) then  !Level 4 - Bottom
!                    GAUSS(i)   = interp_linear1_internal((/ FCH*(1.0_rk/5.0_rk),FCH*(2.0_rk/5.0_rk) /), &
!                        (/ 0.239314_rk,0.284444_rk /),ZK(i))
!                else if (ZK(i) .le. FCH*(1.0_rk/5.0_rk) ) then  !Level 4 - Bottom
!                    GAUSS(i)   = interp_linear1_internal((/ ZK(1),FCH*(1.0_rk/5.0_rk) /), &
!                        (/ 0.118464_rk,0.239314_rk /),ZK(i))
!                end if
!            end do
!
!            do i=1, SIZE(ZK)
!                VPGWT(i) = GAUSS(i)/sum(GAUSS(1:LAYERS))
!            end do
!            EMI_OUT(SIZE(ZK)) = LAI * EF * SUM(GammaTLEAF_AVE(1:LAYERS) * GammaPPFD_AVE(1:LAYERS) * &
!                VPGWT(1:LAYERS)) * GAMMACO2 * CCE * GAMMALEAFAGE * GAMMASOIM * &
!                GAMMAAQ * GAMMAHT * GAMMALT * GAMMAHW * CANLOSS_FAC    !put into top model layer (ug m-2 hr-1)
!            EMI_OUT = EMI_OUT * 2.7777777777778E-13_rk    !convert emissions output to    (kg m-2 s-1)
!        else if (VERT .eq. 3) then       !"MEGANv3-like": Add weighted sum of activity coefficients equally
!            !across canopy layers
!            !--warning: weights are not consistent with FCLAI distribution
!            !used for biomass distribution used for sunlit/shaded in Gamma TLEAF and GammaPPFD.
!            LAYERS = min((floor(FCH/MODRES) + 1),MODLAYS)
!            do i=1,  SIZE(ZK)
!                VPGWT(i) = 1.0_rk/LAYERS
!            end do
!            EMI_OUT(SIZE(ZK)) = LAI * EF * SUM(GammaTLEAF_AVE(1:LAYERS) * GammaPPFD_AVE(1:LAYERS) * &
!                VPGWT(1:LAYERS)) * GAMMACO2 * CCE * GAMMALEAFAGE * GAMMASOIM * &
!                GAMMAAQ * GAMMAHT * GAMMALT * GAMMAHW  * CANLOSS_FAC    !put into top model layer (ug m-2 hr-1)
!            EMI_OUT = EMI_OUT * 2.7777777777778E-13_rk    !convert emissions output to    (kg m-2 s-1)
!        else
!            write(*,*)  'Wrong BIOVERT_OPT choice of ', VERT, ' in namelist...exiting'
!            call exit(2)
!        end if

    END SUBROUTINE CANOPY_GAS_DRYDEP_ZHANG

end module canopy_drydep_mod
