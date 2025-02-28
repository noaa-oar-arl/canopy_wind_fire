module canopy_drydep_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_GAS_DRYDEP_ZHANG( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, &
        ZK, FCH, TEMPA, PRESSA, &
        RELHUMA, FSUN, PPFD_SUN, PPFD_SHADE, UBAR, &
        SRAD, DEP_IND, DEP_OUT)

!-----------------------------------------------------------------------

! Description:
!     computes parameterized canopy gas dry deposition based on Zhang et al. (2003)

! Revision History:
!     Prototype 02/25 by PCC, based on Zhang et al. (2003) algorithms as implemented
!                             in ACCESS (Saylor 2013)
! Citation:
! Zhang, L., Brook, J. R., and Vet, R.: A revised parameterization for gaseous dry
! deposition in air-quality models, Atmos. Chem. Phys., 3, 2067–2082,
! https://doi.org/10.5194/acp-3-2067-2003, 2003.
!
! Saylor, R. D.: The Atmospheric Chemistry and Canopy Exchange Simulation System (ACCESS):
! model description and application to a temperate deciduous forest canopy,
! Atmos. Chem. Phys., 13, 693–715, https://doi.org/10.5194/acp-13-693-2013, 2013.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Feb 2025 P.C. Campbell: Initial Zhang et al. gas dry deposition
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod,  ONLY: rk                       !constants for canopy models
        use canopy_utils_mod,  ONLY: MolecDiff,rs_zhang_gas,EffHenrysLawCoeff,&
            ReactivityParam,rbl,rcl,rml

! Arguments:
!     IN/OUT
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_OPT    ! Select chemical mechanism
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_TOT    ! Select chemical mechanism gas species list
        REAL(RK),    INTENT( IN )       :: ZK(:)              ! Model heights (m)
        REAL(RK),    INTENT( IN )       :: FCH                ! Canopy height (m)
        REAL(RK),    INTENT( IN )       :: FSUN(:)            ! Sunlit/Shaded fraction from photolysis correction factor
        REAL(RK),    INTENT( IN )       :: PPFD_SUN(:)        ! PPFD for sunlit leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: PPFD_SHADE(:)      ! PPFD for shaded leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: TEMPA(:)           ! Ambient Temperature profile in canopy [K]
        REAL(RK),    INTENT( IN )       :: PRESSA(:)          ! Ambient Pressure profile in canopy [mb]
        REAL(RK),    INTENT( IN )       :: RELHUMA(:)         ! Ambient Relative Humidity profile in canopy [%]
        REAL(RK),    INTENT( IN )       :: UBAR(:)            ! Mean above/in-canopy wind speed [m/s]
        REAL(RK),    INTENT( IN )       :: SRAD               ! Incoming solar irradiation top of canopy (W/m^2)
        INTEGER,     INTENT( IN )       :: DEP_IND            ! Gas deposition species index (depends on gas mech, set in constants)
        REAL(RK),    INTENT( OUT )      :: DEP_OUT(:)         ! Output canopy layer gas dry deposition rate for each DEP_IND

! Local Variables
        REAL(RK) ::  PPFD(SIZE(ZK))                           ! PPFD ave sun and shade (umol/m2 s)
        REAL(RK) ::  mdiffl(SIZE(ZK))                         ! Molecular diffusivity for species l based on DEP_IND [cm2/s]
        REAL(RK) ::  rs(SIZE(ZK))                             ! Stomatal resistance for species l based on DEP_IND [s/cm]
        REAL(RK) ::  rb(SIZE(ZK))                             ! Boundary layer resistance for species l based on DEP_IND [s/cm]
        REAL(RK) ::  rc(SIZE(ZK))                             ! Cuticular resistance for species l based on DEP_IND [s/cm]
        REAL(RK) ::  rm(SIZE(ZK))                             ! Mesophyll resistance for species l based on DEP_IND [s/cm]
        REAL(RK) ::  hstarl                                   ! effective Henry's law coefficient based on DEP_IND (M/atm)
        REAL(RK) ::  f01                                      ! reactivity parameter based on DEP_IND (0-1)

        REAL(RK) :: rnum,rden,rlx,vdlx
        INTEGER i

        PPFD = (PPFD_SUN*FSUN) + (PPFD_SHADE*(1.0-FSUN)) ! average = sum sun and shade weighted by sunlit fraction

!! Calculate molecular diffusivity (cm^2/s) and resistances (cm/s) of species l using from DEP_IND
        hstarl  = EffHenrysLawCoeff(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)
        f01     = ReactivityParam(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)
        do i=1, SIZE(ZK)
            if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then           ! above ground level and at/below canopy top
                mdiffl(i)  = MolecDiff(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND,TEMPA(i),PRESSA(i))
                rs(i)      = rs_zhang_gas(mdiffl(i),TEMPA(i),PRESSA(i),PPFD(i),SRAD,RELHUMA(i)) !stomatal resistance (s/cm)
                rb(i)      = rbl(mdiffl(i), UBAR(i)*100.0_rk) !leaf boundary layer resistance (s/cm)
                rc(i)      = rcl(hstarl, f01)                              !leaf cuticular resistance (s/cm)
                rm(i)      = rml(hstarl, f01)                              !leaf mesophyll resistance (s/cm)
                rnum = rc(i) * (rs(i) + rm(i))
                rden = rc(i) + 2.0_rk * (rs(i) + rm(i))
                rlx   = rb(i) + (rnum/rden)
                vdlx  = 1.0_rk/rlx
                dep_out(i) = vdlx                                          !calculate deposition velocity (cm/s)
            else
                rb(i) = 0.0_rk
                rc(i) = 0.0_rk
                rm(i) = 0.0_rk
                rs(i) = 0.0_rk
                dep_out(i) = 0.0_rk
            endif
        end do

    END SUBROUTINE CANOPY_GAS_DRYDEP_ZHANG


    SUBROUTINE CANOPY_GAS_DRYDEP_SOIL( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, &
        TEMPA, PRESSA, UBAR, SOCAT, SOTYP, DSOIL, STHETA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk                       !constants for canopy models
        use canopy_utils_mod,  ONLY: MolecDiff,SoilResist,SoilRbg

        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_OPT    ! Select chemical mechanism
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_TOT    ! Select chemical mechanism gas species list
        REAL(RK),    INTENT( IN )       :: TEMPA(:)           ! Ambient Temperature profile in canopy [K]
        REAL(RK),    INTENT( IN )       :: PRESSA(:)          ! Ambient Pressure profile in canopy [mb]
        REAL(RK),    INTENT( IN )       :: UBAR(:)            ! Mean above/in-canopy wind speed [m/s]
        INTEGER,     INTENT( IN )       :: SOCAT              ! input soil category datset used
        INTEGER,     INTENT( IN )       :: SOTYP              ! input soil type integer associated with soilcat
        REAL(RK),    INTENT( IN )       :: DSOIL              ! depth of topsoil (cm)
        REAL(RK),    INTENT( IN )       :: STHETA             ! volumetric soil water content in topsoil(m^3/m^3)

        INTEGER,     INTENT( IN )       :: DEP_IND            ! Gas deposition species index (depends on gas mech)
        REAL(RK),    INTENT( OUT )      :: DEP_OUT            ! Output soil layer gas dry deposition rate for each DEP_IND

        real(rk)                        :: mdiffl             ! molecular diffusivity (cm^2/s)
        real(rk)                        :: rsoill             ! resistance to diffusion thru soil pore space for chemical species (s/cm)
        real(rk)                        :: rbg                ! ground boundary layer resistance (s/cm)

        mdiffl = MolecDiff(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND,TEMPA(1),PRESSA(1))  !Use surface temperature and pressure

        rsoill = SoilResist(mdiffl,SOCAT,SOTYP,DSOIL,STHETA)

        rbg = SoilRbg(UBAR(2)*100.0_rk) !convert wind to cm/s   !Rbg(ground boundary layer resistance, s/cm)
        !Rbg is invariant to species not layers
        !Must use second model layer as no-slip boundary condition for wind
        !speed at z = 0

        DEP_OUT = 1.0/(rbg+rsoill)                               !deposition velocity to ground surface under canopy (cm/s)

        return
    END SUBROUTINE CANOPY_GAS_DRYDEP_SOIL

end module canopy_drydep_mod
