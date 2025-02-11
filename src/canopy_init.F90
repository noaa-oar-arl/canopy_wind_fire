
SUBROUTINE canopy_init

!-------------------------------------------------------------------------------
! Name:     Initialize Arrays for Canopy Inputs and Outputs
! Purpose:  Initialize arrays for Canopy Inputs and Outputs
! Revised:  03 Oct 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------

    USE canopy_const_mod, ONLY: rk, fillreal
    USE canopy_canopts_mod
    USE canopy_coord_mod
    USE canopy_canmet_mod
    USE canopy_canvars_mod

    IMPLICIT NONE

!-------------------------------------------------------------------------------
! Initialize arrays for Canopy Distribution
!-------------------------------------------------------------------------------

    if(allocated(zhc))                zhc(:)                 = fillreal
    if(allocated(fafraczInt))         fafraczInt(:)          = fillreal
    if(allocated(fsun))               fsun(:)                = fillreal
    if(allocated(tleaf_sun))          tleaf_sun(:)           = fillreal
    if(allocated(tleaf_shade))        tleaf_shade(:)         = fillreal
    if(allocated(tleaf_ave))          tleaf_ave(:)           = fillreal
    if(allocated(ppfd_sun))           ppfd_sun(:)            = fillreal
    if(allocated(ppfd_shade))         ppfd_shade(:)          = fillreal
    if(allocated(ppfd_ave))           ppfd_ave(:)            = fillreal
    if(allocated(lad))                lad(:,:)               = fillreal
    if(allocated(lad_3d))             lad_3d(:,:,:)          = fillreal
    if(allocated(tka))                tka(:,:)               = fillreal
    if(allocated(tka_3d))             tka_3d(:,:,:)          = fillreal
    if(allocated(pressa))             pressa(:,:)            = fillreal
    if(allocated(pressa_3d))          pressa_3d(:,:,:)       = fillreal
    if(allocated(relhuma))            relhuma(:,:)           = fillreal
    if(allocated(relhuma_3d))         relhuma_3d(:,:,:)      = fillreal
    if(allocated(spechuma))           spechuma(:,:)          = fillreal
    if(allocated(spechuma_3d))        spechuma_3d(:,:,:)     = fillreal

    if(allocated(tleaf_sun24_tmp))    tleaf_sun24_tmp(:,:,:)     = 0.0_rk
    if(allocated(tleaf_shade24_tmp))  tleaf_shade24_tmp(:,:,:)   =  0.0_rk
    if(allocated(tleaf_ave24_tmp))    tleaf_ave24_tmp(:,:,:)     = 0.0_rk
    if(allocated(ppfd_sun24_tmp))     ppfd_sun24_tmp(:,:,:)      = 0.0_rk
    if(allocated(ppfd_shade24_tmp))   ppfd_shade24_tmp(:,:,:)    = 0.0_rk
    if(allocated(tleaf_sun240_tmp))   tleaf_sun240_tmp(:,:,:)    = 0.0_rk
    if(allocated(tleaf_shade240_tmp)) tleaf_shade240_tmp(:,:,:)  = 0.0_rk
    if(allocated(tleaf_ave240_tmp))   tleaf_ave240_tmp(:,:,:)    = 0.0_rk
    if(allocated(ppfd_sun240_tmp))    ppfd_sun240_tmp(:,:,:)     = 0.0_rk
    if(allocated(ppfd_shade240_tmp))  ppfd_shade240_tmp(:,:,:)   = 0.0_rk
    if(allocated(tmp2mref_tmp))       tmp2mref_tmp(:,:)          = 0.0_rk
    if(allocated(ubzref_tmp))         ubzref_tmp(:,:)            = 0.0_rk
    if(allocated(tleaf_sun24))        tleaf_sun24(:,:)         = fillreal
    if(allocated(tleaf_shade24))      tleaf_shade24(:,:)       = fillreal
    if(allocated(tleaf_ave24))        tleaf_ave24(:,:)         = fillreal
    if(allocated(ppfd_sun24))         ppfd_sun24(:,:)          = fillreal
    if(allocated(ppfd_shade24))       ppfd_shade24(:,:)        = fillreal
    if(allocated(tleaf_sun240))       tleaf_sun240(:,:)        = fillreal
    if(allocated(tleaf_shade240))     tleaf_shade240(:,:)      = fillreal
    if(allocated(tleaf_ave240))       tleaf_ave240(:,:)        = fillreal
    if(allocated(ppfd_sun240))        ppfd_sun240(:,:)         = fillreal
    if(allocated(ppfd_shade240))      ppfd_shade240(:,:)       = fillreal
    if(allocated(daily_maxt2m))       daily_maxt2m(:)          = fillreal
    if(allocated(daily_mint2m))       daily_mint2m(:)          = fillreal
    if(allocated(daily_maxws10m))     daily_maxws10m(:)          = fillreal

    if(allocated(tleaf_sun24_tmp_3d))    tleaf_sun24_tmp_3d(:,:,:,:)     = 0.0_rk
    if(allocated(tleaf_shade24_tmp_3d))  tleaf_shade24_tmp_3d(:,:,:,:)   =  0.0_rk
    if(allocated(tleaf_ave24_tmp_3d))    tleaf_ave24_tmp_3d(:,:,:,:)     = 0.0_rk
    if(allocated(ppfd_sun24_tmp_3d))     ppfd_sun24_tmp_3d(:,:,:,:)      = 0.0_rk
    if(allocated(ppfd_shade24_tmp_3d))   ppfd_shade24_tmp_3d(:,:,:,:)    = 0.0_rk
    if(allocated(tleaf_sun240_tmp_3d))   tleaf_sun240_tmp_3d(:,:,:,:)    = 0.0_rk
    if(allocated(tleaf_shade240_tmp_3d)) tleaf_shade240_tmp_3d(:,:,:,:)  = 0.0_rk
    if(allocated(tleaf_ave240_tmp_3d))   tleaf_ave240_tmp_3d(:,:,:,:)    = 0.0_rk
    if(allocated(ppfd_sun240_tmp_3d))    ppfd_sun240_tmp_3d(:,:,:,:)     = 0.0_rk
    if(allocated(ppfd_shade240_tmp_3d))  ppfd_shade240_tmp_3d(:,:,:,:)   = 0.0_rk
    if(allocated(tmp2mref_tmp_3d))       tmp2mref_tmp_3d(:,:,:)          = 0.0_rk
    if(allocated(ubzref_tmp_3d))         ubzref_tmp_3d(:,:,:)            = 0.0_rk
    if(allocated(tleaf_sun24_3d))        tleaf_sun24_3d(:,:,:)         = fillreal
    if(allocated(tleaf_shade24_3d))      tleaf_shade24_3d(:,:,:)       = fillreal
    if(allocated(tleaf_ave24_3d))        tleaf_ave24_3d(:,:,:)         = fillreal
    if(allocated(ppfd_sun24_3d))         ppfd_sun24_3d(:,:,:)          = fillreal
    if(allocated(ppfd_shade24_3d))       ppfd_shade24_3d(:,:,:)        = fillreal
    if(allocated(tleaf_sun240_3d))       tleaf_sun240_3d(:,:,:)        = fillreal
    if(allocated(tleaf_shade240_3d))     tleaf_shade240_3d(:,:,:)      = fillreal
    if(allocated(tleaf_ave240_3d))       tleaf_ave240_3d(:,:,:)        = fillreal
    if(allocated(ppfd_sun240_3d))        ppfd_sun240_3d(:,:,:)         = fillreal
    if(allocated(ppfd_shade240_3d))      ppfd_shade240_3d(:,:,:)       = fillreal
    if(allocated(daily_maxt2m_2d))       daily_maxt2m_2d(:,:)          = fillreal
    if(allocated(daily_mint2m_2d))       daily_mint2m_2d(:,:)          = fillreal
    if(allocated(daily_maxws10m_2d))     daily_maxws10m_2d(:,:)        = fillreal
!-------------------------------------------------------------------------------
! Initialize arrays for Canopy Wind
!-------------------------------------------------------------------------------

    if (ifcanwind .or. ifcanwaf) then
        if(allocated(canBOT))        canBOT(:)            = fillreal
        if(allocated(canTOP))        canTOP(:)            = fillreal
        if(allocated(canWIND))       canWIND(:,:)         = fillreal
        if(allocated(canWIND_3d))    canWIND_3d(:,:,:)    = fillreal
        if(allocated(dx))            dx(:)                = fillreal
        if(allocated(dx_2d))         dx_2d(:,:)           = fillreal
        if(allocated(flameh))        flameh(:)            = fillreal
        if(allocated(flameh_2d))     flameh_2d(:,:)       = fillreal
        if(allocated(waf))           waf(:)               = fillreal
        if(allocated(waf_2d))        waf_2d(:,:)          = fillreal
    end if

!-------------------------------------------------------------------------------
! Initialize arrays for Canopy Diffusivity Profile
!-------------------------------------------------------------------------------

    if (ifcaneddy) then
        if(allocated(Kz))            Kz(:,:)      = fillreal
        if(allocated(Kz_3d))         Kz_3d(:,:,:) = fillreal
    end if


!-------------------------------------------------------------------------------
! Initialize arrays for Canopy Photolysis Correction Factor
!-------------------------------------------------------------------------------

    if (ifcanphot) then
        if(allocated(rjcf))            rjcf(:,:)      = fillreal
        if(allocated(rjcf_3d))         rjcf_3d(:,:,:) = fillreal
    end if

!-------------------------------------------------------------------------------
! Initialize arrays for Canopy Biogenic Emissions
!-------------------------------------------------------------------------------

    if (ifcanbio) then
        if(allocated(emi_isop))            emi_isop(:,:)      = fillreal
        if(allocated(emi_isop_3d))         emi_isop_3d(:,:,:) = fillreal
        if(allocated(emi_myrc))            emi_myrc(:,:)      = fillreal
        if(allocated(emi_myrc_3d))         emi_myrc_3d(:,:,:) = fillreal
        if(allocated(emi_sabi))            emi_sabi(:,:)      = fillreal
        if(allocated(emi_sabi_3d))         emi_sabi_3d(:,:,:) = fillreal
        if(allocated(emi_limo))            emi_limo(:,:)      = fillreal
        if(allocated(emi_limo_3d))         emi_limo_3d(:,:,:) = fillreal
        if(allocated(emi_care))            emi_care(:,:)      = fillreal
        if(allocated(emi_care_3d))         emi_care_3d(:,:,:) = fillreal
        if(allocated(emi_ocim))            emi_ocim(:,:)      = fillreal
        if(allocated(emi_ocim_3d))         emi_ocim_3d(:,:,:) = fillreal
        if(allocated(emi_bpin))            emi_bpin(:,:)      = fillreal
        if(allocated(emi_bpin_3d))         emi_bpin_3d(:,:,:) = fillreal
        if(allocated(emi_apin))            emi_apin(:,:)      = fillreal
        if(allocated(emi_apin_3d))         emi_apin_3d(:,:,:) = fillreal
        if(allocated(emi_mono))            emi_mono(:,:)      = fillreal
        if(allocated(emi_mono_3d))         emi_mono_3d(:,:,:) = fillreal
        if(allocated(emi_farn))            emi_farn(:,:)      = fillreal
        if(allocated(emi_farn_3d))         emi_farn_3d(:,:,:) = fillreal
        if(allocated(emi_cary))            emi_cary(:,:)      = fillreal
        if(allocated(emi_cary_3d))         emi_cary_3d(:,:,:) = fillreal
        if(allocated(emi_sesq))            emi_sesq(:,:)      = fillreal
        if(allocated(emi_sesq_3d))         emi_sesq_3d(:,:,:) = fillreal
        if(allocated(emi_mbol))            emi_mbol(:,:)      = fillreal
        if(allocated(emi_mbol_3d))         emi_mbol_3d(:,:,:) = fillreal
        if(allocated(emi_meth))            emi_meth(:,:)      = fillreal
        if(allocated(emi_meth_3d))         emi_meth_3d(:,:,:) = fillreal
        if(allocated(emi_acet))            emi_acet(:,:)      = fillreal
        if(allocated(emi_acet_3d))         emi_acet_3d(:,:,:) = fillreal
        if(allocated(emi_co))              emi_co(:,:)        = fillreal
        if(allocated(emi_co_3d))           emi_co_3d(:,:,:)   = fillreal
        if(allocated(emi_bvoc))            emi_bvoc(:,:)      = fillreal
        if(allocated(emi_bvoc_3d))         emi_bvoc_3d(:,:,:) = fillreal
        if(allocated(emi_svoc))            emi_svoc(:,:)      = fillreal
        if(allocated(emi_svoc_3d))         emi_svoc_3d(:,:,:) = fillreal
        if(allocated(emi_ovoc))            emi_ovoc(:,:)      = fillreal
        if(allocated(emi_ovoc_3d))         emi_ovoc_3d(:,:,:) = fillreal
    end if

!-------------------------------------------------------------------------------
! Initialize arrays for Canopy Gas Dry Deposition
!-------------------------------------------------------------------------------

    if (ifcanddepgas) then
        if(allocated(ddep_o3))             ddep_o3(:,:)      = fillreal
    end if

END SUBROUTINE canopy_init
