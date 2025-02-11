
SUBROUTINE canopy_calcs(nn)

!-------------------------------------------------------------------------------
! Name:     Contains main canopy calculations dependent on canopy conditions
! Purpose:  Contains main canopy calculations dependent on canopy conditions
! Revised:  06 Oct 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------

    use canopy_const_mod      !constants for canopy models
    use canopy_coord_mod      !main canopy coordinate descriptions
    use canopy_canopts_mod    !main canopy option descriptions
    use canopy_canmet_mod     !main canopy met/sfc input descriptions
    use canopy_canvars_mod    !main canopy variables descriptions
    use canopy_utils_mod      !main canopy utilities
    use canopy_dxcalc_mod     !main canopy dx calculation
    use canopy_profile_mod    !main canopy foliage profile routines
    use canopy_var3din_mod    !main canopy 3D variable in routines
    use canopy_rad_mod        !main canopy radiation sunlit/shaded routines
    use canopy_tleaf_mod      !main canopy leaf temperature sunlit/shaded routines
    use canopy_wind_mod       !main canopy components
    use canopy_fire_mod
    use canopy_phot_mod
    use canopy_eddy_mod
    use canopy_bioemi_mod
    use canopy_drydep_mod

    IMPLICIT NONE

    INTEGER,     INTENT( IN )       :: nn         ! Input time step

    !Local variables
    integer i,j,k,loc
    ! LAI variables for Leaf Age factor calculations
    REAL(rk) :: pastlai             !Past LAI [cm2/cm2]
    REAL(rk), save :: currentlai    ! Current LAI [cm2/cm2]  (saved from one timestep to the next)
    REAL(rk) :: tsteplai            !Number of days between the past and current LAI
    !Historical Averaging variables for biogenics
    REAL(rk) :: dnewfrac,doldfrac,hnewfrac,holdfrac
    !Other
    REAL(rk) :: lat2d(nlon,nlat), lon2d(nlon,nlat), lat1d(nlon*nlat), lon1d(nlon*nlat)

    write(*,*)  'Calculating Canopy Parameters'
    write(*,*)  '-------------------------------'

    lat2d = variables_2d%lat
    lon2d = variables_2d%lon
    lat1d = variables%lat
    lon1d = variables%lon

    if (infmt_opt .eq. 0) then !Main input format is 2D NetCDF and output will be 2D NetCDf

        if (ifcanwind .or. ifcanwaf) then !only calculate if canopy wind or WAF option
            call canopy_calcdx_2d(dx_opt, dx_set, nlat, nlon, lat2d, &
                lon2d, dx_2d)
        end if

        if (href_opt .eq. 0 ) then !setting entire array = href_set value from user NL
            variables_2d%href = href_set
        else if (href_opt .eq. 1 ) then !from file array
            variables_2d%href =  variables_2d%href
        else
            write(*,*)  'Wrong HREF_OPT choice of ', href_opt, ' in namelist...exiting'
            call exit(2)
        end if

! ... Main loop through model grid cells for 2D input/output

        do i=1, nlon
            do j=1, nlat

                hcmref         = variables_2d(i,j)%ch
                uref           = variables_2d(i,j)%ugrd10m
                vref           = variables_2d(i,j)%vgrd10m
                cluref         = variables_2d(i,j)%clu
                lairef         = variables_2d(i,j)%lai
                vtyperef       = variables_2d(i,j)%vtype
                canfracref     = variables_2d(i,j)%canfrac
                ustref         = variables_2d(i,j)%fricv
                cszref         = variables_2d(i,j)%csz
                z0ref          = variables_2d(i,j)%sfcr
                molref         = variables_2d(i,j)%mol
                frpref         = variables_2d(i,j)%frp
                hgtref         = variables_2d(i,j)%href
                sotypref       = variables_2d(i,j)%sotyp
                pressfcref     = variables_2d(i,j)%pressfc
                dswrfref       = variables_2d(i,j)%dswrf
                shtflref       = variables_2d(i,j)%shtfl
                tmpsfcref      = variables_2d(i,j)%tmpsfc
                tmp2mref       = variables_2d(i,j)%tmp2m
                spfh2mref      = variables_2d(i,j)%spfh2m
                hpblref        = variables_2d(i,j)%hpbl
                prate_averef   = variables_2d(i,j)%prate_ave
                soilw1ref      = variables_2d(i,j)%soilw1
                soilw2ref      = variables_2d(i,j)%soilw2
                soilw3ref      = variables_2d(i,j)%soilw3
                soilw4ref      = variables_2d(i,j)%soilw4
                wiltref        = variables_2d(i,j)%wilt
                ozone_w126ref  = variables_2d(i,j)%ozone_w126

! ... calculate wind speed from u and v
                ubzref   = sqrt((uref**2.0) + (vref**2.0))

! ... get scaled canopy model profile and sub-canopy layers
                if (hcmref > 0) then
                    zhc = zk / hcmref
                else
                    zhc = 0
                end if
                cansublays  = min(floor(hcmref/modres),modlays)
                if (cansublays .lt. 1) then !case where model resolution >= canopy height
                    cansublays=1            !only one layer allowed
                end if
! ... check for valid model vegetation types
                if (lu_opt .eq. 0 .or. lu_opt .eq. 1 ) then !VIIRS or MODIS
                    if (vtyperef .gt. 0 .and. vtyperef .le. 10 .or. vtyperef .eq. 12 &
                        .or. vtyperef .eq. 14 .or. vtyperef .eq. 18 &
                        .or. vtyperef .eq. 19) then

! ... check for ssg_opt from user namelist
                        if (vtyperef .ge. 6 .and. vtyperef .le. 10) then !VIIRS/MODIS shrubs/savannas/grasses (SSG) type
                            if (ssg_opt .eq. 0) then !use GEDI inputs for SSG heights (possibly not measured...)
                                hcmref = hcmref
                            else if (ssg_opt .eq. 1) then !user set constant shrubs/savannas/grasslands height
                                hcmref = ssg_set
                                !recalculate
                                zhc         = zk/hcmref
                                cansublays  = min(floor(hcmref/modres),modlays)
                                if (cansublays .lt. 1) then !case where model resolution >= canopy height
                                    cansublays=1            !only one layer allowed
                                end if
                            else
                                write(*,*)  'Wrong SSG_OPT choice of ', ssg_opt, &
                                    ' in namelist...exiting'
                                call exit(2)
                            end if
                        end if

! ... check for crop_opt from user namelist
                        if (vtyperef .eq. 12 .or. vtyperef .eq. 14) then !VIIRS/MODIS crop types
                            if (crop_opt .eq. 0) then !use GEDI inputs for crop height (possibly not measured...)
                                hcmref = hcmref
                            else if (crop_opt .eq. 1) then !user set constant crop height
                                hcmref = crop_set
                                !recalculate
                                zhc         = zk/hcmref
                                cansublays  = min(floor(hcmref/modres),modlays)
                                if (cansublays .lt. 1) then !case where model resolution >= canopy height
                                    cansublays=1            !only one layer allowed
                                end if
                            else
                                write(*,*)  'Wrong CROP_OPT choice of ', crop_opt, &
                                    ' in namelist...exiting'
                                call exit(2)
                            end if
                        end if

! ... check for contiguous canopy conditions at each model grid cell
                        if (hcmref .gt. ch_thresh .and. canfracref .gt. cf_thresh &
                            .and. lairef .gt. lai_thresh) then

! ... call canopy parameters to get canopy, fire info, and shape distribution parameters

                            call canopy_parm(vtyperef, hcmref, canfracref, lairef, &
                                pai_opt, pai_set, lu_opt, firetype, cdrag, &
                                pai, zcanmax, sigmau, sigma1)

! ... Choose between prescribed canopy/foliate shape profile or observed GEDI PAVD profile
                            if (pavd_opt .eq. 0) then
! ... calculate canopy/foliage distribution shape profile - bottom up total in-canopy and fraction at z
                                call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                    fafraczInt)
                            else
                                pavdref = variables_3d(i,j,:)%pavd
                                levref  = variables_1d%lev
! ... derive canopy/foliage distribution shape profile from interpolated GEDI PAVD profile - bottom up total in-canopy and fraction at z
                                if (variables_2d(i,j)%lat .gt. (-1.0_rk*pavd_set) .and. &
                                    variables_2d(i,j)%lat .lt. pavd_set) then !use GEDI PAVD
                                    call canopy_pavd2fafrac(zcanmax, sigmau, sigma1, hcmref, zhc, &
                                        pavdref, levref, fafraczInt)
                                    !check if there is observed canopy height but no PAVD profile
                                    if (hcmref .gt. 0.0 .and. maxval(fafraczInt) .le. 0.0) then !revert to prescribed shape profile
                                        call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                            fafraczInt)
                                    end if
                                else !revert back to using prescribed shape profile
                                    call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                        fafraczInt)
                                end if
                            end if

! ... calculate initial canopy temp/pressure/humidity/density profiles from first guess approximations (i.e., no leaf energy balance)

                            do k=1, modlays
                                tka_3d(i,j,k)=CalcTemp(zk(k)*100.0_rk, 200.0_rk, tmp2mref-273.15_rk, tmpsfcref-273.15_rk) ! temp    [K]
                                pressa_3d(i,j,k)=CalcPressure(zk(k)*100.0_rk, 200.0_rk, pressfcref*0.01_rk, &
                                    tmp2mref, tmpsfcref)                                                                  ! press   [mb]
                                relhuma_3d(i,j,k)=CalcRelHum(tka_3d(i,j,k),pressa_3d(i,j,k),spfh2mref*1000.0_rk)    ! relhum  [%]
                                spechuma_3d(i,j,k)=CalcSpecHum(relhuma_3d(i,j,k),tka_3d(i,j,k),pressa_3d(i,j,k))     ! spechum [g/kg]
                            end do

! ... calculate leaf area density profile from foliage shape function for output (m2/m3)
                            do k=1, modlays
                                if (zk(k) .gt. 0.0 .and. zk(k) .le. hcmref) then ! above ground level and at/below canopy top
                                    if (k .lt. modlays)  then
                                        lad_3d(i,j,k) = ((fafraczInt(k+1) - fafraczInt(k))*lairef)/modres
                                    else
                                        lad_3d(i,j,k) = lad_3d(i,j,modlays-1)
                                    end if
                                else
                                    lad_3d(i,j,k) = 0.0_rk
                                end if
                            end do
! ... calculate zero-plane displacement height/hc and surface (soil+veg) roughness lengths/hc
                            call canopy_zpd(zhc(1:cansublays), fafraczInt(1:cansublays), &
                                ubzref, z0ghc, lambdars, cdrag, pai, hcmref, hgtref, &
                                z0ref, vtyperef, lu_opt, z0_opt, d_h, zo_h)

! ... calculate canopy radiation (sunlit and shade) profile

                            call canopy_fsun_clu( fafraczInt, lairef, cluref, cszref, fsun)

! ... calculate canopy leaf temperature (sun/shade) profile

                            call canopy_tleaf_lin(zk, hcmref, tmp2mref, fsun, &
                                tleaf_sun, tleaf_shade, tleaf_ave)

! ... calculate canopy Photosynthetic Photon Flux Density (PPFD) (sun/shade) profile

                            call canopy_ppfd_exp(zk, hcmref, dswrfref, lairef, fsun, &
                                ppfd_sun, ppfd_shade, ppfd_ave)

! ... user option to calculate in-canopy wind speeds at height z and midflame WAF

                            if (ifcanwind .or. ifcanwaf) then
                                if (rsl_opt .eq. 0) then
                                    do k=1, modlays
                                        call canopy_wind_most(hcmref, zk(k), fafraczInt(k), ubzref, &
                                            z0ghc, cdrag, pai, hgtref, d_h, zo_h, &
                                            lambdars, canBOT(k), canTOP(k), canWIND_3d(i,j,k))
                                    end do
                                else
                                    write(*,*) 'wrong RSL_OPT namelist option = ', rsl_opt, 'only option = 0 right now'
                                    call exit(2)
                                end if

! ... determine midflamepoint and flame height from user or FRP calculation
                                call canopy_flameh(flameh_opt, flameh_set, dx_2d(i,j), modres, &
                                    frpref, frp_fac, hcmref, lu_opt, vtyperef, flameh_cal, &
                                    midflamepoint, flameh_2d(i,j))
                                if (firetype .eq. 0) then !forest/sub-canopy firetype
                                    if (flameh_2d(i,j) .gt. 0.0) then !flameh must be > 0
                                        if (flameh_2d(i,j) .le. hcmref) then !only calculate when flameh <= FCH
                                            call canopy_waf(hcmref, lambdars, hgtref, flameh_2d(i,j), &
                                                firetype, d_h, zo_h, canBOT(midflamepoint), &
                                                canTOP(midflamepoint), waf_2d(i,j))
                                        else
                                            write(*,*) 'warning...sub-canopy type fire, but flameh > FCH, setting WAF=1'
                                            waf_2d(i,j) = 1.0_rk
                                        end if
                                    end if
                                else  !grass/crops, above-canopy firetype
                                    if (flameh_2d(i,j) .gt. 0.0) then !flameh still must be > 0
                                        call canopy_waf(hcmref, lambdars, hgtref, flameh_2d(i,j), &
                                            firetype, d_h, zo_h, canBOT(midflamepoint), &
                                            canTOP(midflamepoint), waf_2d(i,j))
                                    end if
                                end if
                                if (waf_2d(i,j) .gt. 1.0_rk) then !Final check of WAF > 1, must be <=1
                                    waf_2d(i,j) = 1.0_rk
                                end if
                            end if

! ... user option to calculate in-canopy eddy diffusivities at height z
                            if (ifcaneddy) then
                                do k=1, modlays
                                    call canopy_eddyx(hcmref, zk(k), ustref, molref, Kz_3d(i,j,k))
                                end do
                            end if

! ... user option to calculate in-canopy eddy photolysis attenuation at height z
                            if (ifcanphot) then
                                if (cszref .ge. 0.0_rk) then !only calculate if cell isn't dark
                                    call canopy_phot(fafraczInt, &
                                        lairef, cluref, cszref, rjcf_3d(i,j,:))
                                end if
                            end if

!.......user option to calculate in-canopy leafage influence and assigning LAI as per timestep


                            if (leafage_opt .eq. 0) then
                                ! Initialize pastlai and currentlai based on current timestep
                                if (nn .eq. 1) then
                                    currentlai = lairef
                                    pastlai = currentlai
                                else
                                    pastlai = currentlai
                                    currentlai = lairef
                                end if

                                !!! Check if the lai_tstep is greater than time_intvl
                                if (lai_tstep .ge. time_intvl) then
                                    tsteplai = lai_tstep/86400.0_rk
                                else
                                    WRITE (*, *) "Error: Input LAI time step cannot be less than model time step...exiting!!!"
                                    CALL EXIT(1)
                                endif
                            end if !leafage_opt = 0 end

!.......user option to calculate historical leaf temperature and PAR for past 24-hours and 240-hours rolling average per timestep
                            !Initialize

                            if (hist_opt .eq. 0) then      !Use instantaneous only
                                ppfd_sun24_3d(i,j,:)     = ppfd_sun
                                ppfd_shade24_3d(i,j,:)   = ppfd_shade
                                tleaf_sun24_3d(i,j,:)    = tleaf_sun
                                tleaf_shade24_3d(i,j,:)  = tleaf_shade
                                tleaf_ave24_3d(i,j,:)    = tleaf_ave
                                ppfd_sun240_3d(i,j,:)    = ppfd_sun
                                ppfd_shade240_3d(i,j,:)  = ppfd_shade
                                tleaf_sun240_3d(i,j,:)   = tleaf_sun
                                tleaf_shade240_3d(i,j,:) = tleaf_shade
                                tleaf_ave240_3d(i,j,:)   = tleaf_ave
                                !for AQ, temp, and wind stress factors in biogenics, use instantaneous only (caution!)
                                daily_maxt2m_2d(i,j)     = tmp2mref
                                daily_mint2m_2d(i,j)     = tmp2mref
                                daily_maxws10m_2d(i,j)   = ubzref
                            else if (hist_opt .eq. 1) then  !Try for historical average values
                                ! Calculate weights for running means of historic variables
                                ! DNEWFRAC and DOLDFRAC are the weights given to the current
                                ! and existing value, respectively, when updating running means
                                ! over the last X days. HNEWFRAC and HOLDFRAC are the same but
                                ! for the 24H means.
                                dnewfrac = time_intvl / ( tau_days * 24.0_rk * 3600.0_rk )
                                doldfrac = 1.0_rk - dnewfrac
                                hnewfrac = time_intvl / ( tau_hours * 3600.0_rk )
                                holdfrac = 1.0_rk - hnewfrac

                                if (nn .le. 24) then
                                    !Track times for moving time window (running) average below
                                    ppfd_sun24_tmp_3d(nn,i,j,:)     = ppfd_sun
                                    ppfd_shade24_tmp_3d(nn,i,j,:)   = ppfd_shade
                                    tleaf_sun24_tmp_3d(nn,i,j,:)    = tleaf_sun
                                    tleaf_shade24_tmp_3d(nn,i,j,:)  = tleaf_shade
                                    tleaf_ave24_tmp_3d(nn,i,j,:)    = tleaf_ave
                                    ppfd_sun240_tmp_3d(nn,i,j,:)    = ppfd_sun
                                    ppfd_shade240_tmp_3d(nn,i,j,:)  = ppfd_shade
                                    tleaf_sun240_tmp_3d(nn,i,j,:)   = tleaf_sun
                                    tleaf_shade240_tmp_3d(nn,i,j,:) = tleaf_shade
                                    tleaf_ave240_tmp_3d(nn,i,j,:)   = tleaf_ave
                                    tmp2mref_tmp_3d(nn,i,j)         = tmp2mref
                                    ubzref_tmp_3d(nn,i,j)           = ubzref

                                    !!TODO:  Restart capability needed to get past leaf temp and PAR if avaialble
                                    !For now, if <= 24 hours then only option is to use current instantaneous values
                                    ppfd_sun24_3d(i,j,:)     = ppfd_sun
                                    ppfd_shade24_3d(i,j,:)   = ppfd_shade
                                    tleaf_sun24_3d(i,j,:)    = tleaf_sun
                                    tleaf_shade24_3d(i,j,:)  = tleaf_shade
                                    tleaf_ave24_3d(i,j,:)    = tleaf_ave
                                    ppfd_sun240_3d(i,j,:)    = ppfd_sun
                                    ppfd_shade240_3d(i,j,:)  = ppfd_shade
                                    tleaf_sun240_3d(i,j,:)   = tleaf_sun
                                    tleaf_shade240_3d(i,j,:) = tleaf_shade
                                    tleaf_ave240_3d(i,j,:)   = tleaf_ave
                                    !for AQ, temp, and wind stress factors in biogenics, still use instantaneous only (caution!)
                                    daily_maxt2m_2d(i,j)        = tmp2mref
                                    daily_mint2m_2d(i,j)        = tmp2mref
                                    daily_maxws10m_2d(i,j)      = ubzref

                                else if (nn .le. 240) then
                                    !Drop nth value into end time slice
                                    ppfd_sun24_tmp_3d(25,i,j,:)     = ppfd_sun
                                    ppfd_shade24_tmp_3d(25,i,j,:)   = ppfd_shade
                                    tleaf_sun24_tmp_3d(25,i,j,:)    = tleaf_sun
                                    tleaf_shade24_tmp_3d(25,i,j,:)  = tleaf_shade
                                    tleaf_ave24_tmp_3d(25,i,j,:)    = tleaf_ave
                                    ppfd_sun240_tmp_3d(nn,i,j,:)    = ppfd_sun
                                    ppfd_shade240_tmp_3d(nn,i,j,:)  = ppfd_shade
                                    tleaf_sun240_tmp_3d(nn,i,j,:)   = tleaf_sun
                                    tleaf_shade240_tmp_3d(nn,i,j,:) = tleaf_shade
                                    tleaf_ave240_tmp_3d(nn,i,j,:)   = tleaf_ave
                                    tmp2mref_tmp_3d(25,i,j)         = tmp2mref
                                    ubzref_tmp_3d(25,i,j)           = ubzref

                                    !Calculate 24 and 240 hour averages
                                    ppfd_sun24_3d(i,j,:)     = sum(ppfd_sun24_tmp_3d(1:24,i,j,:),1)/24.0
                                    ppfd_shade24_3d(i,j,:)   = sum(ppfd_shade24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_sun24_3d(i,j,:)    = sum(tleaf_sun24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_shade24_3d(i,j,:)  = sum(tleaf_shade24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_ave24_3d(i,j,:)    = sum(tleaf_ave24_tmp_3d(1:24,i,j,:),1)/24.0
                                    !Use 24 hour data until 240 data become available
                                    ppfd_sun240_3d(i,j,:)    = sum(ppfd_sun24_tmp_3d(1:24,i,j,:),1)/24.0
                                    ppfd_shade240_3d(i,j,:)  = sum(ppfd_shade24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_sun240_3d(i,j,:)   = sum(tleaf_sun24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_shade240_3d(i,j,:) = sum(tleaf_shade24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_ave240_3d(i,j,:)   = sum(tleaf_ave24_tmp_3d(1:24,i,j,:),1)/24.0
                                    !Updated rolling 24 hour (hourly, short term) and 240 hour (daily, long-term) averages
                                    ppfd_sun24_3d(i,j,:)     = ( holdfrac * ppfd_sun24_3d(i,j,:) ) + ( hnewfrac * ppfd_sun )
                                    ppfd_shade24_3d(i,j,:)   = ( holdfrac * ppfd_shade24_3d(i,j,:) ) + ( hnewfrac * ppfd_shade )
                                    tleaf_sun24_3d(i,j,:)    = ( holdfrac * tleaf_sun24_3d(i,j,:) ) + ( hnewfrac * tleaf_sun )
                                    tleaf_shade24_3d(i,j,:)  = ( holdfrac * tleaf_shade24_3d(i,j,:) ) + ( hnewfrac * tleaf_shade )
                                    tleaf_ave24_3d(i,j,:)    = ( holdfrac * tleaf_ave24_3d(i,j,:) ) + ( hnewfrac * tleaf_ave )
                                    ppfd_sun240_3d(i,j,:)     = ( holdfrac * ppfd_sun24_3d(i,j,:) ) + ( hnewfrac * ppfd_sun )
                                    ppfd_shade240_3d(i,j,:)   = ( holdfrac * ppfd_shade24_3d(i,j,:) ) + ( hnewfrac * ppfd_shade )
                                    tleaf_sun240_3d(i,j,:)    = ( holdfrac * tleaf_sun24_3d(i,j,:) ) + ( hnewfrac * tleaf_sun )
                                    tleaf_shade240_3d(i,j,:)  = ( holdfrac * tleaf_shade24_3d(i,j,:) ) + ( hnewfrac * tleaf_shade )
                                    tleaf_ave240_3d(i,j,:)    = ( holdfrac * tleaf_ave24_3d(i,j,:) ) + ( hnewfrac * tleaf_ave )
                                    !Take daily max/min for temperature and daily max for wind speed
                                    daily_maxt2m_2d(i,j)        = maxval(tmp2mref_tmp_3d(1:24,i,j))
                                    daily_mint2m_2d(i,j)        = minval(tmp2mref_tmp_3d(1:24,i,j))
                                    daily_maxws10m_2d(i,j)      = maxval(ubzref_tmp_3d(1:24,i,j))

                                    !Shift time window to make room for new nth value
                                    ppfd_sun24_tmp_3d(1:24,i,j,:)     = ppfd_sun24_tmp_3d(2:25,i,j,:)
                                    ppfd_shade24_tmp_3d(1:24,i,j,:)   = ppfd_shade24_tmp_3d(2:25,i,j,:)
                                    tleaf_sun24_tmp_3d(1:24,i,j,:)    = tleaf_sun24_tmp_3d(2:25,i,j,:)
                                    tleaf_shade24_tmp_3d(1:24,i,j,:)  = tleaf_shade24_tmp_3d(2:25,i,j,:)
                                    tleaf_ave24_tmp_3d(1:24,i,j,:)    = tleaf_ave24_tmp_3d(2:25,i,j,:)
                                    tmp2mref_tmp_3d(1:24,i,j)         = tmp2mref_tmp_3d(2:25,i,j)
                                    ubzref_tmp_3d(1:24,i,j)           = ubzref_tmp_3d(2:25,i,j)

                                else
                                    !Drop nth value into end time slice
                                    ppfd_sun24_tmp_3d(25,i,j,:)     = ppfd_sun
                                    ppfd_shade24_tmp_3d(25,i,j,:)   = ppfd_shade
                                    tleaf_sun24_tmp_3d(25,i,j,:)    = tleaf_sun
                                    tleaf_shade24_tmp_3d(25,i,j,:)  = tleaf_shade
                                    tleaf_ave24_tmp_3d(25,i,j,:)    = tleaf_ave
                                    ppfd_sun240_tmp_3d(241,i,j,:)    = ppfd_sun
                                    ppfd_shade240_tmp_3d(241,i,j,:)  = ppfd_shade
                                    tleaf_sun240_tmp_3d(241,i,j,:)   = tleaf_sun
                                    tleaf_shade240_tmp_3d(241,i,j,:) = tleaf_shade
                                    tleaf_ave240_tmp_3d(241,i,j,:)   = tleaf_ave
                                    tmp2mref_tmp_3d(25,i,j)         = tmp2mref
                                    ubzref_tmp_3d(25,i,j)           = ubzref

                                    !Calculate 24 and 240 hour averages
                                    ppfd_sun24_3d(i,j,:)     = sum(ppfd_sun24_tmp_3d(1:24,i,j,:),1)/24.0
                                    ppfd_shade24_3d(i,j,:)   = sum(ppfd_shade24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_sun24_3d(i,j,:)    = sum(tleaf_sun24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_shade24_3d(i,j,:)  = sum(tleaf_shade24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_ave24_3d(i,j,:)    = sum(tleaf_ave24_tmp_3d(1:24,i,j,:),1)/24.0
                                    ppfd_sun240_3d(i,j,:)    = sum(ppfd_sun240_tmp_3d(1:240,i,j,:),1)/240.0
                                    ppfd_shade240_3d(i,j,:)  = sum(ppfd_shade240_tmp_3d(1:240,i,j,:),1)/240.0
                                    tleaf_sun240_3d(i,j,:)   = sum(tleaf_sun240_tmp_3d(1:240,i,j,:),1)/240.0
                                    tleaf_shade240_3d(i,j,:) = sum(tleaf_shade240_tmp_3d(1:240,i,j,:),1)/240.0
                                    tleaf_ave240_3d(i,j,:)   = sum(tleaf_ave240_tmp_3d(1:240,i,j,:),1)/240.0
                                    !Update for current time value using efolding (holdfrac and hnewfrac)
                                    ppfd_sun24_3d(i,j,:)     = ( holdfrac * ppfd_sun24_3d(i,j,:) ) + ( hnewfrac * ppfd_sun )
                                    ppfd_shade24_3d(i,j,:)   = ( holdfrac * ppfd_shade24_3d(i,j,:) ) + ( hnewfrac * ppfd_shade )
                                    tleaf_sun24_3d(i,j,:)    = ( holdfrac * tleaf_sun24_3d(i,j,:) ) + ( hnewfrac * tleaf_sun )
                                    tleaf_shade24_3d(i,j,:)  = ( holdfrac * tleaf_shade24_3d(i,j,:) ) + ( hnewfrac * tleaf_shade )
                                    tleaf_ave24_3d(i,j,:)    = ( holdfrac * tleaf_ave24_3d(i,j,:) ) + ( hnewfrac * tleaf_ave )
                                    ppfd_sun240_3d(i,j,:)    = ( doldfrac * ppfd_sun240_3d(i,j,:) ) + ( dnewfrac * ppfd_sun )
                                    ppfd_shade240_3d(i,j,:)  = ( doldfrac * ppfd_shade240_3d(i,j,:) ) + ( dnewfrac * ppfd_shade )
                                    tleaf_sun240_3d(i,j,:)   = ( doldfrac * tleaf_sun240_3d(i,j,:) ) + ( dnewfrac * tleaf_sun )
                                    tleaf_shade240_3d(i,j,:) = ( doldfrac * tleaf_shade240_3d(i,j,:) ) + ( dnewfrac * tleaf_shade )
                                    tleaf_ave240_3d(i,j,:)   = ( doldfrac * tleaf_ave240_3d(i,j,:) ) + ( dnewfrac * tleaf_ave )
                                    !Take daily max/min for temperature and daily max for wind speed
                                    daily_maxt2m_2d(i,j)        = maxval(tmp2mref_tmp_3d(1:24,i,j))
                                    daily_mint2m_2d(i,j)        = minval(tmp2mref_tmp_3d(1:24,i,j))
                                    daily_maxws10m_2d(i,j)      = maxval(ubzref_tmp_3d(1:24,i,j))

                                    !Shift time window to make room for new nth value
                                    ppfd_sun24_tmp_3d(1:24,i,j,:)     = ppfd_sun24_tmp_3d(2:25,i,j,:)
                                    ppfd_shade24_tmp_3d(1:24,i,j,:)   = ppfd_shade24_tmp_3d(2:25,i,j,:)
                                    tleaf_sun24_tmp_3d(1:24,i,j,:)    = tleaf_sun24_tmp_3d(2:25,i,j,:)
                                    tleaf_shade24_tmp_3d(1:24,i,j,:)  = tleaf_shade24_tmp_3d(2:25,i,j,:)
                                    tleaf_ave24_tmp_3d(1:24,i,j,:)    = tleaf_ave24_tmp_3d(2:25,i,j,:)
                                    ppfd_sun240_tmp_3d(1:240,i,j,:)    = ppfd_sun240_tmp_3d(2:241,i,j,:)
                                    ppfd_shade240_tmp_3d(1:240,i,j,:)  = ppfd_shade240_tmp_3d(2:241,i,j,:)
                                    tleaf_sun240_tmp_3d(1:240,i,j,:)   = tleaf_sun240_tmp_3d(2:241,i,j,:)
                                    tleaf_shade240_tmp_3d(1:240,i,j,:) = tleaf_shade240_tmp_3d(2:241,i,j,:)
                                    tleaf_ave240_tmp_3d(1:240,i,j,:)   = tleaf_ave240_tmp_3d(2:241,i,j,:)
                                    tmp2mref_tmp_3d(1:24,i,j)         = tmp2mref_tmp_3d(2:25,i,j)
                                    ubzref_tmp_3d(1:24,i,j)           = ubzref_tmp_3d(2:25,i,j)
                                end if
                            else
                                write(*,*) 'wrong HIST_OPT namelist option = ', hist_opt, 'only option = 0 or 1'
                                call exit(2)
                            end if
! ... user option to calculate in-canopy biogenic emissions
                            if (ifcanbio) then
                                if (cszref .ge. 0.0_rk .and. dswrfref .gt. 0.0_rk &
                                    .and. cluref .gt. 0.0_rk) then
                                    !ISOP
                                    if (biospec_opt == 0 .or. biospec_opt == 1) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 1, emi_isop_3d(i,j,:))
                                    end if
                                    !MYRC
                                    if (biospec_opt == 0 .or. biospec_opt == 2) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 2, emi_myrc_3d(i,j,:))
                                    end if
                                    !SABI
                                    if (biospec_opt == 0 .or. biospec_opt == 3) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 3, emi_sabi_3d(i,j,:))
                                    end if
                                    !LIMO
                                    if (biospec_opt == 0 .or. biospec_opt == 4) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 4, emi_limo_3d(i,j,:))
                                    end if
                                    !CARE
                                    if (biospec_opt == 0 .or. biospec_opt == 5) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 5, emi_care_3d(i,j,:))
                                    end if
                                    !OCIM
                                    if (biospec_opt == 0 .or. biospec_opt == 6) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 6, emi_ocim_3d(i,j,:))
                                    end if
                                    !BPIN
                                    if (biospec_opt == 0 .or. biospec_opt == 7) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 7, emi_bpin_3d(i,j,:))
                                    end if
                                    !APIN
                                    if (biospec_opt == 0 .or. biospec_opt == 8) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 8, emi_apin_3d(i,j,:))
                                    end if
                                    !MONO
                                    if (biospec_opt == 0 .or. biospec_opt == 9) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 9, emi_mono_3d(i,j,:))
                                    end if
                                    !FARN
                                    if (biospec_opt == 0 .or. biospec_opt == 10) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 10, emi_farn_3d(i,j,:))
                                    end if
                                    !CARY
                                    if (biospec_opt == 0 .or. biospec_opt == 11) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 11, emi_cary_3d(i,j,:))
                                    end if
                                    !SESQ
                                    if (biospec_opt == 0 .or. biospec_opt == 12) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 12, emi_sesq_3d(i,j,:))
                                    end if
                                    !MBOL
                                    if (biospec_opt == 0 .or. biospec_opt == 13) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 13, emi_mbol_3d(i,j,:))
                                    end if
                                    !METH
                                    if (biospec_opt == 0 .or. biospec_opt == 14) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 14, emi_meth_3d(i,j,:))
                                    end if
                                    !ACET
                                    if (biospec_opt == 0 .or. biospec_opt == 15) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 15, emi_acet_3d(i,j,:))
                                    end if
                                    !CO
                                    if (biospec_opt == 0 .or. biospec_opt == 16) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 16, emi_co_3d(i,j,:))
                                    end if
                                    !BIDI VOC
                                    if (biospec_opt == 0 .or. biospec_opt == 17) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 17, emi_bvoc_3d(i,j,:))
                                    end if
                                    !Stress VOC
                                    if (biospec_opt == 0 .or. biospec_opt == 18) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 18, emi_svoc_3d(i,j,:))
                                    end if
                                    !Other VOC
                                    if (biospec_opt == 0 .or. biospec_opt == 19) then
                                        call canopy_bio(zk, fafraczInt, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 19, emi_ovoc_3d(i,j,:))
                                    end if
                                end if
                            end if

                        end if !Contiguous Canopy

                    else
!                        write(*,*)  'Warning VIIRS/MODIS VTYPE ', vtyperef, ' is not supported...continue'
                    end if   !Vegetation types
                else
                    write(*,*)  'Wrong LU_OPT choice of ', lu_opt, ' in namelist...exiting'
                    call exit(2)

                end if       !Landuse Options
            end do !lat
        end do   !lon

!----------------------------------------------------------->

    else if (infmt_opt .eq. 1) then !Main input format is 1D/2D text and output will be 1D/2D text

        if (ifcanwind .or. ifcanwaf) then !only calculate if canopy wind or WAF option
            call canopy_calcdx(dx_opt, dx_set, nlat, nlon, lat1d, &
                lon1d, dx)
        end if

        if (href_opt .eq. 0 ) then !setting entire array = href_set value from user NL
            variables%href = href_set
        else if (href_opt .eq. 1 ) then !from file array
            variables%href =  variables%href
        else
            write(*,*)  'Wrong HREF_OPT choice of ', href_opt, ' in namelist...exiting'
            call exit(2)
        end if

! ... Main loop through model grid cells
        do loc=1, nlat*nlon
            hcmref         = variables(loc)%ch
            uref           = variables(loc)%ugrd10m
            vref           = variables(loc)%vgrd10m
            cluref         = variables(loc)%clu
            lairef         = variables(loc)%lai
            vtyperef       = variables(loc)%vtype
            canfracref     = variables(loc)%canfrac
            ustref         = variables(loc)%fricv
            cszref         = variables(loc)%csz
            z0ref          = variables(loc)%sfcr
            molref         = variables(loc)%mol
            frpref         = variables(loc)%frp
            hgtref         = variables(loc)%href
            sotypref       = variables(loc)%sotyp
            pressfcref     = variables(loc)%pressfc
            dswrfref       = variables(loc)%dswrf
            shtflref       = variables(loc)%shtfl
            tmpsfcref      = variables(loc)%tmpsfc
            tmp2mref       = variables(loc)%tmp2m
            spfh2mref      = variables(loc)%spfh2m
            hpblref        = variables(loc)%hpbl
            prate_averef   = variables(loc)%prate_ave
            soilw1ref      = variables(loc)%soilw1
            soilw2ref      = variables(loc)%soilw2
            soilw3ref      = variables(loc)%soilw3
            soilw4ref      = variables(loc)%soilw4
            wiltref        = variables(loc)%wilt
            ozone_w126ref  = variables(loc)%ozone_w126

            if (var3d_opt .eq. 1) then !allocated so set
                pavd_arr     = (/variables_can(loc)%pavd01, &
                    variables_can(loc)%pavd02, &
                    variables_can(loc)%pavd03, &
                    variables_can(loc)%pavd04, &
                    variables_can(loc)%pavd05, &
                    variables_can(loc)%pavd06, &
                    variables_can(loc)%pavd07, &
                    variables_can(loc)%pavd08, &
                    variables_can(loc)%pavd09, &
                    variables_can(loc)%pavd10, &
                    variables_can(loc)%pavd11, &
                    variables_can(loc)%pavd12, &
                    variables_can(loc)%pavd13, &
                    variables_can(loc)%pavd14/)
                lev_arr      = (/variables_can(loc)%lev01, &
                    variables_can(loc)%lev02, &
                    variables_can(loc)%lev03, &
                    variables_can(loc)%lev04, &
                    variables_can(loc)%lev05, &
                    variables_can(loc)%lev06, &
                    variables_can(loc)%lev07, &
                    variables_can(loc)%lev08, &
                    variables_can(loc)%lev09, &
                    variables_can(loc)%lev10, &
                    variables_can(loc)%lev11, &
                    variables_can(loc)%lev12, &
                    variables_can(loc)%lev13, &
                    variables_can(loc)%lev14/)
            end if

! ... calculate wind speed from u and v
            ubzref   = sqrt((uref**2.0) + (vref**2.0))

! ... get scaled canopy model profile and sub-canopy layers
            if (hcmref > 0) then
                zhc = zk / hcmref
            else
                zhc = 0
            end if
            cansublays  = min(floor(hcmref/modres),modlays)
            if (cansublays .lt. 1) then !case where model resolution >= canopy height
                cansublays=1            !only one layer allowed
            end if

! ... check for valid model vegetation types
            if (lu_opt .eq. 0 .or. lu_opt .eq. 1 ) then !VIIRS or MODIS
                if (vtyperef .gt. 0 .and. vtyperef .le. 10 .or. vtyperef .eq. 12 &
                    .or. vtyperef .eq. 14 .or. vtyperef .eq. 18 &
                    .or. vtyperef .eq. 19) then

! ... check for ssg_opt from user namelist
                    if (vtyperef .ge. 6 .and. vtyperef .le. 10) then !VIIRS/MODIS shrubs/savannas/grasses (SSG) type
                        if (ssg_opt .eq. 0) then !use GEDI inputs for SSG heights (possibly not measured...)
                            hcmref = hcmref
                        else if (ssg_opt .eq. 1) then !user set constant shrubs/savannas/grasslands height
                            hcmref = ssg_set
                            !recalculate
                            zhc         = zk/hcmref
                            cansublays  = min(floor(hcmref/modres),modlays)
                            if (cansublays .lt. 1) then !case where model resolution >= canopy height
                                cansublays=1            !only one layer allowed
                            end if
                        else
                            write(*,*)  'Wrong SSG_OPT choice of ', ssg_opt, &
                                ' in namelist...exiting'
                            call exit(2)
                        end if
                    end if

! ... check for crop_opt from user namelist
                    if (vtyperef .eq. 12 .or. vtyperef .eq. 14) then !VIIRS/MODIS crop type
                        if (crop_opt .eq. 0) then !use GEDI inputs for crop height
                            hcmref = hcmref
                        else if (crop_opt .eq. 1) then !user set constant crop height
                            hcmref = crop_set
                            !recalculate
                            zhc         = zk/hcmref
                            cansublays  = min(floor(hcmref/modres),modlays)
                            if (cansublays .lt. 1) then !case where model resolution >= canopy height
                                cansublays=1            !only one layer allowed
                            end if
                        else
                            write(*,*)  'Wrong CROP_OPT choice of ', crop_opt, &
                                ' in namelist...exiting'
                            call exit(2)
                        end if
                    end if

! ... check for contiguous canopy conditions at each model grid cell
                    if (hcmref .gt. ch_thresh .and. canfracref .gt. cf_thresh &
                        .and. lairef .gt. lai_thresh) then

! ... call canopy parameters to get canopy, fire info, and shape distribution parameters

                        call canopy_parm(vtyperef, hcmref, canfracref, lairef, &
                            pai_opt, pai_set, lu_opt, firetype, cdrag, &
                            pai, zcanmax, sigmau, sigma1)

! ... Choose between prescribed canopy/foliate shape profile or observed GEDI PAVD profile

                        if (pavd_opt .eq. 0) then
! ... calculate canopy/foliage distribution shape profile - bottom up total in-canopy and fraction at z
                            call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                fafraczInt)
                        else
                            if (var3d_opt .ne. 1) then
                                write(*,*) 'wrong VAR3D_OPT namelist option = ', var3d_opt, &
                                    'change to 1 for supporting PAVD text read'
                                call exit(2)
                            end if
! ... derive canopy/foliage distribution shape profile from interpolated GEDI PAVD profile - bottom up total in-canopy and fraction at z
                            if (variables(loc)%lat .gt. (-1.0_rk*pavd_set) .and. &
                                variables(loc)%lat .lt. pavd_set) then !use GEDI PAVD
                                call canopy_pavd2fafrac(zcanmax, sigmau, sigma1, hcmref, zhc, &
                                    pavd_arr, lev_arr, fafraczInt)
                                !check if there is observed canopy height but no PAVD profile
                                if (hcmref .gt. 0.0 .and. maxval(fafraczInt) .le. 0.0) then !revert to prescribed shape profile
                                    call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                        fafraczInt)
                                end if
                            else !revert back to using prescribed shape profile
                                call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                    fafraczInt)
                            end if
                        end if

! ... calculate initial canopy temp/pressure/humidity/density profiles from first guess approximations (i.e., no leaf energy balance)

                        do k=1, modlays
                            tka(loc,k)=CalcTemp(zk(k)*100.0_rk, 200.0_rk, tmp2mref-273.15_rk, tmpsfcref-273.15_rk) ! temp    [K]
                            pressa(loc,k)=CalcPressure(zk(k)*100.0_rk, 200.0_rk, pressfcref*0.01_rk, &
                                tmp2mref, tmpsfcref)                                                               ! press   [mb]
                            relhuma(loc,k)=CalcRelHum(tka(loc,k),pressa(loc,k),spfh2mref*1000.0_rk)          ! relhum  [%]
                            spechuma(loc,k)=CalcSpecHum(relhuma(loc,k),tka(loc,k),pressa(loc,k))              ! spechum [g/kg]
                        end do

! ... calculate leaf area density profile from foliage shape function for output (m2/m3)
                        do k=1, modlays
                            if (zk(k) .gt. 0.0 .and. zk(k) .le. hcmref) then ! above ground level and at/below canopy top
                                if (k .lt. modlays)  then
                                    lad(loc,k) = ((fafraczInt(k+1) - fafraczInt(k))*lairef)/modres
                                else
                                    lad(loc,k) = lad(loc,modlays-1)
                                end if
                            else
                                lad(loc,k) = 0.0_rk
                            end if
                        end do

! ... calculate zero-plane displacement height/hc and surface (soil+veg) roughness lengths/hc

                        call canopy_zpd(zhc(1:cansublays), fafraczInt(1:cansublays), &
                            ubzref, z0ghc, lambdars, cdrag, pai, hcmref, hgtref, &
                            z0ref, vtyperef, lu_opt, z0_opt, d_h, zo_h)

! ... calculate canopy radiation (sunlit and shade) profile

                        call canopy_fsun_clu( fafraczInt, lairef, cluref, cszref, fsun)

! ... calculate canopy leaf temperature (sun/shade) profile

                        call canopy_tleaf_lin(zk, hcmref, tmp2mref, fsun, &
                            tleaf_sun, tleaf_shade, tleaf_ave)

! ... calculate canopy Photosynthetic Photon Flux Density (PPFD) (sun/shade) profile

                        call canopy_ppfd_exp(zk, hcmref, dswrfref, lairef, fsun, &
                            ppfd_sun, ppfd_shade, ppfd_ave)

! ... user option to calculate in-canopy wind speeds at height z and midflame WAF

                        if (ifcanwind .or. ifcanwaf) then
                            if (rsl_opt .eq. 0) then
                                do k=1, modlays
                                    call canopy_wind_most(hcmref, zk(k), fafraczInt(k), ubzref, &
                                        z0ghc, cdrag, pai, hgtref, d_h, zo_h, &
                                        lambdars, canBOT(k), canTOP(k), canWIND(loc,k))
                                end do
                            else
                                write(*,*) 'wrong RSL_OPT namelist option = ', rsl_opt, 'only option = 0 right now'
                                call exit(2)
                            end if

! ... determine midflamepoint and flame height from user or FRP calculation
                            call canopy_flameh(flameh_opt, flameh_set, dx(loc), modres, &
                                frpref, frp_fac, hcmref, lu_opt, vtyperef, flameh_cal, &
                                midflamepoint, flameh(loc))
                            if (firetype .eq. 0) then !forest/sub-canopy firetype
                                if (flameh(loc) .gt. 0.0) then !flameh must be > 0
                                    if (flameh(loc) .le. hcmref) then !only calculate when flameh <= FCH
                                        call canopy_waf(hcmref, lambdars, hgtref, flameh(loc), &
                                            firetype, d_h, zo_h, canBOT(midflamepoint), &
                                            canTOP(midflamepoint), waf(loc))
                                    else
                                        write(*,*) 'warning...sub-canopy type fire, but flameh > FCH, setting WAF=1'
                                        waf(loc) = 1.0_rk
                                    end if
                                end if
                            else  !grass/crops, above-canopy firetype
                                if (flameh(loc) .gt. 0.0) then !flameh still must be > 0
                                    call canopy_waf(hcmref, lambdars, hgtref, flameh(loc), &
                                        firetype, d_h, zo_h, canBOT(midflamepoint), &
                                        canTOP(midflamepoint), waf(loc))
                                end if
                            end if
                            if (waf(loc) .gt. 1.0_rk) then !Final check of WAF > 1, must be <=1
                                waf(loc) = 1.0_rk
                            end if
                        end if

! ... user option to calculate in-canopy eddy diffusivities at height z
                        if (ifcaneddy) then
                            do k=1, modlays
                                call canopy_eddyx(hcmref, zk(k), ustref, molref, Kz(loc, k))
                            end do
                        end if

! ... user option to calculate in-canopy eddy photolysis attenuation at height z
                        if (ifcanphot) then
                            if (cszref .ge. 0.0_rk) then !only calculate if cell isn't dark
                                call canopy_phot(fafraczInt, &
                                    lairef, cluref, cszref, rjcf(loc, :))
                            end if
                        end if

!.......user option to calculate in-canopy leafage influence and assigning LAI as per timestep
                        if (leafage_opt .eq. 0) then
                            ! Initialize pastlai and currentlai based on current timestep
                            if (nn .eq. 1) then
                                currentlai = lairef
                                pastlai = currentlai
                            else
                                pastlai = currentlai
                                currentlai = lairef
                            end if

                            !!! Check if the lai_tstep is greater than time_intvl
                            if (lai_tstep .ge. time_intvl) then
                                tsteplai = lai_tstep/86400.0_rk
                            else
                                WRITE (*, *) "Error: Input LAI time step cannot be less than model time step...exiting!!!"
                                CALL EXIT(1)
                            endif
                        end if !leafage_opt = 0 end

!.......user option to calculate historical leaf temperature and PAR for past 24-hours and 240-hours rolling average per timestep
                        !Initialize
                        if (hist_opt .eq. 0) then      !Use instantaneous only
                            ppfd_sun24(loc,:)     = ppfd_sun
                            ppfd_shade24(loc,:)   = ppfd_shade
                            tleaf_sun24(loc,:)    = tleaf_sun
                            tleaf_shade24(loc,:)  = tleaf_shade
                            tleaf_ave24(loc,:)    = tleaf_ave
                            ppfd_sun240(loc,:)    = ppfd_sun
                            ppfd_shade240(loc,:)  = ppfd_shade
                            tleaf_sun240(loc,:)   = tleaf_sun
                            tleaf_shade240(loc,:) = tleaf_shade
                            tleaf_ave240(loc,:)   = tleaf_ave
                            !for AQ, temp, and wind stress factors in biogenics, use instantaneous only (caution!)
                            daily_maxt2m(loc)     = tmp2mref
                            daily_mint2m(loc)     = tmp2mref
                            daily_maxws10m(loc)   = ubzref
                        else if (hist_opt .eq. 1) then  !Try for historical average values
                            ! Calculate weights for running means of historic variables
                            ! DNEWFRAC and DOLDFRAC are the weights given to the current
                            ! and existing value, respectively, when updating running means
                            ! over the last X days. HNEWFRAC and HOLDFRAC are the same but
                            ! for the 24H means.
                            dnewfrac = time_intvl / ( tau_days * 24.0_rk * 3600.0_rk )
                            doldfrac = 1.0_rk - dnewfrac
                            hnewfrac = time_intvl / ( tau_hours * 3600.0_rk )
                            holdfrac = 1.0_rk - hnewfrac

                            if (nn .le. 24) then
                                !Track times for moving time window (running) average below
                                ppfd_sun24_tmp(nn,loc,:)     = ppfd_sun
                                ppfd_shade24_tmp(nn,loc,:)   = ppfd_shade
                                tleaf_sun24_tmp(nn,loc,:)    = tleaf_sun
                                tleaf_shade24_tmp(nn,loc,:)  = tleaf_shade
                                tleaf_ave24_tmp(nn,loc,:)    = tleaf_ave
                                ppfd_sun240_tmp(nn,loc,:)    = ppfd_sun
                                ppfd_shade240_tmp(nn,loc,:)  = ppfd_shade
                                tleaf_sun240_tmp(nn,loc,:)   = tleaf_sun
                                tleaf_shade240_tmp(nn,loc,:) = tleaf_shade
                                tleaf_ave240_tmp(nn,loc,:)   = tleaf_ave
                                tmp2mref_tmp(nn,loc)         = tmp2mref
                                ubzref_tmp(nn,loc)           = ubzref

                                !TODO:  Restart capability needed to get past leaf temp and PAR if avaialble
                                !For now, if <= 24 hours then use current instantaneous
                                ppfd_sun24(loc,:)     = ppfd_sun
                                ppfd_shade24(loc,:)   = ppfd_shade
                                tleaf_sun24(loc,:)    = tleaf_sun
                                tleaf_shade24(loc,:)  = tleaf_shade
                                tleaf_ave24(loc,:)    = tleaf_ave
                                ppfd_sun240(loc,:)    = ppfd_sun
                                ppfd_shade240(loc,:)  = ppfd_shade
                                tleaf_sun240(loc,:)   = tleaf_sun
                                tleaf_shade240(loc,:) = tleaf_shade
                                tleaf_ave240(loc,:)   = tleaf_ave
                                !for AQ, temp, and wind stress factors in biogenics, still use instantaneous only (caution!)
                                daily_maxt2m(loc)     = tmp2mref
                                daily_mint2m(loc)     = tmp2mref
                                daily_maxws10m(loc)   = ubzref

                            else if (nn .le. 240) then
                                !Drop nth value into end time slice
                                ppfd_sun24_tmp(25,loc,:)     = ppfd_sun
                                ppfd_shade24_tmp(25,loc,:)   = ppfd_shade
                                tleaf_sun24_tmp(25,loc,:)    = tleaf_sun
                                tleaf_shade24_tmp(25,loc,:)  = tleaf_shade
                                tleaf_ave24_tmp(25,loc,:)    = tleaf_ave
                                ppfd_sun240_tmp(nn,loc,:)    = ppfd_sun
                                ppfd_shade240_tmp(nn,loc,:)  = ppfd_shade
                                tleaf_sun240_tmp(nn,loc,:)   = tleaf_sun
                                tleaf_shade240_tmp(nn,loc,:) = tleaf_shade
                                tleaf_ave240_tmp(nn,loc,:)   = tleaf_ave
                                tmp2mref_tmp(25,loc)         = tmp2mref
                                ubzref_tmp(25,loc)           = ubzref

                                !Calculate 24 and 240 hour averages
                                ppfd_sun24(loc,:)     = sum(ppfd_sun24_tmp(1:24,loc,:),1)/24.0
                                ppfd_shade24(loc,:)   = sum(ppfd_shade24_tmp(1:24,loc,:),1)/24.0
                                tleaf_sun24(loc,:)    = sum(tleaf_sun24_tmp(1:24,loc,:),1)/24.0
                                tleaf_shade24(loc,:)  = sum(tleaf_shade24_tmp(1:24,loc,:),1)/24.0
                                tleaf_ave24(loc,:)    = sum(tleaf_ave24_tmp(1:24,loc,:),1)/24.0
                                !Use 24 hour restarts until 240 restarts become available
                                ppfd_sun240(loc,:)    = sum(ppfd_sun24_tmp(1:24,loc,:),1)/24.0
                                ppfd_shade240(loc,:)  = sum(ppfd_shade24_tmp(1:24,loc,:),1)/24.0
                                tleaf_sun240(loc,:)   = sum(tleaf_sun24_tmp(1:24,loc,:),1)/24.0
                                tleaf_shade240(loc,:) = sum(tleaf_shade24_tmp(1:24,loc,:),1)/24.0
                                tleaf_ave240(loc,:)   = sum(tleaf_ave24_tmp(1:24,loc,:),1)/24.0
                                !Updated rolling 24 hour (hourly, short term) and 240 hour (daily, long-term) averages
                                ppfd_sun24(loc,:)     = ( holdfrac * ppfd_sun24(loc,:) )     + ( hnewfrac * ppfd_sun )
                                ppfd_shade24(loc,:)   = ( holdfrac * ppfd_shade24(loc,:) )   + ( hnewfrac * ppfd_shade )
                                tleaf_sun24(loc,:)    = ( holdfrac * tleaf_sun24(loc,:) )    + ( hnewfrac * tleaf_sun )
                                tleaf_shade24(loc,:)  = ( holdfrac * tleaf_shade24(loc,:) )  + ( hnewfrac * tleaf_shade )
                                tleaf_ave24(loc,:)    = ( holdfrac * tleaf_ave24(loc,:) )    + ( hnewfrac * tleaf_ave )
                                ppfd_sun240(loc,:)     = ( holdfrac * ppfd_sun24(loc,:) )     + ( hnewfrac * ppfd_sun )
                                ppfd_shade240(loc,:)   = ( holdfrac * ppfd_shade24(loc,:) )   + ( hnewfrac * ppfd_shade )
                                tleaf_sun240(loc,:)    = ( holdfrac * tleaf_sun24(loc,:) )    + ( hnewfrac * tleaf_sun )
                                tleaf_shade240(loc,:)  = ( holdfrac * tleaf_shade24(loc,:) )  + ( hnewfrac * tleaf_shade )
                                tleaf_ave240(loc,:)    = ( holdfrac * tleaf_ave24(loc,:) )    + ( hnewfrac * tleaf_ave )
                                !Take daily max/min for temperature and daily max for wind speed
                                daily_maxt2m(loc)        = maxval(tmp2mref_tmp(1:24,loc))
                                daily_mint2m(loc)        = minval(tmp2mref_tmp(1:24,loc))
                                daily_maxws10m(loc)      = maxval(ubzref_tmp(1:24,loc))

                                !Shift time window to make room for new nth value
                                ppfd_sun24_tmp(1:24,loc,:)     = ppfd_sun24_tmp(2:25,loc,:)
                                ppfd_shade24_tmp(1:24,loc,:)   = ppfd_shade24_tmp(2:25,loc,:)
                                tleaf_sun24_tmp(1:24,loc,:)    = tleaf_sun24_tmp(2:25,loc,:)
                                tleaf_shade24_tmp(1:24,loc,:)  = tleaf_shade24_tmp(2:25,loc,:)
                                tleaf_ave24_tmp(1:24,loc,:)    = tleaf_ave24_tmp(2:25,loc,:)
                                tmp2mref_tmp(1:24,loc)         = tmp2mref_tmp(2:25,loc)
                                ubzref_tmp(1:24,loc)           = ubzref_tmp(2:25,loc)

                            else
                                !Drop nth value into end time slice
                                ppfd_sun24_tmp(25,loc,:)     = ppfd_sun
                                ppfd_shade24_tmp(25,loc,:)   = ppfd_shade
                                tleaf_sun24_tmp(25,loc,:)    = tleaf_sun
                                tleaf_shade24_tmp(25,loc,:)  = tleaf_shade
                                tleaf_ave24_tmp(25,loc,:)    = tleaf_ave
                                ppfd_sun240_tmp(241,loc,:)    = ppfd_sun
                                ppfd_shade240_tmp(241,loc,:)  = ppfd_shade
                                tleaf_sun240_tmp(241,loc,:)   = tleaf_sun
                                tleaf_shade240_tmp(241,loc,:) = tleaf_shade
                                tleaf_ave240_tmp(241,loc,:)   = tleaf_ave
                                tmp2mref_tmp(25,loc)         = tmp2mref
                                ubzref_tmp(25,loc)           = ubzref

                                !Calculate 24 and 240 hour averages
                                ppfd_sun24(loc,:)     = sum(ppfd_sun24_tmp(1:24,loc,:),1)/24.0
                                ppfd_shade24(loc,:)   = sum(ppfd_shade24_tmp(1:24,loc,:),1)/24.0
                                tleaf_sun24(loc,:)    = sum(tleaf_sun24_tmp(1:24,loc,:),1)/24.0
                                tleaf_shade24(loc,:)  = sum(tleaf_shade24_tmp(1:24,loc,:),1)/24.0
                                tleaf_ave24(loc,:)    = sum(tleaf_ave24_tmp(1:24,loc,:),1)/24.0
                                ppfd_sun240(loc,:)    = sum(ppfd_sun240_tmp(1:240,loc,:),1)/240.0
                                ppfd_shade240(loc,:)  = sum(ppfd_shade240_tmp(1:240,loc,:),1)/240.0
                                tleaf_sun240(loc,:)   = sum(tleaf_sun240_tmp(1:240,loc,:),1)/240.0
                                tleaf_shade240(loc,:) = sum(tleaf_shade240_tmp(1:240,loc,:),1)/240.0
                                tleaf_ave240(loc,:)   = sum(tleaf_ave240_tmp(1:240,loc,:),1)/240.0
                                !Updated rolling 24 hour (hourly, short term) and 240 hour (daily, long-term) averages
                                ppfd_sun24(loc,:)     = ( holdfrac * ppfd_sun24(loc,:) )     + ( hnewfrac * ppfd_sun )
                                ppfd_shade24(loc,:)   = ( holdfrac * ppfd_shade24(loc,:) )   + ( hnewfrac * ppfd_shade )
                                tleaf_sun24(loc,:)    = ( holdfrac * tleaf_sun24(loc,:) )    + ( hnewfrac * tleaf_sun )
                                tleaf_shade24(loc,:)  = ( holdfrac * tleaf_shade24(loc,:) )  + ( hnewfrac * tleaf_shade )
                                tleaf_ave24(loc,:)    = ( holdfrac * tleaf_ave24(loc,:) )    + ( hnewfrac * tleaf_ave )
                                ppfd_sun240(loc,:)    = ( doldfrac * ppfd_sun240(loc,:) )    + ( dnewfrac * ppfd_sun )
                                ppfd_shade240(loc,:)  = ( doldfrac * ppfd_shade240(loc,:) )  + ( dnewfrac * ppfd_shade )
                                tleaf_sun240(loc,:)   = ( doldfrac * tleaf_sun240(loc,:) )   + ( dnewfrac * tleaf_sun )
                                tleaf_shade240(loc,:) = ( doldfrac * tleaf_shade240(loc,:) ) + ( dnewfrac * tleaf_shade )
                                tleaf_ave240(loc,:)   = ( doldfrac * tleaf_ave240(loc,:) )   + ( dnewfrac * tleaf_ave )
                                !Take daily max/min for temperature and daily max for wind speed
                                daily_maxt2m(loc)        = maxval(tmp2mref_tmp(1:24,loc))
                                daily_mint2m(loc)        = minval(tmp2mref_tmp(1:24,loc))
                                daily_maxws10m(loc)      = maxval(ubzref_tmp(1:24,loc))

                                !Shift time window to make room for new nth value
                                ppfd_sun24_tmp(1:24,loc,:)     = ppfd_sun24_tmp(2:25,loc,:)
                                ppfd_shade24_tmp(1:24,loc,:)   = ppfd_shade24_tmp(2:25,loc,:)
                                tleaf_sun24_tmp(1:24,loc,:)    = tleaf_sun24_tmp(2:25,loc,:)
                                tleaf_shade24_tmp(1:24,loc,:)  = tleaf_shade24_tmp(2:25,loc,:)
                                tleaf_ave24_tmp(1:24,loc,:)    = tleaf_ave24_tmp(2:25,loc,:)
                                ppfd_sun240_tmp(1:240,loc,:)    = ppfd_sun240_tmp(2:241,loc,:)
                                ppfd_shade240_tmp(1:240,loc,:)  = ppfd_shade240_tmp(2:241,loc,:)
                                tleaf_sun240_tmp(1:240,loc,:)   = tleaf_sun240_tmp(2:241,loc,:)
                                tleaf_shade240_tmp(1:240,loc,:) = tleaf_shade240_tmp(2:241,loc,:)
                                tleaf_ave240_tmp(1:240,loc,:)   = tleaf_ave240_tmp(2:241,loc,:)
                                tmp2mref_tmp(1:24,loc)          = tmp2mref_tmp(2:25,loc)
                                ubzref_tmp(1:24,loc)            = ubzref_tmp(2:25,loc)
                            end if
                        else
                            write(*,*) 'wrong HIST_OPT namelist option = ', hist_opt, 'only option = 0 or 1'
                            call exit(2)
                        end if
! ... user option to calculate in-canopy biogenic emissions
                        if (ifcanbio) then
                            if (cszref .ge. 0.0_rk .and. dswrfref .gt. 0.0_rk &
                                .and. cluref .gt. 0.0_rk) then
                                !ISOP
                                if (biospec_opt == 0 .or. biospec_opt == 1) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 1, emi_isop(loc,:))
                                end if
                                !MYRC
                                if (biospec_opt == 0 .or. biospec_opt == 2) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 2, emi_myrc(loc,:))
                                end if
                                !SABI
                                if (biospec_opt == 0 .or. biospec_opt == 3) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 3, emi_sabi(loc,:))
                                end if
                                !LIMO
                                if (biospec_opt == 0 .or. biospec_opt == 4) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 4, emi_limo(loc,:))
                                end if
                                !CARE
                                if (biospec_opt == 0 .or. biospec_opt == 5) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 5, emi_care(loc,:))
                                end if
                                !OCIM
                                if (biospec_opt == 0 .or. biospec_opt == 6) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 6, emi_ocim(loc,:))
                                end if
                                !BPIN
                                if (biospec_opt == 0 .or. biospec_opt == 7) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 7, emi_bpin(loc,:))
                                end if
                                !APIN
                                if (biospec_opt == 0 .or. biospec_opt == 8) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 8, emi_apin(loc,:))
                                end if
                                !MONO
                                if (biospec_opt == 0 .or. biospec_opt == 9) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 9, emi_mono(loc,:))
                                end if
                                !FARN
                                if (biospec_opt == 0 .or. biospec_opt == 10) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 10, emi_farn(loc,:))
                                end if
                                !CARY
                                if (biospec_opt == 0 .or. biospec_opt == 11) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 11, emi_cary(loc,:))
                                end if
                                !SESQ
                                if (biospec_opt == 0 .or. biospec_opt == 12) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 12, emi_sesq(loc,:))
                                end if
                                !MBOL
                                if (biospec_opt == 0 .or. biospec_opt == 13) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 13, emi_mbol(loc,:))
                                end if
                                !METH
                                if (biospec_opt == 0 .or. biospec_opt == 14) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 14, emi_meth(loc,:))
                                end if
                                !ACET
                                if (biospec_opt == 0 .or. biospec_opt == 15) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 15, emi_acet(loc,:))
                                end if
                                !CO
                                if (biospec_opt == 0 .or. biospec_opt == 16) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 16, emi_co(loc,:))
                                end if
                                !BIDI VOC
                                if (biospec_opt == 0 .or. biospec_opt == 17) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 17, emi_bvoc(loc,:))
                                end if
                                !Stress VOC
                                if (biospec_opt == 0 .or. biospec_opt == 18) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 18, emi_svoc(loc,:))
                                end if
                                !Other VOC
                                if (biospec_opt == 0 .or. biospec_opt == 19) then
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 19, emi_ovoc(loc,:))
                                end if
                            end if
                        end if

! ... user option to calculate in-canopy dry deposition velocity
                        if (ifcanddepgas ) then
                            if (ifcanwind) then   !ubar needed for rbl
                                !if(ddepspec_opt == 0 .or. ddepspec_opt == 3) then TBD for all species =
                                call canopy_gas_drydep_zhang(zk, hcmref, tka(loc,:), pressa(loc,:), &
                                    relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canWIND(loc,:),  &
                                    dswrfref, 3, ddep_o3(loc,:))   ! [cm/s]
                                !endif
                                print*,ddep_o3(loc,:)
                            else
                                write(*,*)  'Wrong IfCanWind choice of ', ifcanwind, ' in namelist...exiting'
                                write(*,*)  'Set IfCanwind to True to use IfCanDDepGas'
                                call exit(2)
                            end if
                        end if
                    end if !Contiguous Canopy

                else
!                    write(*,*)  'Warning VIIRS/MODIS VTYPE ', vtyperef, ' is not supported...continue'
                end if   !Vegetation types
            else
                write(*,*)  'Wrong LU_OPT choice of ', lu_opt, ' in namelist...exiting'
                call exit(2)

            end if       !Landuse Options
        end do  !Lat*Lon

    else
        write(*,*)  'Wrong INFMT_OPT choice of ', infmt_opt, ' in namelist...exiting'
        call exit(2)

    end if !Input Format (1D or 2D text or NetCDF)

END SUBROUTINE canopy_calcs
