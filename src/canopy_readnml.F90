SUBROUTINE canopy_readnml

!-------------------------------------------------------------------------------
! Name:     Read Canopy Namelist
! Purpose:  Reads input namelist to get user control variables.
!           15 Jul 2022  Original Version (P.C. Campbell)
!
!-------------------------------------------------------------------------------
    USE canopy_files_mod
    USE canopy_canopts_mod
    USE canopy_coord_mod
    use canopy_canvars_mod, ONLY: zk

    IMPLICIT NONE

    !local variables
    INTEGER                            :: istat
    INTEGER                            :: n,i
    CHARACTER(LEN=*),      PARAMETER   :: pname = 'CANOPY_READNML'

    NAMELIST /filenames/ file_vars, file_canvars, file_out

    NAMELIST /userdefs/  infmt_opt, time_start, time_end, time_intvl, ntime, &
        nlat, nlon, modlays, modres, href_opt, href_set, z0ghc, lambdars, &
        var3d_opt, var3d_set, pavd_opt, pavd_set, &
        flameh_opt, flameh_cal, flameh_set, frp_fac, ifcanwind, &
        ifcanwaf, ifcaneddy, ifcanphot, ifcanbio, ifcanddepgas, pai_opt, pai_set, lu_opt, &
        z0_opt, dx_opt, dx_set, lai_thresh, cf_thresh, ch_thresh, rsl_opt, bio_cce, &
        biospec_opt, biovert_opt, ssg_opt, ssg_set, crop_opt, crop_set, co2_opt, co2_set, &
        leafage_opt, lai_tstep, soim_opt, soild1, soild2, soild3, soild4, aq_opt, w126_set, &
        ht_opt, lt_opt, hw_opt, hist_opt, loss_opt, loss_set, loss_ind, lifetime, &
        ddepspecgas_opt, chemmechgas_opt, chemmechgas_tot


!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

    CHARACTER(LEN=256), PARAMETER :: f9000 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR OPENING CANOPY NAMELIST FILE ON UNIT ', i3, &
    & /, 1x, '***   NAMELIST FILE NAME = ', a, &
    & /, 1x, '***   IOSTAT = ', i4, &
    & /, 1x, 70('*'))"


    CHARACTER(LEN=256), PARAMETER :: f9050 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR READING NAMELIST FILE ON UNIT ', i3, &
    & /, 1x, '***   NAMELIST FILE NAME = ', a, &
    & /, 1x, '***   NAMELIST = ', a, &
    & /, 1x, '***   IOSTAT = ', i4, &
    & /, 1x, 70('*'))"

!-------------------------------------------------------------------------------
! Open canopy namelist file.
!-------------------------------------------------------------------------------

    OPEN (iutnml, FILE=file_nml, STATUS='OLD', IOSTAT=istat)

    IF ( istat > 0 ) THEN
        WRITE (*,f9000) TRIM(pname), iutnml, TRIM(file_nml), istat
        CALL EXIT(istat)
    ENDIF

!-------------------------------------------------------------------------------
! Initialize canopy file names.
!-------------------------------------------------------------------------------

    file_vars(:)    = " "
    file_canvars(:) = " "
    file_out(:)     = " "

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for 1D (txt or ncf) or 2D (ncf only) input file format
! (default = 0, i.e., 2D)
    infmt_opt = 0
!-------------------------------------------------------------------------------


! Set default value for number of timesteps, start/end, and interval (seconds)
    ntime      =  0
    time_start = '0000-00-00-00:00:00.0000'
    time_end   = '0000-00-00-00:00:00.0000'
    time_intvl =  0
!-------------------------------------------------------------------------------

! Set default value for number of latitude/longitude cells (default = 1D point)
    nlat = 1
    nlon = 1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer value for number of canopy layers (default = 100 layers)
    modlays = 100
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for canopy vertical resolution (m) (Default = 0.5 m)
    modres = 0.5_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for using input 3D variables from file (default = 0)
    var3d_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default value for number of input 3D levels in variables from file (Default = 14)
    var3d_set = 14
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for using 3D GEDI PAVD inputs from file (default = 0)
    pavd_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for latitude threshold when using 3D GEDI PAVD inputs from file (default = 52 degrees)
    pavd_set = 52.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for reference height set values or array from file (default = 0)
    href_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default value for reference height above canopy (m) (Default = 10 m)
    href_set = 10.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for ratio of ground roughness length to canopy top height
    z0ghc = 0.0025_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for Influence function associated with roughness sublayer
    lambdars = 1.25_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for flame height set values or calculation (default = 0)
    flameh_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for FRP to flame height relationships used (default = 0)
    flameh_cal = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for flame height (m) (Default = 2.0 m)
    flameh_set = 2.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for FRP tuning factor for flameh (m) (Default = 1.0 )
    frp_fac = 1.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default logical for canopy wind option (default = .FALSE.)
    ifcanwind = .FALSE.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default logical for canopy WAF option (default = .FALSE.)
    ifcanwaf = .FALSE.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default logical for canopy eddy diffusivity (default = .FALSE.)
    ifcaneddy = .FALSE.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default logical for canopy photolysis attenuation (default = .FALSE.)
    ifcanphot = .FALSE.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default logical for canopy biogenic emissions (default = .FALSE.)
    ifcanbio = .FALSE.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default logical for canopy gas dry deposition (default = .FALSE.)
    ifcanddepgas = .FALSE.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for PAI set values or calculation (default = 0)
    pai_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for PAI set value (default = 4.0)
    pai_set = 4.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set integer for LU type used from model mapped to Massman et al. (default = 0)
    lu_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set integer for z0 estimate either from model or vegtype dependent (default = 0)
    z0_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for DX set values or calculation (default = 0)
    dx_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for DX Cell resolution set value (default = 1.0 m)
    dx_set = 1.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for user LAI threshold value for canopy (default = 0.1)
    lai_thresh = 0.1_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for user FRT threshold value for canopy (default = 0.5)
    cf_thresh = 0.5_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for user FCH threshold value for canopy (default = 0.5 m)
    ch_thresh = 0.1_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set integer for unified RSL option used in model (default = 0, off)
    rsl_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for MEGAN biogenic canopy environment coeficient (default = 0.21; Silva et al. (2020)
    bio_cce = 0.21_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer value to select species for biogenic emissions output (0, all)
    biospec_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer value for MEGAN vertical integration of emissions (0, off full leaf-level emissions)
    biovert_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for shrubs/savanaa/grasslands vegtype option from GEDI or user (default = 0)
    ssg_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default value for shrubs/savanaa/grasslands vegtype heights used in model (m) (Default = 1 m)
    ssg_set = 1.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for crop vegtype option from GEDI or user (default = 0)
    crop_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default value for crop vegtype heights used in model (m) (Default = 3 m)
    crop_set = 3.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for co2 inhibition option for biogenic isoprene emissions (default = 0; Possell & Hewitt (2011))
    co2_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default value for atmospheric co2 concentration for co2_opt (Default = 400.0 ppmv)
    co2_set = 400.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for Leaf Age response option for biogenic (all) emissions
! (default is OFF i.e., leafage_opt=1 making GAMMA_LEAFAGE=1 i.e. leaf age response to Biogenic VOCs is off.)
!  Otherwise switched ON i.e., leafage_opt= 0 for which lai_tstep is defined to  enable GAMMA_LEAFAGE calculation
    leafage_opt = 1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default timestep for LAI input as daily = 24*3600 seconds, otherwise specified in namelist
    lai_tstep = 86400 !Daily LAI inputs
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for turning on canopy loss ratios calculation or constant value for
! adjusting top of canopy net emissions (default = 0; Off)
    loss_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for constant loss factor (used only when loss_opt=2)
    loss_set = 0.96_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default value for above-canopy BVOC lifetime (s) used only with loss_opt=1 (Default = 3600 s)
    lifetime = 3600.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for applying loss factor to all species (default = 0) or specific biogenic indices only (> 0 )
    loss_ind = 0

! Set default integer for using historically averaged leaf temp and PAR for biogenic emissions (default=0; Off)
    hist_opt = 0
!-------------------------------------------------------------------------------

! Set default integer for soil moisture response option for biogenic (isoprene only) emissions
! (default is OFF i.e., leafage_opt=1 making GAMMA_SOIM=1 i.e. soil moisture response to Biogenic Isoprene VOCs is OFF.)
!  Otherwise switched ON i.e., soim_opt= 0 to  enable GAMMA_SOIM calculation
    soim_opt = 1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default value for depth of soil layer 1 (Default = 5.0 cm)
    soild1 = 5.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default value for depth of soil layer 2 (Default = 25.0 cm)
    soild2 = 25.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default value for depth of soil layer 3 (Default = 70.0 cm)
    soild3 = 70.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default value for depth of soil layer 4 (Default = 150.0 cm)
    soild4 = 150.0_rk
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Set default integer for using air quality stress gamma for biogenic emissions (default=2; Off)
    aq_opt = 2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default value for spatially constant w126 value (ppm-hours)
    w126_set = 20.0_rk
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Set default integer for using high temperature stress gamma for biogenic emissions (default=1; Off)
    ht_opt = 1
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Set default integer for using low temperature stress gamma for biogenic emissions (default=1; Off)
    lt_opt = 1
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Set default integer for using high wind stress gamma for biogenic emissions (default=1; Off)
    hw_opt = 1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer value to select species for drydep gas output (0, all)
    ddepspecgas_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer value to select chemical mechanism (0, RACM2)
    chemmechgas_opt = 0
!-------------------------------------------------------------------------------

! Set default integer value to select chemical mechanism gas species list including transported (31, RACM2)
    chemmechgas_tot = 31
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Read namelist to get user definitions.  Rewind namelist file after each
! read in case namelists are not in the correct order in the namelist.
!-------------------------------------------------------------------------------

    READ (iutnml, filenames, IOSTAT=istat)
    IF ( istat > 0 ) THEN
        WRITE (*,f9050) TRIM(pname), iutnml, TRIM(file_nml), "filenames", istat
        CALL EXIT(istat)
    ENDIF
    REWIND (iutnml)

    READ (iutnml, userdefs, IOSTAT=istat)
    IF ( istat > 0 ) THEN
        WRITE (*,f9050) TRIM(pname), iutnml, TRIM(file_nml), "userdefs", istat
        CALL EXIT(istat)
    ENDIF
    REWIND (iutnml)

!-------------------------------------------------------------------------------
! Crop blank spaces off ends of file names.
!-------------------------------------------------------------------------------

    DO n = 1, SIZE(file_vars)
        file_vars(n)= TRIM( ADJUSTL( file_vars(n) ) )
    ENDDO

    DO n = 1, SIZE(file_canvars)
        file_canvars(n)= TRIM( ADJUSTL( file_canvars(n) ) )
    ENDDO

    DO n = 1, SIZE(file_out)
        file_out(n)= TRIM( ADJUSTL( file_out(n) ) )
    ENDDO

!-------------------------------------------------------------------------------
! Verify values of user-defined options (need conditions added...)
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Derive canopy model profile heights from user NL input
!-------------------------------------------------------------------------------

    if(.not.allocated(zk))         allocate(zk(modlays))
    zk(1) = 0.0
    do i=2, modlays
        zk(i)   = zk(i-1) + modres
    end do

!-------------------------------------------------------------------------------
! Close namelist file.
!-------------------------------------------------------------------------------
!  write(*,*)'namelist intvl=',intvl
    CLOSE (iutnml)

END SUBROUTINE canopy_readnml
