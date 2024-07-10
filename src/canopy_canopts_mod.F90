MODULE canopy_canopts_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Option Variable Descriptions
! Purpose:  Contains canopy option variable descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!-------------------------------------------------------------------------------
    use canopy_const_mod, ONLY: rk
    IMPLICIT NONE

!! .... defines canopy options (read from user namelist)
    integer             ::    infmt_opt   !Integer for choosing 1D or 2D input file format (default = 0, 2D)
    integer             ::    var3d_opt   !Integer for choosing if 3D variables will be read from file (default = 0, off)
    integer             ::    var3d_set   !Integer for number fo 3D levels in input file, if var3d_opt=1 (default = 14)
    integer             ::    pavd_opt    !Integer for turning on 3D GEDI PAVD profiles read from input file (default = 0, off)
    real(rk)            ::    pavd_set    !Real value for latitude +/- threshold when 3D GEDI PAVD profiles read from input file (default = 52.0)
    integer             ::    href_opt    !Integer for using set href in namelist (=0) or array from file(=1) (default = 0)
    real(rk)            ::    href_set    !Set reference Height above canopy @ 10 m  (m)
    logical             ::    ifcanwind   !logical canopy wind option (default = .FALSE.)
    logical             ::    ifcanwaf    !logical canopy WAF option (default = .FALSE.)
    logical             ::    ifcaneddy   !logical canopy eddy Kz option (default = .FALSE.)
    logical             ::    ifcanphot   !logical canopy photolsyis atten option (default = .FALSE.)
    logical             ::    ifcanbio    !logical canopy biogenic emissions option (default = .FALSE.)
    integer             ::    pai_opt     !integer for PAI values used or calculated (default = 0)
    real(rk)            ::    pai_set     !real value for PAI set values used (default = 4.0)
    integer             ::    lu_opt      !integer for LU type from model mapped to Massman et al. (default = 0/VIIRS)
    integer             ::    z0_opt      !integer for setting first estimate of z0 (default = 0 for Z0_MOD)
    integer             ::    flameh_opt  !Integer for flameh values used or calculated (default = 0)
    integer             ::    flameh_cal  !Integer for FRP to flame height relationships used (default = 0)
    real(rk)            ::    flameh_set  !User Set Flame Height (m)
    real(rk)            ::    frp_fac     !FRP tuning factor for flame height calculation (default = 1.0)
    integer             ::    dx_opt      !Integer for dx resolution values used or calculated (default = 0)
    real(rk)            ::    dx_set      !User Set Grid Cell Resolution (m)
    real(rk)            ::    lai_thresh  !User set grid cell LAI threshold to apply canopy conditions (m2/m2)
    real(rk)            ::    cf_thresh   !User set grid cell canopy fraction threshold to apply canopy conditions ()
    real(rk)            ::    ch_thresh   !User set grid cell canopy height threshold to apply canopy conditions (m)
    integer             ::    rsl_opt     !RSL option used in model from Rosenzweig et al. 2021 (default = 0, off)
    real(rk)            ::    z0ghc       !ratio of ground roughness length to canopy top height
    real(rk)            ::    lambdars    !Value representing influence of roughness sublayer (nondimensional)
    real(rk)            ::    bio_cce     !MEGAN biogenic emission canopy environment coefficient.
    integer             ::    biospec_opt !Set default integer for species output option (default = 0, all)
    integer             ::    biovert_opt !MEGAN vertical integration of emissions option (default = 0, off)
    integer             ::    ssg_opt     !Set default integer for shrubs/savanna/grassland vegtype option from GEDI or user (default = 0)
    real(rk)            ::    ssg_set     !Set default value for shrubs/savanna/grassland vegtype heights used in model (m) (Default = 1 m)
    integer             ::    crop_opt    !Set default integer for crop vegtype option from GEDI or user (default = 0)
    real(rk)            ::    crop_set    !Set default value for crop vegtype heights used in model (m) (Default = 3 m)
    integer             ::    co2_opt     !Set default integer for co2 inhibition option for biogenic isoprene emissions (default=0; Possell & Hewitt (2011))
    real(rk)            ::    co2_set     !Set default value for atmospheric co2 concentration for co2_opt (m) (Default = 400.0 ppmv)
    integer             ::    leafage_opt !Set default = 1 for Leaf Age factor option for BVOCs, =1 give GAMMA_LEAFGAE =1
    integer             ::    lai_tstep   !Set default = 24*3600 seconds (Daily timestep for LAI input, = 30*24*3600 for say Monthly time ste or other user-defined in seconds)
    integer             ::    loss_opt    !Set default integer for turning on canopy loss ratios for adjusting top of canopy net emissions (default = 0; Off)
    real(rk)            ::    lifetime    !Set default value for above-canopy BVOC lifetime (s) used with loss_opt=1 (Default = 3600 s)
    real(rk)            ::    loss_set    !Set default value for constant canopy loss factor applied used with loss_opt=2 (Default = 0.96)
    integer             ::    loss_ind    !Set default integer for applying canopy loss factor to all species (=0) or only specific biogenics specie indices (> 0)
    integer             ::    hist_opt    !Set default integer for using historically averaged leaf temp and PAR for biogenic emissions (default=0; Off)
    integer             ::    soim_opt    !Set default = 1 for Soil Moisture factor option for BVOCs, =1 give GAMMA_SOIM =1
    real(rk)            ::    soild1      !User set real value of depth of soil layer 1 (default = 5 cm @ centerpoint, based on Noah/Noah-MP)
    real(rk)            ::    soild2      !User set real value of depth of soil layer 2 (default = 25 cm @ centerpoint, based on Noah/Noah-MP)
    real(rk)            ::    soild3      !User set real value of depth of soil layer 3 (default = 70 cm @ centerpoint, based on Noah/Noah-MP)
    real(rk)            ::    soild4      !User set real value of depth of soil layer 4 (default = 150 cm @ centerpoint, based on Noah/Noah-MP)
    integer             ::    aq_opt      !Set default integer for air quality stress index for gamma_aq in biogenic emissions
    real(rk)            ::    w126_set    !Set default value for constant ozone W126 value (ppm-hours)
    integer             ::    ht_opt      !Set default integer for high temperature stress index for gamma_ht in biogenic emissions
    integer             ::    lt_opt      !Set default integer for low temperature stress index for gamma_lt in biogenic emissions
    integer             ::    hw_opt      !Set default integer for high wind stress index for gamma_hw in biogenic emissions
END MODULE canopy_canopts_mod
