module canopy_utils_mod

!    use canopy_const_mod, ONLY: ntotal, pi, rk, rearth    !constants for canopy models
    use canopy_const_mod, ONLY: pi, rk, rearth

    implicit none

    private
    public IntegrateTrapezoid,interp_linear1_internal,CalcPAI, &
        CalcDX,CalcFlameH,GET_GAMMA_CO2,GET_GAMMA_LEAFAGE, &
        GET_GAMMA_SOIM,GET_GAMMA_AQ,GET_GAMMA_HT,GET_GAMMA_LT, &
        GET_GAMMA_HW,GET_CANLOSS_BIO,CalcTemp,CalcPressure,esat, &
        CalcRelHum,CalcSpecHum,CalcCair,Convert_qh_to_h2o, &
        SetMolecDiffSTP,MolecDiff,rs_zhang_gas,rbl,rcl,rml, &
        SetEffHenrysLawCoeffs,SetReactivityParams,ReactivityParam, &
        EffHenrysLawCoeff

contains

    function IntegrateTrapezoid(x, y)
        !! Calculates the integral of an array y with respect to x using the trapezoid
        !! approximation. Note that the mesh spacing of x does not have to be uniform.
        real(rk), intent(in)  :: x(:)                !! Variable x
        real(rk), intent(in)  :: y(size(x))          !! Function y(x)
        real(rk)              :: IntegrateTrapezoid  !! Integral ∫y(x)·dx
        ! Integrate using the trapezoidal rule
        associate(n => size(x))
            IntegrateTrapezoid = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
        end associate
    end function
!-------------------------------------------------------------------------------------

    function interp_linear1_internal(x,y,xout) result(yout)
        !! Interpolates for the y value at the desired x value,
        !! given x and y values around the desired point.

        implicit none

        real(rk), intent(IN)  :: x(2), y(2), xout
        real(rk) :: yout
        real(rk) :: alph

        if ( xout .lt. x(1) .or. xout .gt. x(2) ) then
            write(*,*) "interp1: xout < x0 or xout > x1 !"
            write(*,*) "xout = ",xout
            write(*,*) "x0   = ",x(1)
            write(*,*) "x1   = ",x(2)
            stop
        end if

        alph = (xout - x(1)) / (x(2) - x(1))
        yout = y(1) + alph*(y(2) - y(1))

        return

    end function interp_linear1_internal
!--------------------------------------------------------------------------------------

    function CalcPAI(ch, canfrac)
        !! Calculates the Plant Area Index as a function of canopy height and canopy/
        !! forest fraction (Based on Eq. 19 of Massman et al., 2017).

        !!  W.J. Massman, J.M. Forthofer, and M.A. Finney. An improved
        !!  canopy wind model for predicting wind adjustment factors
        !!  and wildland fire behavior. Canadian Journal of Forest Research.
        !!  47(5): 594-603. https://doi.org/10.1139/cjfr-2016-0354

        !! Assume Canopy Cover Fraction, C,  = CANFRAC
        !! Assume Canopy Crown Ratio, F, = CC/3.0 = CANFRAC/3.0 (Eq. 9 in  Andrews, 2012).
        !! Andrews, P.L. 2012. Modeling wind adjustment factor and midflame wind speed
        !!  for Rothermel’s surface fire spread model. USDA For. Serv. Gen. Tech. Rep. RMRS-GTR-266.


        real(rk), intent(in)  :: ch           !! Input Grid cell canopy height (m
        real(rk), intent(in)  :: canfrac      !! Input Grid cell canopy fraction
        real(rk)              :: CalcPAI      !! Calculated Plant area index (PAI)

        CalcPAI=( (ch*(canfrac/3.0_rk)*10.6955_rk) / (2.0_rk * pi) ) * canfrac !Massman PAI calculation (Eq. 19)

    end function
    !--------------------------------------------------------------------------------------

    real(rk) function CalcDX(lat, dlon) result(dx)
        !! Compute the zonal distance, dx, corresponding to longitude increment `dlon`.

        real(rk), intent(in) :: lat   !! Latitude (degrees)
        real(rk), intent(in) :: dlon  !! Longitude increment (degrees)
        real(rk) :: lat_rad, dlon_rad

        lat_rad = lat * pi / 180._rk
        dlon_rad = dlon * pi / 180._rk

        dx = rearth * cos(lat_rad) * dlon_rad

    end function
    !--------------------------------------------------------------------------------------

    ! real(rk) function CalcGCDist(lat1, lat2, lon1, lon2) result(d)
    !     !! Compute great-circle distance between two points using the spherical law of cosines formula.

    !     real(rk), intent(in)  :: lat1,lat2                  !! Two model latitudes
    !     real(rk), intent(in)  :: lon1,lon2                  !! Two model longitudes
    !     real(rk) :: lat_rad1, lat_rad2, lon_rad1, lon_rad2  !! radians

    !     lat_rad1 = lat1/(180.0_rk/pi)
    !     lon_rad1 = lon1/(180.0_rk/pi)
    !     lat_rad2 = lat2/(180.0_rk/pi)
    !     lon_rad2 = lon2/(180.0_rk/pi)

    !     d = rearth*acos( &
    !         sin(lat_rad1)*sin(lat_rad2) &
    !         + cos(lat_rad1)*cos(lat_rad2)*cos(lon_rad2-lon_rad1) &
    !     )

    ! end function
    ! !--------------------------------------------------------------------------------------

    real(rk) function CalcFlameH(frp, dx, lu_opt, vtype, flameh_cal) result( CalcFlameH_set )
        !! Approximates the Flame Height as a function of FRP intensity, grid cell distance (dx),
        !! and vegetation type (Based on Alexander and Cruz 2012).

        !!  Alexander, Martin E.; Cruz, Miguel G. 2012. Interdependencies between flame length and
        !!  fireline intensity in predicting crown fire initiation and crown scorch height.
        !!  International Journal of Wildland Fire 21(2):95-113.

        !! Assume Flame Length = Flame Height under calm winds

        real(rk), intent(in)  :: frp              !! Input Grid cell Fire Radiative Power (MW/cell)
        real(rk), intent(in)  :: dx               !! Input Grid cell length (m)
        integer,  intent(in)  :: lu_opt           !! Supported land use classifications
        integer,  intent(in)  :: vtype            !! Grid cell dominant vegetation type
        integer,  intent(in)  :: flameh_cal       !! Option of vegtype dependent FRP to Flame Height relationships used
        real(rk)              :: CalcFlameH_B59   !! Fireline intensity to flame length (=height) (m): Byram (1959)
        real(rk)              :: CalcFlameH_A66_L !! Fireline intensity to flame length (=height) (m): Anderson et al. (1966)
        real(rk)              :: CalcFlameH_A66_D !! Fireline intensity to flame length (=height) (m): Anderson et al. (1966)
        real(rk)              :: CalcFlameH_N80   !! Fireline intensity to flame length (=height) (m): Nelson (1980)
        real(rk)              :: CalcFlameH_C83_H !! Fireline intensity to flame length (=height) (m): Clark (1983)
        real(rk)              :: CalcFlameH_C83_B !! Fireline intensity to flame length (=height) (m): Clark (1983)
        real(rk)              :: CalcFlameH_N86   !! Fireline intensity to flame length (=height) (m): Nelson and Adkins (1986)
        real(rk)              :: CalcFlameH_V86   !! Fireline intensity to flame length (=height) (m): van Wilgen (1986)
        real(rk)              :: CalcFlameH_B94   !! Fireline intensity to flame length (=height) (m): Burrows (1994)
        real(rk)              :: CalcFlameH_W96   !! Fireline intensity to flame length (=height) (m): Weise and Biging (1996)
        real(rk)              :: CalcFlameH_V98   !! Fireline intensity to flame length (=height) (m): Vega et al. (1998)
        real(rk)              :: CalcFlameH_C98   !! Fireline intensity to flame length (=height) (m): Catchpole et al. (1998)
        real(rk)              :: CalcFlameH_F00   !! Fireline intensity to flame length (=height) (m): Fernandes et al. (2000)
        real(rk)              :: CalcFlameH_B04   !! Fireline intensity to flame length (=height) (m): Butler et al. (2004)
        real(rk)              :: CalcFlameH_F09_H !! Fireline intensity to flame length (=height) (m): Fernandes et al. (2009)
        real(rk)              :: CalcFlameH_F09_B !! Fireline intensity to flame length (=height) (m): Fernandes et al. (2009)
        real(rk)              :: CalcFlameS_M71   !! Fireline intensity to flame scorch (m): McArthur (1971)
        real(rk)              :: CalcFlameS_V73   !! Fireline intensity to flame scorch (m): Van Wagner (1973)
        real(rk)              :: CalcFlameS_C78   !! Fireline intensity to flame scorch (m): Cheyney (1978)
        real(rk)              :: CalcFlameS_L78   !! Fireline intensity to flame scorch (m): Luke and McArthur (1978)
        real(rk)              :: CalcFlameS_B88   !! Fireline intensity to flame scorch (m): Burrows et al. (1988)
        real(rk)              :: CalcFlameS_B89   !! Fireline intensity to flame scorch (m): Burrows et al. (1989)
        real(rk)              :: CalcFlameS_S90   !! Fireline intensity to flame scorch (m): Saveland et al. (1990)
        real(rk)              :: CalcFlameS_F93   !! Fireline intensity to flame scorch (m): Finney and Martin (1993)
        real(rk)              :: CalcFlameS_B94_1 !! Fireline intensity to flame scorch (m): Burrows (1994)
        real(rk)              :: CalcFlameS_B94_2 !! Fireline intensity to flame scorch (m): Burrows (1994)
        real(rk)              :: CalcFlameS_W98   !! Fireline intensity to flame scorch (m): Williams et al. (1998)
        real(rk)              :: CalcFlameS_F02   !! Fireline intensity to flame scorch (m): Fernandes (2002)
        real(rk)              :: CalcFlameS_set   !! Vegetation dependent average crown scorch height (m)


        if (flameh_cal .eq. 0) then !We assume that flame height=flame length,
            !and use the relationships in Table 1 of Alexander and Cruz (2012)

            !!  -------------------------------------------------------------------------------------
            !!  Byram, GM (1959). Combustion of Forest Fuels. In Forest Fire: Control and Use.
            !!  (Ed. KP David) pp. 61-89.  McGraw Hill, New York, NY
            !Pine litter with grass understorey (Evergreen and Grasses)
            CalcFlameH_B59=0.0775_rk*((frp*1000.0_rk)/dx)**0.46_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Anderson HE, Brackebusch AP, Mutch RW, Rothermel RC (1966) Mechanisms of fire spread
            !!  research progress report 2. USDA Forest Service,
            !!  Intermountain Forest and Range Experiment Station, Research Paper
            !!  INT-28. (Ogden, UT)
            !Lodgepole pine slash (Evergreen)
            CalcFlameH_A66_L=0.074_rk*((frp*1000.0_rk)/dx)**0.651_rk
            !Douglas-fir slash  (Evergreen)
            CalcFlameH_A66_D=0.0447_rk*((frp*1000.0_rk)/dx)**0.67_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Nelson RM Jr (1980) Flame characteristics for fires in southern fuels. USDA
            !!  Forest Service, Southeastern Forest Experiment Station, Research Paper
            !!  SE-205. (Asheville, NC)
            !Southern USA Fuels (Mixed Forest)
            CalcFlameH_N80=0.0377_rk*((frp*1000.0_rk)/dx)**0.50_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Clark RG (1983) Threshold requirements for fire spread in grassland fuels.
            !!  PhD dissertation, Texas Tech University, Lubbock.
            !Grasslands (head fire) (Grasslands)
            CalcFlameH_C83_H=0.00015_rk*((frp*1000.0_rk)/dx)**1.75_rk
            !Grasslands (backfire) (Grasslands)
            CalcFlameH_C83_B=0.000722_rk*((frp*1000.0_rk)/dx)**0.99_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Nelson RM Jr, Adkins CW (1986) Flame characteristics of wind-driven
            !!  surface fires. Canadian Journal of Forest Research 16, 1293–1300.
            !!  doi:10.1139/X86-229
            !Litter and shrubs (Shrublands)
            CalcFlameH_N86=0.0475_rk*((frp*1000.0_rk)/dx)**0.493_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  van Wilgen BW (1986) A simple relationship for estimating the intensity of
            !!  fires in natural vegetation. South African Journal of Botany 52, 384–385
            !Fynbos shrublands (Shrublands)
            CalcFlameH_V86=0.0075_rk*((frp*1000.0_rk)/dx)**0.46_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Burrows ND (1994) Experimental development of a fire management model
            !!  for jarrah (Eucalyptus marginata Donn ex Sm.) forest. PhD thesis,
            !!  Australian National University, Canberra.
            !Eucalypt Forest (Evergreens)
            CalcFlameH_B94=0.0147_rk*((frp*1000.0_rk)/dx)**0.767_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Weise DR, Biging GS (1996) Effects of wind velocity and slope on flame
            !!  properties. Canadian Journal of Forest Research 26, 1849–1858.
            !!  doi:10.1139/X26-210
            !Excelsior (Deciduous)
            CalcFlameH_W96=0.016_rk*((frp*1000.0_rk)/dx)**0.7_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Vega JA, Cuinas P, Fonturbel T, Perez-Gorostiaga P, Fernandez C (1998)
            !!  Predicting fire behaviour in Galician (NW Spain) shrubland fuel
            !!  complexes. In ‘Proceedings of 3rd International Conference on Forest
            !!  Flame length and fireline intensity interdependences.
            !Shrublands
            CalcFlameH_V98=0.087_rk*((frp*1000.0_rk)/dx)**0.493_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  CatchpoleWR, Bradstock RA, Choate J, Fogarty LG, Gellie N, McCarthy G,
            !!  McCaw WL, Marsden-Smedley JB, Pearce G (1998) Cooperative
            !!  development of equations for heathland fire behaviour. In ‘Proceedings
            !!  of 3rd International Conference on Forest Fire Research and 14th
            !!  Conference on Fire and Forest Meteorology, Volume II’, 16–20
            !!  November 1998, Luso–Coimbra, Portugal. (Ed. DX Viegas)
            !!  pp. 631–645. (University of Coimbra: Coimbra, Portugal)
            !Shrublands
            CalcFlameH_C98=0.0325_rk*((frp*1000.0_rk)/dx)**0.56_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Fernandes PM, Catchpole WR, Rego FC (2000) Shrubland fire behaviour
            !!  modelling with microplot data.Canadian Journal of Forest Research 30,
            !!  889–899. doi:10.1139/X00-012
            !Shrublands
            CalcFlameH_F00=0.0516_rk*((frp*1000.0_rk)/dx)**0.453_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Butler BW, Finney MA, Andrews PL, Albini FA (2004) A radiation-driven
            !!  model of crown fire spread. Canadian Journal of Forest Research 34,
            !!  1588–1599. doi:10.1139/X04-074
            !Jack Pine Forest - Crown Fire (Evergreen)
            CalcFlameH_B04=0.0175_rk*((frp*1000.0_rk)/dx)**0.667_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Fernandes PM, Botelho HS, Rego FC, Loureiro C (2009) Empirical
            !!  modelling of surface fire behaviour in maritime pine stands.
            !!  International Journal of Wildland Fire 18, 698–710. doi:10.1071/
            !!  WF08023
            !Maritime Pine (head fire) (Evergreen)
            CalcFlameH_F09_H=0.045_rk*((frp*1000.0_rk)/dx)**0.543_rk
            !Maritime Pine (backfire) (Evergreen)
            CalcFlameH_F09_B=0.029_rk*((frp*1000.0_rk)/dx)**0.724_rk
            !!  -------------------------------------------------------------------------------------

            if (lu_opt .eq. 0 .or. lu_opt .eq. 1) then !VIIRS or MODIS LU types
                if (vtype .ge. 1 .and. vtype .le. 2) then !VIIRS/MODIS Cat 1-2/Evergreen Needleleaf & Broadleaf
                    CalcFlameH_set=(CalcFlameH_B59+CalcFlameH_A66_L+CalcFlameH_A66_D+CalcFlameH_B94+ &
                        CalcFlameH_F09_H+CalcFlameH_F09_B) / 6.0_rk
                    if (((frp*1000.0_rk)/dx) .ge. 1700.0_rk ) then !Evergreen forest crowning likely
                        CalcFlameH_set=CalcFlameH_B04!for completeness in Alexander and Cruz, will leave here
                        !but is overidden in WAF calculation anyway for same crowning
                        !check and flame heigh set to canopy height.
                    end if
                else if (vtype .ge. 3 .and. vtype .le. 4) then !VIIRS/MODIS Cat 3-4 Deciduous Needleleaf and  Broadleaf
                    CalcFlameH_set=(CalcFlameH_N80+CalcFlameH_W96) / 2.0_rk
                else if (vtype .eq. 5) then !VIIRS/MODIS Cat 5 Mixed Forests
                    CalcFlameH_set=CalcFlameH_N80
                else if (vtype .ge. 6 .and. vtype .le. 7) then !VIIRS/MODIS Cat 6-7 Shrublands
                    CalcFlameH_set=(CalcFlameH_N86+CalcFlameH_V86+CalcFlameH_V98+CalcFlameH_C98+ &
                        CalcFlameH_F00) / 5.0_rk
                else if (vtype .ge. 8 .and. vtype .le. 10) then !VIIRS/MODIS Cat 8-10 Savannas and Grasslands
                    CalcFlameH_set=(CalcFlameH_B59+CalcFlameH_C83_H+CalcFlameH_C83_B) / 3.0_rk
                else if (vtype .ge. 12 .and. vtype .lt. 13) then !VIIRS/MODIS Cat 12 Croplands
                    CalcFlameH_set=(CalcFlameH_B59+CalcFlameH_C83_H+CalcFlameH_C83_B) / 3.0_rk
                else if (vtype .ge. 14 .and. vtype .lt. 15) then !VIIRS/MODIS Cat 14 Cropland/Natural Mosaic
                    CalcFlameH_set=(CalcFlameH_B59+CalcFlameH_C83_H+CalcFlameH_C83_B) / 3.0_rk
                else if (vtype .ge. 18 .and. vtype .le. 19) then !VIIRS/MODIS Cat 18-19 Wooded and Mixed Tundra
                    CalcFlameH_set=(CalcFlameH_N86+CalcFlameH_V86+CalcFlameH_V98+CalcFlameH_C98+ &
                        CalcFlameH_F00) / 5.0_rk
                else
                    CalcFlameH_set=(CalcFlameH_B59+CalcFlameH_A66_L+CalcFlameH_A66_D+CalcFlameH_N80+ &
                        CalcFlameH_C83_H+CalcFlameH_C83_B+CalcFlameH_N86+CalcFlameH_V86+ &
                        CalcFlameH_B94+CalcFlameH_W96+CalcFlameH_V98+CalcFlameH_C98+ &
                        CalcFlameH_F00+CalcFlameH_B04+CalcFlameH_F09_H+ &
                        CalcFlameH_F09_B) / 16.0_rk
                end if
            else
                write(*,*)  'Wrong LU_OPT choice of ', LU_OPT, 'in namelist, only VIIRS/MODIS available right now...exiting'
                call exit(2)
            end if



        else if (flameh_cal .eq. 1) then!We use the crown scorch height relationships in
            !Table 2 of Alexander and Cruz (2012) and Equation 14
            !that directly relates flame height to crown scorch height


            !!  -------------------------------------------------------------------------------------
            !!  McArthur AG (1971) Aspects of fire control in the P. caribaea and
            !!  P. elliottii plantations of north-western Viti Levu, Fiji Islands.
            !!  Commonwealth of Australia, Forest and Timber Bureau, Forest
            !!  Research Institute. (Canberra, ACT)
            !Slash and Carribean Pine (Evergreen)
            CalcFlameS_M71=0.1226_rk*((frp*1000.0_rk)/dx)**0.667_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Van Wagner CE (1973) Height of crown scorch in forest fires. Canadian
            !!  Journal of Forest Research 3, 373–378. doi:10.1139/X73-055
            !Red Pine, White Pine, Jack Pine, and Northern Red Oak (Mixed Forest)
            CalcFlameS_V73=0.1483_rk*((frp*1000.0_rk)/dx)**0.667_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Cheney NP (1978) Guidelines for fire management on forested watersheds,
            !!  based on Australian experience. In ‘Special Readings in Conservation’,
            !!  FAO Conservation Guide 4, pp. 1–37. (Food and Agriculture Organization
            !!  of the United Nations: Rome, Italy)
            !Eucalypt Forest (Evergreen)
            CalcFlameS_C78=0.1297_rk*((frp*1000.0_rk)/dx)**0.667_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Luke RH, McArthur AG (1978) ‘Bushfires in Australia.’ (Australian
            !!  Government Publishing Service: Canberra, ACT)
            !Eucalypt Forest (Evergreen)
            CalcFlameS_L78=0.1523_rk*((frp*1000.0_rk)/dx)**0.667_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Burrows ND, Smith RH, RobinsonAD(1988) Prescribed burning slash fuels
            !!  in Pinus radiata plantations in Western Australia. Western Australia
            !!  Department of Conservation and Land Management, Technical Report
            !!  20. (Perth, WA)
            !Radiata Pine Thinning Slash (Evergreen)
            CalcFlameS_B88=0.1579_rk*((frp*1000.0_rk)/dx)**0.667_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Burrows ND, Woods YC, Ward BG, Robinson AD (1989) Prescribing low
            !!  intensity fire to kill wildings in Pinus radiata plantations in Western
            !!  Australia. Australian Forestry 52, 45–52.
            !Radiata Pine Wildings (Evergreen)
            if (frp .gt. 0) then
                CalcFlameS_B89=max(0.0_rk,(0.248_rk*((frp*1000.0_rk)/dx)**0.667_rk) - 0.41)
            else
                CalcFlameS_B89=0.0_rk
            end if
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Saveland JM, Bakken SR, Neuenschwander LF (1990) Predicting mortality
            !!  and scorch height from prescribed burning for ponderosa pine in
            !!  northern Idaho. University of Idaho, College of Forestry, Wildlife and
            !!  Range Sciences, Idaho Forest, Wildlife and Range Experiment Station,
            !!  Station Bulletin 53. (Moscow, ID)
            !Ponderosa Pine (Evergreen)
            CalcFlameS_S90=0.063_rk*((frp*1000.0_rk)/dx)**0.667_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Finney MA, Martin RE (1993) Modeling effects of prescribed fire on younggrowth
            !!  coast redwood trees. Canadian Journal of Forest Research 23,
            !!  1125–1135. doi:10.1139/X93-143
            !Coast Redwood (Evergreen)
            CalcFlameS_F93=0.228_rk*((frp*1000.0_rk)/dx)**0.667_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  BurrowsND(1994) Experimental development of a fire management model
            !!  for jarrah (Eucalyptus marginata Donn ex Sm.) forest. PhD thesis,
            !!  Australian National University, Canberra.
            !Jarrah Forest - spring (Evergreen)
            CalcFlameS_B94_1=0.28_rk*((frp*1000.0_rk)/dx)**0.58_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  BurrowsND(1994) Experimental development of a fire management model
            !!  for jarrah (Eucalyptus marginata Donn ex Sm.) forest. PhD thesis,
            !!  Australian National University, Canberra.
            !Jarrah Forest - summer (Evergreen)
            CalcFlameS_B94_2=0.36_rk*((frp*1000.0_rk)/dx)**0.59_rk
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Williams RJ, Gill AM, Moore PHR (1998) Seasonal changes in fire
            !!  behaviour in a tropical savanna in Northern Australia. International
            !!  Journal of Wildland Fire 8, 227–239. doi:10.1071/WF9980227
            !Grassland-eucalypt savanna (Savanna, Grassland, Crops)
            if (frp .gt. 0) then
                CalcFlameS_W98=(21.2_rk - 17.6_rk)*exp(0.000287_rk*((frp*1000.0_rk)/dx))
            else
                CalcFlameS_W98=0.0_rk
            end if
            !!  -------------------------------------------------------------------------------------

            !!  -------------------------------------------------------------------------------------
            !!  Fernandes PM (2002) Desenvolvimento de relacoes predictivas para uso no
            !!  planeamento de fogo controlado em povoamentos de Pinus pinaster Ait.
            !!  [Development of predictive relationships for use in planning prescribed
            !!  fire in Pinus pinaster Ait. stands]. PhD thesis, Universidade de Tras os
            !!  Montes e Alto Douro, Vilas Real, Portugal. [In Portugese]
            !Maritime Pine Head Fire (Evergreen)
            CalcFlameS_F02=0.125_rk*((frp*1000.0_rk)/dx)**0.724_rk
            !!  -------------------------------------------------------------------------------------

            if (lu_opt .eq. 0 .or. lu_opt .eq. 1) then !VIIRS or MODIS LU types
                if (vtype .ge. 1 .and. vtype .le. 2) then !VIIRS/MODIS Cat 1-2/Evergreen Needleleaf & Broadleaf
                    CalcFlameS_set=(CalcFlameS_M71+CalcFlameS_C78+CalcFlameS_L78+CalcFlameS_B88+ &
                        CalcFlameS_B89+CalcFlameS_S90+CalcFlameS_F93+CalcFlameS_B94_1+ &
                        CalcFlameS_B94_2+CalcFlameS_F02) / 10.0_rk
                    CalcFlameH_set=(CalcFlameS_set/5.232_rk)**(1.0_rk/0.756_rk)
                else if (vtype .ge. 3 .and. vtype .le. 4) then !VIIRS/MODIS Cat 3-4 Deciduous Needleleaf and  Broadleaf
                    CalcFlameS_set=CalcFlameS_V73
                    CalcFlameH_set=(CalcFlameS_set/5.232_rk)**(1.0_rk/0.756_rk)
                else if (vtype .eq. 5) then !VIIRS/MODIS Cat 5 Mixed Forests
                    CalcFlameS_set=CalcFlameS_V73
                    CalcFlameH_set=(CalcFlameS_set/5.232_rk)**(1.0_rk/0.756_rk)
                else if (vtype .ge. 6 .and. vtype .le. 7) then !VIIRS/MODIS Cat 6-7 Shrublands
                    CalcFlameS_set=CalcFlameS_W98
                    CalcFlameH_set=(CalcFlameS_set/5.232_rk)**(1.0_rk/0.756_rk)
                else if (vtype .ge. 8 .and. vtype .le. 10) then !VIIRS/MODIS Cat 8-10 Savannas and Grasslands
                    CalcFlameS_set=CalcFlameS_W98
                    CalcFlameH_set=(CalcFlameS_set/5.232_rk)**(1.0_rk/0.756_rk)
                else if (vtype .ge. 12 .and. vtype .lt. 13) then !VIIRS/MODIS Cat 12 Croplands
                    CalcFlameS_set=CalcFlameS_W98
                    CalcFlameH_set=(CalcFlameS_set/5.232_rk)**(1.0_rk/0.756_rk)
                else if (vtype .ge. 14 .and. vtype .lt. 15) then !VIIRS/MODIS Cat 14 Cropland/Natural Mosaic
                    CalcFlameS_set=CalcFlameS_W98
                    CalcFlameH_set=(CalcFlameS_set/5.232_rk)**(1.0_rk/0.756_rk)
                else if (vtype .ge. 18 .and. vtype .le. 19) then !VIIRS/MODIS Cat 18-19 Wooded and Mixed Tundra
                    CalcFlameS_set=CalcFlameS_W98
                    CalcFlameH_set=(CalcFlameS_set/5.232_rk)**(1.0_rk/0.756_rk)
                else
                    CalcFlameS_set=(CalcFlameS_M71+CalcFlameS_C78+CalcFlameS_L78+CalcFlameS_B88+ &
                        CalcFlameS_B89+CalcFlameS_S90+CalcFlameS_F93+CalcFlameS_B94_1+ &
                        CalcFlameS_B94_2+CalcFlameS_F02+CalcFlameS_V73+CalcFlameS_W98) / 12.0_rk
                    CalcFlameH_set=(CalcFlameS_set/5.232_rk)**(1.0_rk/0.756_rk)
                end if
            else
                write(*,*)  'Wrong LU_OPT choice of ', LU_OPT, 'in namelist, only VIIRS/MODIS available right now...exiting'
                call exit(2)
            end if
        else
            write(*,*)  'Wrong FLAMEH_CAL choice of ', FLAMEH_CAL, 'in namelist...exiting'
            call exit(2)
        end if

    end function

    real(rk) function GET_GAMMA_CO2(co2_opt, co2_set)  result( GAMMA_CO2 )
        ! !IROUTINE: get_gamma_co2
        !
        ! !DESCRIPTION: Function GET\_GAMMA\_CO2 computes the CO2 activity factor
        !  associated with CO2 inhibition of isoprene emission. Called from
        !  GET\_MEGAN\_EMISSIONS only.
        !\\
        !\\
        ! !INTERFACE:

        ! !INPUT PARAMETERS:
        integer,  INTENT(IN) :: co2_opt       ! Option for co2 inhibition calculation
        ! 0=Possell & Hewitt (2011);
        ! 1=Wilkinson et al. (2009)
        ! >1=off
        real(rk), INTENT(IN) :: co2_set       ! User set atmospheric CO2 conc [ppmv]
        !
        ! !RETURN VALUE:
!        REAL(rk)             :: GAMMA_CO2  ! CO2 activity factor [unitless]
        !
        ! !LOCAL VARIABLES:
        REAL(rk)             :: CO2i       ! Intercellular CO2 conc [ppmv]
        REAL(rk)             :: ISMAXi     ! Asymptote for intercellular CO2
        REAL(rk)             :: HEXPi      ! Exponent for intercellular CO2
        REAL(rk)             :: CSTARi     ! Scaling coef for intercellular CO2
        REAL(rk)             :: ISMAXa     ! Asymptote for atmospheric CO2
        REAL(rk)             :: HEXPa      ! Exponent for atmospheric CO2
        REAL(rk)             :: CSTARa     ! Scaling coef for atmospheric CO2
        !
        ! !REMARKS:
        !  References:
        !  ============================================================================
        !  (1 ) Heald, C. L., Wilkinson, M. J., Monson, R. K., Alo, C. A.,
        !       Wang, G. L., and Guenther, A.: Response of isoprene emission
        !       to ambient co(2) changes and implications for global budgets,
        !       Global Change Biology, 15, 1127-1140, 2009.
        !  (2 ) Wilkinson, M. J., Monson, R. K., Trahan, N., Lee, S., Brown, E.,
        !       Jackson, R. B., Polley, H. W., Fay, P. A., and Fall, R.: Leaf
        !       isoprene emission rate as a function of atmospheric CO2
        !       concentration, Global Change Biology, 15, 1189-1200, 2009.
        !  (3 ) Possell, M., and Hewitt, C. N.: Isoprene emissions from plants
        !       are mediated by atmospheric co2 concentrations, Global Change
        !       Biology, 17, 1595-1610, 2011.

        ! !REVISION HISTORY:
        !  (1 ) Implemented in the standard code by A. Tai (Jun 2012).
        !  See https://github.com/geoschem/hemco for complete history
        !EOP
        !------------------------------------------------------------------------------
        !BOC

        !-----------------------
        ! Compute GAMMA_CO2
        !-----------------------

        !----------------------------------------------------------
        ! Choose between two alternative CO2 inhibition schemes
        !----------------------------------------------------------

        IF ( co2_opt .eq. 0 ) THEN

            ! Use empirical relationship of Possell & Hewitt (2011):
            ! Empirical relationship of Possell & Hewitt (2011) based on nine
            ! experimental studies including Wilkinson et al. (2009).

            GAMMA_CO2 = 8.9406_rk / ( 1.0_rk + 8.9406_rk * 0.0024_rk * co2_set )


        ELSEIF ( co2_opt .eq. 1 ) THEN

            ! Use parameterization of Wilkinson et al. (2009):
            ! Semi-process-based parameterization of Wilkinson et al. (2009),
            ! taking into account of sensitivity to intercellular CO2
            ! fluctuation, which is here set as a constant fraction of
            ! atmospheric CO2. This is especially recommended for sub-ambient
            ! CO2 concentrations::


            ! Parameters for intercellular CO2 using linear interpolation:
            IF ( co2_set <= 600.0_rk ) THEN
                ISMAXi = 1.036_rk  - (1.036_rk - 1.072_rk) / &
                    (600.0_rk - 400.0_rk) * (600.0_rk - co2_set)
                HEXPi  = 2.0125_rk - (2.0125_rk - 1.7000_rk) / &
                    (600.0_rk - 400.0_rk) * (600.0_rk - co2_set)
                CSTARi = 1150.0_rk - (1150.0_rk - 1218.0_rk) / &
                    (600.0_rk - 400.0_rk) * (600.0_rk - co2_set)
            ELSEIF ( co2_set > 600.0_rk .AND. co2_set < 800.0_rk ) THEN
                ISMAXi = 1.046_rk  - (1.046_rk - 1.036_rk) / &
                    (800.0_rk - 600.0_rk) * (800.0_rk - co2_set)
                HEXPi  = 1.5380_rk - (1.5380_rk - 2.0125_rk) / &
                    (800.0_rk - 600.0_rk) * (800.0_rk - co2_set)
                CSTARi = 2025.0_rk - (2025.0_rk - 1150.0_rk) / &
                    (800.0_rk - 600.0_rk) * (800.0_rk - co2_set)
            ELSE
                ISMAXi = 1.014_rk - (1.014_rk - 1.046_rk) / &
                    (1200.0_rk - 800.0_rk) * (1200.0_rk - co2_set)
                HEXPi  = 2.8610_rk - (2.8610_rk - 1.5380_rk) / &
                    (1200.0_rk - 800.0_rk) * (1200.0_rk - co2_set)
                CSTARi = 1525.0_rk - (1525.0_rk - 2025.0_rk) / &
                    (1200.0_rk - 800.0_rk) * (1200.0_rk - co2_set)
            ENDIF

            ! Parameters for atmospheric CO2:
            ISMAXa    = 1.344_rk
            HEXPa     = 1.4614_rk
            CSTARa    = 585.0_rk

            ! For now, set CO2_Ci = 0.7d0 * CO2_Ca as recommended by Heald
            ! et al. (2009):
            CO2i      = 0.7_rk * co2_set

            ! Compute GAMMA_CO2:
            GAMMA_CO2 = ( ISMAXi -  ISMAXi * CO2i**HEXPi / &
                ( CSTARi**HEXPi + CO2i**HEXPi ) )  &
                * ( ISMAXa - ISMAXa * ( 0.7_rk * co2_set )**HEXPa / &
                ( CSTARa**HEXPa + ( 0.7_rk * co2_set )**HEXPa ) )

        ELSE

            ! No CO2 inhibition scheme is used; GAMMA_CO2 set to unity:
            GAMMA_CO2 = 1.0_rk

        ENDIF

    end function GET_GAMMA_CO2

    function GET_GAMMA_LEAFAGE(leafage_opt,LAIpast,LAIcurrent,tsteplai,TABOVE,Anew,Agro,Amat,Aold)  result( GAMMA_LEAFAGE )
        ! ROUTINE: GET_GAMMA_LEAFAGE
        !
        ! !DESCRIPTION: Function GET_GAMMA_LEAFAGE computes the leaf age activity factor
        !  associated with foliage fraction calculation.
        !
        !     leaf age response to Biogenic VOCs
        ! Revision: Sept 2023 Quazi Z. Rasool NOAA CSL/CIRES
        !----------------------------------------------------------------
        !
        !       GAMLA = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
        !       where Fnew = new foliage fraction
        !             Fgro = growing foliage fraction
        !             Fmat = mature foliage fraction
        !             Fold = old foliage fraction
        !             Anew = emission activity for new foliage
        !             Agro = emission activity for growing foliage
        !             Amat = emission activity for mature foliage
        !             Aold = emission activity for old foliage
        !           "Age class fractions are determined from LAI changes"
        !             LAIcurrent = current Month's LAI (asuuming monthly LAI but can be
        !             customized as per tsteplai)
        !             LAIpast = past Month's LAI
        !             tsteplai  = length of the time step (days) i.e. days in between LAIpast and LAIcurrent
        !             ti = days between budbreak and emission induction (calculated below)
        !             tm = days between budbreak and peak emission (Calculated below)
        !             TABOVE = 2-meter temperature (K) TEMP2 or tmp2mref from the input model/obs
        !------------------------------------------------------------------------------
        ! !INTERFACE:

        ! !INPUT PARAMETERS:
        INTEGER,  INTENT(IN) :: leafage_opt       ! Option for leaf age emission factor calculation
        ! 0=On;
        ! 1 or  >1 =off i.e. GAMMA_LEAFAGE =1

        REAl(rk), INTENT(IN) :: tsteplai    ! time step, number of days between Past and Current LAI inputs
        REAL(rk), INTENT(IN) :: TABOVE      ! Above canopy temperature (K), t2m or tmpsfc
        REAL(rk), INTENT(IN) :: LAIpast     ! Past LAI [cm2/cm2]
        REAL(rk), INTENT(IN) :: LAIcurrent  ! Current LAI [cm2/cm2]
        REAL(rk), INTENT(IN) :: Anew        ! Relative emiss factor (new leaves)
        REAL(rk), INTENT(IN) :: Agro        ! Relative emiss factor (growing leaves)
        REAL(rk), INTENT(IN) :: Amat        ! Relative emiss factor (mature leaves)
        REAL(rk), INTENT(IN) :: Aold        ! Relative emiss factor (old leaves)
        !
        ! !RETURN VALUE:
        REAL(rk) :: GAMMA_LEAFAGE           ! Leaf age activity factor [unitless]
        !
        ! !LOCAL VARIABLES:
        REAL(rk) :: LAIdelta = 1.0e-10_rk   ! Tolerance for LAI comparison
        !INTEGER :: tsteplai                ! time step
        REAL(rk) :: Fnew, Fgro              ! foliage fractions
        REAL(rk) :: Fmat, Fold
        REAL(rk) :: ti, tm

        !
        ! !REMARKS:
        !  The function computes the leaf age activity factor based on BVOC's leaf age response.
        !
        ! !REVISION HISTORY:
        !  Adapted from HEMCO: https://github.com/geoschem/hemco for complete history
        !EOP
        !------------------------------------------------------------------------------
        !BOC

        !TABOVE (Tt in MEGAN code) -> Above Canopy TEMP (tmp2mref or temp2)

        ! Calculate foliage fraction
        !-----------------------
        !Also, Compute ti and tm
        ! ti: number of days after budbreak required to induce emissions
        ! tm: number of days after budbreak required to reach peak emissions

        !-----------------------
        ! Compute GAMMA_AGE
        !-----------------------
        IF (leafage_opt .eq. 0) THEN
            IF (LAIcurrent - LAIpast > LAIdelta) THEN !(i.e. LAI has Increased)
                IF (TABOVE .le. 303.0_rk) THEN
                    ti = 5.0_rk + 0.7_rk*(300.0_rk-TABOVE)
                ELSE
                    ti = 2.9_rk
                ENDIF
                tm = 2.3_rk*ti
                !Fnew calculated
                IF (ti .ge. tsteplai) THEN
                    Fnew = 1.0_rk - (LAIpast/LAIcurrent)
                ELSE
                    Fnew = (ti/tsteplai) * ( 1.0_rk-(LAIpast/LAIcurrent) )
                ENDIF
                !Fmat calculated
                IF (tm .ge. tsteplai) THEN
                    Fmat = LAIpast/LAIcurrent
                ELSE
                    Fmat = (LAIpast/LAIcurrent) + ( (tsteplai-tm)/tsteplai ) * ( 1.0_rk-(LAIpast/LAIcurrent) )
                ENDIF

                Fgro = 1.0_rk - Fnew - Fmat
                Fold = 0.0_rk
            ELSEIF (LAIpast - LAIcurrent > LAIdelta) THEN !(i.e. LAI has decreased)
                Fnew = 0.0_rk
                Fgro = 0.0_rk
                Fold = ( LAIpast-LAIcurrent ) / LAIpast
                Fmat = 1.0_rk-Fold
            ELSE !(LAIpast == LAIcurrent) THEN !If LAI remains same
                Fnew = 0.0_rk
                Fgro = 0.1_rk
                Fmat = 0.8_rk
                Fold = 0.1_rk
            ENDIF
            ! Compute GAMMA_LEAFAGE
            GAMMA_LEAFAGE = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
            ! Prevent negative values
            GAMMA_LEAFAGE = MAX( GAMMA_LEAFAGE , 0.0_rk )
        ELSEIF (leafage_opt .eq. 1) THEN
            GAMMA_LEAFAGE = 1.0_rk
        ELSE
            GAMMA_LEAFAGE = 1.0_rk
        ENDIF

    end function GET_GAMMA_LEAFAGE

    function GET_CANLOSS_BIO(loss_opt,lifetime,ustar,ch)  result( CANLOSS_BIO )
        ! ROUTINE: CANLOSS_BIO
        !
        ! !DESCRIPTION: Function to calculate BVOC canopy loss ratio due to approximate chemistry and dep loss
        !               Useful for comparing individual BVOC primary emissions to above canopy flux measurements
        !
        ! Based on Guenther et al. (2006) www.atmos-chem-phys.net/6/3181/2006/
        ! Note:  Formulation and emprical parameters based on isoprene oxidation chemistry Only -- Exercise caution applying to other BVOCs
        !
        ! Revision: Feb 2024 Patrick C. Campbell GMU/NOAA-ARL
        !----------------------------------------------------------------
        !
        !       LOSS RATIO = 1 -  D/(lambda*ustar*tau+D)
        !       where D = Canopy Depth (Assumed 1/3 the canopy height) (m)
        !             lambda = 0.3 (Assumed and usually PFT dependent)
        !             ustar = above canopy friction velocity (m/s)
        !             tau = above canopy chemical lifetime (default = 3600 s, isoprene) (s)
        !------------------------------------------------------------------------------
        ! !INTERFACE:

        ! !INPUT PARAMETERS:
        INTEGER,  INTENT(IN) :: loss_opt       ! Option for canopy loss function
        ! 0 = Calculation Off; CANLOSS_BIO = 1
        ! 1 = Calculation On i.e. CANLOSS_BIO < 1

        REAl(rk), INTENT(IN) :: lifetime    ! Above canopy BVOC chemical lifetime [s]
        REAL(rk), INTENT(IN) :: ustar       ! Above canopy friction velocity [m/s]
        REAL(rk), INTENT(IN) :: ch          ! Canopy Height [m]
        !
        ! !RETURN VALUE:
        REAL(rk) :: CANLOSS_BIO             ! Canopy loss factor [unitless]
        !
        ! !LOCAL VARIABLES:
        REAL(rk) :: lambda   = 0.3_rk         ! Assumed value for empirical parameter (Guenther et al., 2006)

        IF (loss_opt .eq. 0) THEN ! Off
            CANLOSS_BIO = 1.0_rk
        ELSE IF (loss_opt .eq. 1) THEN
            CANLOSS_BIO = 1.0_rk - ((ch*0.3333_rk)/(lambda*ustar*lifetime+(ch*0.3333_rk)))
        ELSE
            write(*,*)  'Wrong LOSS_OPT choice of ', loss_opt, 'in namelist...exiting'
            call exit(2)
        END IF

    end function GET_CANLOSS_BIO

    function GET_GAMMA_SOIM(soim_opt,soim1,soim2,soim3,soim4, &
        soid1,soid2,soid3,soid4,wilt,roota,rootb)  result( GAMMA_SOIM )
        ! ROUTINE: GET_GAMMA_SOIM
        !
        ! !DESCRIPTION: Function GET_GAMMA_SOIM computes the soil moisture activity factor
        !  associated with volumetric soil moisture and wilting point
        !
        !     Soil moisture response to Biogenic VOC Emissions
        ! Revision: February 08, 2024:  Patrick C. Campbell
        !----------------------------------------------------------------
        ! Based on Guenther et al. (2006),
        ! https://acp.copernicus.org/articles/6/3181/2006/
        ! And uses Zeng (2001) for PFT dependent root depth fractions,
        ! https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
        !------------------------------------------------------------------------------
!        ! !INPUT PARAMETERS:
        INTEGER,  INTENT(IN) :: soim_opt     ! Option for soil moisture emission factor calculation
!        ! 0=On;
!        ! 1 or  >1 =off i.e. GAMMA_SOIM =1

        REAL(rk), INTENT(IN)  :: soim1     ! Volumetric soil moisture layer 1 [m3/m3]
        REAL(rk), INTENT(IN)  :: soim2     ! Volumetric soil moisture layer 2 [m3/m3]
        REAL(rk), INTENT(IN)  :: soim3     ! Volumetric soil moisture layer 3 [m3/m3]
        REAL(rk), INTENT(IN)  :: soim4     ! Volumetric soil moisture layer 4 [m3/m3]
        REAL(rk), INTENT(IN)  :: soid1     ! Soil depth layer 1 [cm]
        REAL(rk), INTENT(IN)  :: soid2     ! Soil depth layer 2 [cm]
        REAL(rk), INTENT(IN)  :: soid3     ! Soil depth layer 3 [cm]
        REAL(rk), INTENT(IN)  :: soid4     ! Soil depth layer 4 [cm]
        REAL(rk), INTENT(IN)  :: wilt      ! Wilting point [proportion]
        REAL(rk), INTENT(IN)  :: roota     ! Coefficient A for PFT dependent cumulative root depth fraction [m-1]
        REAL(rk), INTENT(IN)  :: rootb     ! Coefficient B for PFT dependent cumulative root depth fraction [m-1]

        ! !RETURN VALUE:
        REAL(rk) :: GAMMA_SOIM            ! Soil moisture activity factor [unitless]

!        ! !LOCAL VARIABLES:
        REAL(rk) :: delta_theta1 = 0.06_rk                  ! Empirical parameter based on Pegoraro et al. (2004)
        ! https://www.publish.csiro.au/fp/FP04142
        REAL(rk) :: theta1                                  ! wilt + delta_theta1

        REAL(rk) :: gamma_soim1,gamma_soim2,gamma_soim3,gamma_soim4 ! gamma soil values for each soim1-4 [unitless]
        REAL(rk) :: rootfrac1,rootfrac2,rootfrac3,rootfrac4 !cumulative root fraction from surface to depth [fraction]
        REAL(rk) :: rootfrac_sum                            !sum of cumulative root fractions from each layer
!        !
!        ! !REMARKS:
!        !  The function computes the soil moisture activity factor
!        !
!        !EOP
!        !------------------------------------------------------------------------------
!        !BOC
!
!        !-----------------------
!        ! Compute GAMMA_SOIM for each layer and take weighted average across fraction of roots
!        !-----------------------

!Calculte the gamma soil values for each soim1-4
        IF (soim_opt .eq. 0) THEN
            theta1 = wilt + delta_theta1
            !Soil Layer 1
            if (soim1 .ge. theta1) then
                gamma_soim1 = 1.0_rk
            else if (soim1 .lt. theta1 .and. soim1 .gt. wilt ) then
                gamma_soim1 = (soim1 - wilt)/delta_theta1
            else
                gamma_soim1 = 0.0_rk
            end if
            !Soil Layer 2
            if (soim2 .ge. theta1) then
                gamma_soim2 = 1.0_rk
            else if (soim2 .lt. theta1 .and. soim2 .gt. wilt ) then
                gamma_soim2 = (soim2 - wilt)/delta_theta1
            else
                gamma_soim2 = 0.0_rk
            end if
            !Soil Layer 3
            if (soim3 .ge. theta1) then
                gamma_soim3 = 1.0_rk
            else if (soim3 .lt. theta1 .and. soim3 .gt. wilt ) then
                gamma_soim3 = (soim3 - wilt)/delta_theta1
            else
                gamma_soim3 = 0.0_rk
            end if
            !Soil Layer 4
            if (soim4 .ge. theta1) then
                gamma_soim4 = 1.0_rk
            else if (soim4 .lt. theta1 .and. soim4 .gt. wilt ) then
                gamma_soim4 = (soim4 - wilt)/delta_theta1
            else
                gamma_soim4 = 0.0_rk
            end if
!Calculate the cumulative root fraction from surface to depth (Y), soild1-4
            rootfrac1 = 1.0_rk - 0.5_rk*(exp(-1.0_rk*roota*(soid1/100.0_rk)) + &
                exp(-1.0_rk*rootb*(soid1/100.0_rk)))
            rootfrac2 = 1.0_rk - 0.5_rk*(exp(-1.0_rk*roota*(soid2/100.0_rk)) + &
                exp(-1.0_rk*rootb*(soid2/100.0_rk)))
            rootfrac3 = 1.0_rk - 0.5_rk*(exp(-1.0_rk*roota*(soid3/100.0_rk)) + &
                exp(-1.0_rk*rootb*(soid3/100.0_rk)))
            rootfrac4 = 1.0_rk - 0.5_rk*(exp(-1.0_rk*roota*(soid4/100.0_rk)) + &
                exp(-1.0_rk*rootb*(soid4/100.0_rk)))
            rootfrac_sum = rootfrac1+rootfrac2+rootfrac3+rootfrac4
!Calculate the weighted average of gamma_soim based on the PFT-dependent cumulative root fraction for each soil layer
            GAMMA_SOIM = ((rootfrac1*gamma_soim1) + (rootfrac2*gamma_soim2) + (rootfrac3*gamma_soim3) + &
                (rootfrac4*gamma_soim4))/rootfrac_sum
        ELSE IF (soim_opt .eq. 1) THEN
            GAMMA_SOIM = 1.0_rk
        ELSE
            GAMMA_SOIM = 1.0_rk
        END IF

    end function GET_GAMMA_SOIM

!----------------------------------------------------------------
!
!   Function GAMMA_AQ
!   EA response to air quality
!
!----------------------------------------------------------------

    real(rk) function GET_GAMMA_AQ(aq_opt, w126_ozone, w126_set, caq, taq, dtaq)  result( GAMMA_AQ )

        ! !IROUTINE: get_gamma_aq
        !
        ! !DESCRIPTION: Function GET\_GAMMA\_AQ computes the  activity factor
        !  associated with air qualtity (ozone W126 stress of biogenic emission. Called from
        !  GET\_MEGAN\_EMISSIONS only.
        !\\
        !\\
        ! !INTERFACE:

        ! !INPUT PARAMETERS:
        integer,  INTENT(IN) :: aq_opt       ! Option for aq stress calculation
        ! 0=MEGANv3.2 implementation with climatological, spatially dependent GFSv16-based W126;
        ! 1=MEGANv3.2 implementation with user-set, spatially constant value of ozone W126
        ! >1=off (gamma_aq=1)
        real(rk), INTENT(IN) :: w126_ozone     ! spatially dependent GFS w126 ozone [ppm-hours]
        real(rk), INTENT(IN) :: w126_set       ! user set spatially constant w126 ozone [ppm-hours]
        real(rk), INTENT(IN) :: taq            ! threshold for poor Air Quality stress (ppm-hours)
        real(rk), INTENT(IN) :: dtaq           ! delta threshold for poor Air Quality stress (ppm-hours)
        real(rk), INTENT(IN) :: caq            ! coefficient for poor Air Quality stress
        !
        ! !RETURN VALUE:
        !REAL(rk)             :: GAMMA_AQ  ! AQ activity factor [unitless]
        !
        ! !LOCAL VARIABLES:
        REAL(rk)             :: w126       ! local w126 value (ppm-hours)
        REAL(rk)             :: t1         ! combined threshold + delta threshold AQ value

        if (aq_opt .eq. 0) then !Use spatial GFS W126 ozone
            w126 = w126_ozone
        else                    !Use user set value of W126 ozone
            w126 = w126_set
        end if


        if (aq_opt .le. 1) then !Calculate GAMMA_AQ
            t1 = taq + dtaq
            if (w126 <= taq) then
                GAMMA_AQ = 1.0_rk
            else if ( w126 > taq .and. w126 < t1) then
                GAMMA_AQ = 1.0_rk + (caq - 1.0_rk)* (w126 - &
                    taq)/dtaq
            else
                GAMMA_AQ = caq
            end if
        else                    ! GAMMA_AQ = 1
            GAMMA_AQ = 1.0_rk
        end if

    end function GET_GAMMA_AQ
!-----------------------------------------------------------------------

    real(rk) function GET_GAMMA_HT(ht_opt, maxt2m, cht, tht, dtht)  result( GAMMA_HT )

        ! !IROUTINE: get_gamma_ht
        !
        ! !DESCRIPTION: Function GET\_GAMMA\_HT computes the  activity factor
        !  associated with high temperature stress of biogenic emission. Called from
        !  GET\_MEGAN\_EMISSIONS only.
        !\\
        !\\
        ! !INTERFACE:

        ! !INPUT PARAMETERS:
        integer,  INTENT(IN) :: ht_opt       ! Option for ht stress calculation
        ! 0=MEGANv3.2 On;
        ! 1=MEGANv3.2 Off
        real(rk), INTENT(IN) :: maxt2m         ! maximum input 2-m temperature [K]
        real(rk), INTENT(IN) :: tht            ! threshold for high temperature [K]
        real(rk), INTENT(IN) :: dtht           ! delta threshold for high temperature [K]
        real(rk), INTENT(IN) :: cht            ! coefficient for high temperature stress
        !
        ! !RETURN VALUE:
        !REAL(rk)             :: GAMMA_HT  ! HT activity factor [unitless]
        !
        ! !LOCAL VARIABLES:
        REAL(rk)             :: t1         ! combined threshold + delta threshold HT value

        if (ht_opt .eq. 0) then !Calculate GAMMA_HT
            t1 = tht + dtht
            if (maxt2m <= tht) then
                GAMMA_HT = 1.0_rk
            else if ( maxt2m > tht .and. maxt2m < t1) then
                GAMMA_HT = 1.0_rk + (cht - 1.0_rk)* (maxt2m - &
                    tht)/dtht
            else
                GAMMA_HT = cht
            end if
        else                    ! GAMMA_HT = 1
            GAMMA_HT = 1.0_rk
        end if

    end function GET_GAMMA_HT
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

    real(rk) function GET_GAMMA_LT(lt_opt, mint2m, clt, tlt, dtlt)  result( GAMMA_LT )

        ! !IROUTINE: get_gamma_lt
        !
        ! !DESCRIPTION: Function GET\_GAMMA\_LT computes the  activity factor
        !  associated with low temperature stress of biogenic emission. Called from
        !  GET\_MEGAN\_EMISSIONS only.
        !\\
        !\\
        ! !INTERFACE:

        ! !INPUT PARAMETERS:
        integer,  INTENT(IN) :: lt_opt       ! Option for lt stress calculation
        ! 0=MEGANv3.2 On;
        ! 1=MEGANv3.2 Off
        real(rk), INTENT(IN) :: mint2m         ! minimum input 2-m temperature [K]
        real(rk), INTENT(IN) :: tlt            ! threshold for low temperature [K]
        real(rk), INTENT(IN) :: dtlt           ! delta threshold for low temperature [K]
        real(rk), INTENT(IN) :: clt            ! coefficient for low temperature stress
        !
        ! !RETURN VALUE:
        !REAL(rk)             :: GAMMA_LT  ! LT activity factor [unitless]
        !
        ! !LOCAL VARIABLES:
        REAL(rk)             :: t1         ! combined threshold + delta threshold LT value

        if (lt_opt .eq. 0) then !Calculate GAMMA_LT
            t1 = tlt - dtlt
            if (mint2m >= tlt) then
                GAMMA_LT = 1.0_rk
            else if ( mint2m < tlt .and. mint2m > t1) then
                GAMMA_LT = 1.0_rk + (clt - 1.0_rk)* (tlt - &
                    mint2m)/dtlt
            else
                GAMMA_LT = clt
            end if
        else                    ! GAMMA_LT = 1
            GAMMA_LT = 1.0_rk
        end if

    end function GET_GAMMA_LT
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

    real(rk) function GET_GAMMA_HW(hw_opt, maxws10m, chw, thw, dthw)  result( GAMMA_HW )

        ! !IROUTINE: get_gamma_hw
        !
        ! !DESCRIPTION: Function GET\_GAMMA\_HW computes the  activity factor
        !  associated with high wind stress of biogenic emission. Called from
        !  GET\_MEGAN\_EMISSIONS only.
        !\\
        !\\
        ! !INTERFACE:

        ! !INPUT PARAMETERS:
        integer,  INTENT(IN) :: hw_opt       ! Option for hw stress calculation
        ! 0=MEGANv3.2 On;
        ! 1=MEGANv3.2 Off
        real(rk), INTENT(IN) :: maxws10m       ! maximum input 10-m wind speed [m/s]
        real(rk), INTENT(IN) :: thw            ! threshold for high wind speed [m/s]
        real(rk), INTENT(IN) :: dthw           ! delta threshold for high wind speed [m/s]
        real(rk), INTENT(IN) :: chw            ! coefficient for high wind speed stress
        !
        ! !RETURN VALUE:
        !REAL(rk)             :: GAMMA_HW  ! HT activity factor [unitless]
        !
        ! !LOCAL VARIABLES:
        REAL(rk)             :: t1         ! combined threshold + delta threshold HW value

        if (hw_opt .eq. 0) then !Calculate GAMMA_HW
            t1 = thw + dthw
            if (maxws10m <= thw) then
                GAMMA_HW = 1.0_rk
            else if ( maxws10m > thw .and. maxws10m < t1) then
                GAMMA_HW = 1.0_rk + (chw - 1.0_rk)* (maxws10m - &
                    thw)/dthw
            else
                GAMMA_HW = chw
            end if
        else                    ! GAMMA_HW = 1
            GAMMA_HW = 1.0_rk
        end if

    end function GET_GAMMA_HW
!-----------------------------------------------------------------------

!Initial guesses for canopy profile air temp, pressure, and humidity based on ACCESS
!**********************************************************************************************************************!
! CalcTemp ... compute air temperature (K) at z
!
!**********************************************************************************************************************!
    function CalcTemp(zk, zref, taref, tsref)
        real(rk), intent(in)  :: zk        ! current z (cm)
        real(rk), intent(in)  :: zref      ! reference z (cm) (e.g., = 200 cm if using 2-m temp)
        real(rk), intent(in)  :: taref     ! air temperature at zref (C)  (e.g., 2-m temperature)
        real(rk), intent(in)  :: tsref     ! surface temperature (C) (use model sfc layer or measured skin temp reference)
        real(rk)              :: dtmp      ! temperature gradient (C/cm)
        real(rk)              :: CalcTemp  ! air temp (K) at z

        ! zref and zn must be same units!!!

        dtmp = (taref-tsref)/zref
        CalcTemp = (tsref + dtmp*zk) + 273.15_rk

        return
    end function CalcTemp

!**********************************************************************************************************************!
! CalcPressure ... compute air pressure (mb) at z
!
!**********************************************************************************************************************!
    function CalcPressure(zk, zref, pmbzref, tref, t0)
        real(rk), intent(in)  :: zk        ! current z (cm)
        real(rk), intent(in)  :: zref      ! zref (cm)      (=200 cm if using 2-m temperature)
        real(rk), intent(in)  :: pmbzref   ! air pressure at zref (mb)  (e.g., surface pressure from model, assume = to 2-m)
        real(rk), intent(in)  :: tref      ! temperature at zref (K)  (=sfc_temp if  using press_sfc)
        real(rk), intent(in)  :: t0        ! temperature at z0 (K) (e.g., use top soil layer temperature)
        real(rk), parameter   :: a0=3.42D-04   ! Mg/R (K/cm)
        real(rk)              :: CalcPressure  ! air pressure at current z (mb)

        CalcPressure = pmbzref*exp(-a0*(zk-zref)/(0.5_rk*(tref+t0)))

        return
    end function CalcPressure

!**********************************************************************************************************************!
! function esat - calculate the saturation vapor pressure of water at a
!                 given temperature
!
!    Rogers, R. R., and Yau, M. K. (1989) A Short Course in Cloud Physics,
!    3rd Ed., Butterworth-Heinemann, Burlington, MA. (p16)
!    Valid over -30C <= T <= 35C
!
!**********************************************************************************************************************!
    function esat(tki)
        real(rk), intent(in) :: tki     ! temperature (K)
        real(rk)             :: esat    ! saturation vapor pressure (kPa)
        real(rk)             :: tc      ! temperature (C)

        tc = tki - 273.15_rk
        esat = 0.6112_rk*exp(17.67_rk*tc/(tc+243.5_rk))
        return
    end function esat


    !**********************************************************************************************************************!
! function CalcRelHum - calculates the relative humidity from the
!                             supplied specific humidity, temperature and
!                             pressure
!**********************************************************************************************************************!
    function CalcRelHum(tki, pmbi, qhi)
        real(rk), intent(in) :: tki        ! air temperature (K)
        real(rk), intent(in) :: pmbi       ! air pressure (mb)
        real(rk), intent(in) :: qhi        ! specific humidity (g/kg)
        real(rk)             :: CalcRelHum ! relative humidity (%)
        real(rk)             :: tc, e, es, qhd, pkpa, rhi
        real(rk), parameter  :: rhmin=0.1
        real(rk), parameter  :: rhmax=99.0
        qhd=0.001_rk*qhi
        pkpa=0.1_rk*pmbi
        tc = tki - 273.15_rk
        e = pkpa*qhd/(0.622_rk+qhd)
        es = 0.6112_rk*exp(17.67_rk*tc/(tc+243.5_rk)) !Rogers et al. (1989)
        rhi = max(rhmin, min(rhmax, 100.0_rk*e/es))   ! bound RH to (rhmin, rhmax)
        CalcRelHum = rhi
        return
    end function CalcRelHum

!**********************************************************************************************************************!
! function CalcSpecHum - calculates specific humidity (g/kg) from the supplied relative humidity,
!                             temperature, and pressure
!**********************************************************************************************************************!
    function CalcSpecHum(rhi, tki, pmbi)
        real(rk), intent(in) :: rhi                 ! relative humidity (%)
        real(rk), intent(in) :: tki                 ! air temperature (K)
        real(rk), intent(in) :: pmbi                ! air pressure (mb)
        real(rk)             :: CalcSpecHum         ! specific humidity (g/kg)
        real(rk)             :: es                  ! saturation vapor pressure at tki (mb)
        real(rk)             :: e                   ! ambient vapor pressure (mb)

        es = esat(tki)*10.0_rk            ! kPa -> mb
        e  = es*rhi*0.01_rk               ! mb

        CalcSpecHum = 622.0_rk*e/(pmbi-0.378_rk*e)

        return

    end function CalcSpecHum

!**********************************************************************************************************************!
! CalcCair ... compute cair at z
!
!**********************************************************************************************************************!
    function CalcCair(pmbi, tki)
        real(rk), intent(in)  :: pmbi       ! air pressure at current z (mb)
        real(rk), intent(in)  :: tki        ! air temperature at current z (K)
        real(rk)              :: CalcCair   ! air concentration (molec/cm3)

        CalcCair  = pmbi*7.2428D+18/tki
        return
    end function CalcCair

!**********************************************************************************************************************!
! function Convert_qh_to_h2o - convert qh (g/kg) to h2o (molecs/cm3)
!
!**********************************************************************************************************************!
    function Convert_qh_to_h2o(qhi, cairi)
        real(rk), intent(in) :: qhi                   ! specific humidity at i (g/kg)
        real(rk), intent(in) :: cairi                 ! air concentration at i (molecs/cm3)
        real(rk)             :: Convert_qh_to_h2o     ! h2o concentration (molecs/cm3)

        Convert_qh_to_h2o = 0.001611_rk*qhi*cairi

        return
    end function Convert_qh_to_h2o


!=====================================================================================
!function mdiffstp - set molecular diffusivity data (cm^2/s) for all species at
!                           0 deg C and 1 atm
!=====================================================================================
    subroutine SetMolecDiffSTP(chemmechgas_opt,chemmechgas_tot,mdiffstp)
        integer, intent(in)  :: chemmechgas_opt, chemmechgas_tot      !chemical mechanism and total gas species including transported
        real(rk),dimension(chemmechgas_tot),intent(out)  :: mdiffstp       !molecular diffusivities of species in air at 0degC and 1 atm [cm^2/s]
!    real(kind = dp), parameter         :: mdiffstp_default = 0.100   !default value of mdiffstp (cm^2/s) with no reliable data
!    integer(kind=i4)                   :: l                          !l is species

        if (chemmechgas_opt .eq. 0) then !RACM2
            !Insert MolecDiffSTP Coefficients for RACM2_plus mechanism
            !species in array are:
            != (/NO,   NO2,    O3,  HONO,  HNO4,   HNO3,   N2O5,     CO,  H2O2,  CH4,
            !MO2,   OP1,   MOH,   NO3,   O3P,    O1D,     HO,    HO2,  ORA1,  HAC,
            !PAA, DHMOB, HPALD,  ISHP, IEPOX, PROPNN, ISOPNB, ISOPND, MACRN, MVKN,
            !ISNP /)
            if (chemmechgas_tot .gt. 31) then ! too many species defined
                write(*,*)  'Too many species of ', chemmechgas_tot, ' in namelist for this chemmechgas_opt....exiting'
                write(*,*)  'Set chemmechgas_tot < 31'
                call exit(2)
            else
                mdiffstp =(/0.1802, 0.1361, 0.1444, 0.1349, 0.1041, 0.1041, 0.0808, 0.1807, 0.1300, 0.1952,  &
                    0.1297, 0.1200, 0.1297, 0.1153, 0.2773, 0.2773, 0.2543, 0.2000, 0.1340, 0.1060,  &
                    0.1040, 0.0892, 0.0845, 0.0837, 0.0837, 0.0834, 0.0750, 0.0750, 0.0745, 0.0745,  &
                    0.0712 /)
            end if
        else
            write(*,*)  'Wrong chemical mechanism option of ', chemmechgas_opt, ' in namelist...exiting'
            write(*,*)  'Set chemmechgas_opt to only 0 (RACM2) for now'
            call exit(2)
        end if
!species (NO2, O3)
!    mdiffstp =(/0.1361, 0.1444/)

        return
    end subroutine SetMolecDiffSTP

!==========================================================================
!function MolecDiff - calculate molecular diffusivities (cm^2/s) at a given
!                     temperature and pressure
!==========================================================================
    function MolecDiff(chemmechgas_opt, chemmechgas_tot,ispec, tkx, pmbx)
        integer, intent(in)  :: chemmechgas_opt, chemmechgas_tot      !chemical mechanism and total gas species including transported
        integer, intent(in)         :: ispec     !dummy id for species
        real(rk), intent(in)        :: tkx       !ambient temp [K]
        real(rk), intent(in)        :: pmbx      !ambient press pmb]
        real(rk)                    :: MolecDiff !cm2/s
        real(rk),dimension(chemmechgas_tot)  ::mdiffstp

        call SetMolecDiffSTP(chemmechgas_opt,chemmechgas_tot,mdiffstp)
        MolecDiff = mdiffstp(ispec)*(1013.25_rk/pmbx)*((tkx/298.15_rk)**1.81)

        return
    end function MolecDiff

!=========================================================================
!function mdiffh2o - calculate molecular diffusivity of water vapor in air
!
!source: Tracy (1980)
!=========================================================================
    function mdiffh2o(tki,pmbi)
        real(rk), intent(in) :: tki
        real(rk), intent(in) :: pmbi
        real(rk)             :: mdiffh2o

        mdiffh2o = 0.226_rk*((tki/273.15_rk)**1.81_rk)*(1000.0_rk/pmbi)

        return
    end function mdiffh2o

!======================================================================
!function rs_zhang_gas - calcualte stomatal resistance for trace species
!
!Source - Zhang et al., (2002 & 2003)
!       - 'rsmin' by Wesely et al (1989)
!       - 'bvpd' by Wolfe and Thornton (2011)
!======================================================================
    function rs_zhang_gas(mdiffl, tki, pmbi, ppfdi, srad, relhumi)
        real(rk), intent(in) :: mdiffl            !molecular diffusivity of trace species in air (cm^2/s)
        real(rk), intent(in) :: tki               !air temperature (K)
        real(rk), intent(in) :: pmbi              !air pressure (mb)
        real(rk), intent(in) :: ppfdi             !photosynthetic photon flux (umol/m^2-s)
        real(rk), intent(in) :: srad              !solar irradiation (W/m^2)
        real(rk), intent(in) :: relhumi           !relative humidity (%)
        real(rk)             :: rs_zhang_gas      !stomatal resistance (s/cm)
        !TODO --- Make rs parameters vegtyp dependent -----
        real(rk), parameter  :: rsmin =1.0        !minimum leaf stomatal resistance (s/cm) for deciduous forest
        real(rk), parameter  :: rsmax =10000.0    !maximum leaf stomatal resistance (s/cm) (stoma are closed)
        real(rk), parameter  :: brsp = 196.5      !empirical constant (umol/m^2-s) for deciduous forest
        real(rk), parameter  :: tmin = 0.0        !temperature correction parameter-deciduous forest
        real(rk), parameter  :: tmax = 45.0       !temperature correction parameter-deciduous forest
        real(rk), parameter  :: topt = 27.0       !temperature correction parameter-deciduous forest
        real(rk), parameter  :: bvpd = 0.10       !empirical constant for VPD correction-deciduous forest
        real(rk), parameter  :: phic1 = -1.9      !empirical constant for water stress correction-deciduous forest
        real(rk), parameter  :: phic2 = -2.5      !empirical constant for water stress correction-deciduous forest
        !TODO --- Make rs parameters vegtyp dependent -----
        real(rk)             :: cft
        real(rk)             :: cfvpd
        real(rk)             :: cfphi
        real(rk)             :: tcel
        real(rk)             :: ft1
        real(rk)             :: ft2
        real(rk)             :: et
        real(rk)             :: vpd
        real(rk)             :: phi

        !temperature correction
        tcel = tki - 273.15_rk !convert to degrees C
        if (tcel .ge. tmax) then  !restrict tcell < tmax
            tcel = tmax-0.1_rk
        else if (tcel .le. tmin) then !restrict tcell > tmin
            tcel = tmin+0.1_rk
        else
            tcel = tcel
        end if
        et   = (tmax-topt)/(topt-tmin)
        ft1  =(tcel-tmin)/(topt-tmin)
        ft2  = (tmax-tcel)/(tmax-topt)
        cft  = ft1*(ft2**et)
        !water vapor pressure defit correction
        vpd  = esat(tki)*(1.0_rk - (relhumi/100.0_rk))
        cfvpd= 1.0_rk - bvpd*vpd
        !water stress correction
        phi  = -0.72_rk - 0.0013_rk*srad
        cfphi= (phi-phic2)/(phic1-phic2)

        if (ppfdi >0.0) then
            rs_zhang_gas = rsmin*(1.0_rk+brsp/ppfdi)*mdiffh2o(tki,pmbi)/(mdiffl*cft*cfvpd*cfphi)
        else
            rs_zhang_gas = rsmax                         !nighttime, stoma are closed
        endif

        return
    end function rs_zhang_gas

!===================================================================
!function rbl - calculate leaf boundary resistance for trace species
!
!Source - rb formulation from Wu et al., (2003)
!       - ustar = 0.14*ubar from Weber (1999)
!===================================================================
    function rbl(mdiffl, ubari)
        real(rk)               :: rbl            !leaf boundary layer resistance (s/cm)
        real(rk), intent(in)   :: mdiffl         !molecular diffusivity of species in air (cm^2/s)
        real(rk), intent(in)   :: ubari          !mean wind speed at layer i (cm/s)

!    rbl = 10.53_rk/((mdiffl**0.666667_rk)*max(1.0D-10,ubari))
        rbl = 10.53_rk/((mdiffl**0.666667_rk)*ubari)
        return
    end function rbl

!===============================================================
!function rcl - calculate cuticular resistance for trace species
!
!Source - Wesely (1989)
!===============================================================
    function rcl(hstarl,f01)
        real(rk), intent(in)   :: hstarl         ! effective Henry's law coefficient (M/atm)
        real(rk), intent(in)   :: f01            ! reactivity parameter (0-1)
        real(rk)               :: rcl            ! cuticular resistance (s/cm)
        real(rk), parameter    :: rcref=20.0     ! rc for ozone (s/cm) for deciduous forest

        rcl = rcref/((hstarl*1.0D-05)+f01)

        return
    end function rcl


!===============================================================
!function rml - calculate mesophyll resistance for trace species
!
!Source - Wesely (1989)
!===============================================================
    function rml(hstarl,f01)
        real(rk), intent(in)   :: hstarl          ! effective Henry's law coefficient (M/atm)
        real(rk), intent(in)   :: f01             ! reactivity parameter (0-1)
        real(rk)               :: rml             ! mesophyll resistance (s/cm)

        rml = 1.0_rk/((hstarl/3000.0_rk)+100.0_rk*f01)

        return
    end function rml

!========================================================================
!function EffHenrysLawCoeff - calculate Henry's Law coefficient (M/atm)
!                             have not include temperature dependence yet
!========================================================================
    function EffHenrysLawCoeff(chemmechgas_opt,chemmechgas_tot,ispec)
        integer, intent(in)         :: chemmechgas_opt,chemmechgas_tot ! chemical mechanism and total gas species including transported
        integer, intent(in)         :: ispec  !dummy id for species
        real(rk)                    :: EffHenrysLawCoeff
        real(rk),dimension(chemmechgas_tot) :: hstar

        call SetEffHenrysLawCoeffs(chemmechgas_opt,chemmechgas_tot,hstar)
        EffHenrysLawCoeff = hstar(ispec)

        return
    end function EffHenrysLawCoeff

!====================================================================
!function hstar - set henry's law coefficient for all species (M/atm)
!
!source - Nguyen et al., (2015)
!====================================================================
    subroutine SetEffHenrysLawCoeffs(chemmechgas_opt,chemmechgas_tot,hstar)
        integer, intent(in)                             :: chemmechgas_opt,chemmechgas_tot ! chemical mechanism and total gas species including transported
        real(rk),dimension(chemmechgas_tot),intent(out) :: hstar

        if (chemmechgas_opt .eq. 0) then !RACM2
            !Insert MolecDiffSTP Coefficients for RACM2_plus mechanism
            !species in array are:
            != (/NO,   NO2,    O3,  HONO,  HNO4,   HNO3,   N2O5,     CO,  H2O2,  CH4,
            !MO2,   OP1,   MOH,   NO3,   O3P,    O1D,     HO,    HO2,  ORA1,  HAC,
            !PAA, DHMOB, HPALD,  ISHP, IEPOX, PROPNN, ISOPNB, ISOPND, MACRN, MVKN,
            !ISNP /)
            if (chemmechgas_tot .gt. 31) then ! too many species defined
                write(*,*)  'Too many species of ', chemmechgas_tot, ' in namelist for this chemmechgas_opt...exiting'
                write(*,*)  'Set chemmechgas_tot < 31'
                call exit(2)
            else
                hstar = (/1.9D-03, 1.2D-02,1.03D-02, 2.6D+05, 1.0D+07, 3.2D+13, 1.0D+14, 9.8D-04, 8.4D+04, 1.4D-03,  &
                    6.9D+02, 3.0D+02, 2.2D+02, 3.8D-02, 3.8D-02, 3.8D-02, 3.9D+01, 6.9D+02, 5.6D+03, 2.0D+03,  &
                    5.2D+02, 2.0D+03, 4.0D+04, 7.0D+07, 7.0D+07, 1.0D+04, 5.0D+03, 5.0D+03, 6.0D+03, 6.0D+03,  &
                    5.0D+03 /)
            end if
        else
            write(*,*)  'Wrong chemical mechanism option of ', chemmechgas_opt, ' in namelist...exiting'
            write(*,*)  'Set chemmechgas_opt to only 0 (RACM2) for now'
            call exit(2)
        end if

        return
    end subroutine SetEffHenrysLawCoeffs

!==========================================================================
!function ReactivityParam - calculate reactivity parameters (dimensionless)
!==========================================================================
    function ReactivityParam(chemmechgas_opt,chemmechgas_tot,ispec)
        integer, intent(in)                             :: chemmechgas_opt,chemmechgas_tot ! chemical mechanism and total gas species including transported
        integer, intent(in)        :: ispec   !dummy id for species
        real(rk)                   :: ReactivityParam
        real(rk),dimension(chemmechgas_tot) :: f0

        call SetReactivityParams(chemmechgas_opt,chemmechgas_tot,f0)
        ReactivityParam = f0(ispec)

        return
    end function ReactivityParam

!========================================================
!function f0 - set reactivity parameters for all species
!
!source - Wesely et al., (1989) and Nguyen et al., (2015)
!========================================================
    subroutine SetReactivityParams(chemmechgas_opt,chemmechgas_tot,f0)
!    integer(kind=i4)                              :: l               !l is species
        integer, intent(in)                             :: chemmechgas_opt,chemmechgas_tot ! chemical mechanism and total gas species including transported
        real(rk),dimension(chemmechgas_tot),intent(out) :: f0

        if (chemmechgas_opt .eq. 0) then !RACM2
            !Insert MolecDiffSTP Coefficients for RACM2_plus mechanism
            !species in array are:
            != (/NO,   NO2,    O3,  HONO,  HNO4,   HNO3,   N2O5,     CO,  H2O2,  CH4,
            !MO2,   OP1,   MOH,   NO3,   O3P,    O1D,     HO,    HO2,  ORA1,  HAC,
            !PAA, DHMOB, HPALD,  ISHP, IEPOX, PROPNN, ISOPNB, ISOPND, MACRN, MVKN,
            !ISNP /)
            if (chemmechgas_tot .gt. 31) then ! too many species defined
                write(*,*)  'Too many species of ', chemmechgas_tot, ' in namelist for this chemmechgas_opt....exiting'
                write(*,*)  'Set chemmechgas_tot < 31'
                call exit(2)
            else
                f0 = (/0.0, 0.1, 1.0, 0.1, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,  &
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.0, 0.0,  &
                    0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  &
                    0.0 /)
            end if
        else
            write(*,*)  'Wrong chemical mechanism option of ', chemmechgas_opt, ' in namelist...exiting'
            write(*,*)  'Set chemmechgas_opt to only 0 (RACM2) for now'
            call exit(2)
        end if
        return
    end subroutine SetReactivityParams

end module canopy_utils_mod
