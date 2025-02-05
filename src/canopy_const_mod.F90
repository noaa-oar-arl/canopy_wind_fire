MODULE canopy_const_mod

!-------------------------------------------------------------------------------
! Name:     Constants
! Purpose:  Contains fundamental constants for air quality modeling.
!           15 Jul 2022  Initial Version (P. C. Campbell)
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! References:
!
!      CRC76,        "CRC Handbook of Chemistry and Physics (76th Ed)",
!                     CRC Press, 1995
!      Hobbs, P.V.   "Basic Physical Chemistry for the Atmospheric Sciences",
!                     Cambridge Univ. Press, 206 pp, 1995.
!      Snyder, J.P., "Map Projections-A Working Manual, U.S. Geological Survey
!                     Paper 1395 U.S.GPO, Washington, DC, 1987.
!      Stull, R. B., "An Introduction to Boundary Layer Meteorology", Kluwer,
!                     Dordrecht, 1988
!-------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(15, 307)
    REAL(rk),    PARAMETER :: fillreal  = -9.0e20  ! netCDF _FillValue

    integer, parameter            :: ntotal = 31        !# of total species based on chemical mechanism, including transported species

! Geometric Constants:

    ! pi  (single precision: 3.141593)
    REAL(RK),       PARAMETER     :: pi = 3.14159265358979324_rk

    ! pi/180 [ rad/deg ]
    REAL(RK),          PARAMETER     :: pi180 = pi / 180.0_rk


! Geodetic Constants:

    ! radius of earth [ m ]
    ! -- radius of sphere having same surface area as Clarke ellipsoid of 1866
    !    (Source: Snyder, 1987)
    REAL(RK),          PARAMETER     :: rearth = 6371008.8_rk  ! WGS84 arithmetic mean radius
!    REAL(RK)                         :: rearth

    ! length of a sidereal day [ sec ]  (Source:  CRC76, pp. 14-6)
    REAL(RK),          PARAMETER     :: siday = 86164.09_rk

    ! mean gravitational acceleration [ m/sec**2 ]
    ! --  mean of polar and equatorial values  (Source:  CRC76, pp. 14-6)
    REAL(RK),          PARAMETER     :: grav = 9.80622_rk

    ! latitude degrees to meters
!!!  REAL,          PARAMETER     :: dg2m = rearth * pi180
    REAL(RK)                         :: dg2m

    ! Solar Constant  [ W/m**2 ]  (Source:  CRC76, pp. 14-2)
    REAL(RK),          PARAMETER     :: solcnst = 1373.0


! Fundamental Constants: (Source: CRC76, pp. 1-1 to 1-6)

    ! Avogadro's Constant [ number/mol ]
    REAL(RK),          PARAMETER     :: avo = 6.02214076d23

    ! universal gas constant [ J/mol-K ]
    REAL(RK),          PARAMETER     :: rgasuniv = 8.31446261815324_rk

    ! standard atmosphere [ Pa ]
    REAL(RK),          PARAMETER     :: stdatmpa = 101325.0

    ! standard temperature [ K ]
    REAL(RK),          PARAMETER     :: stdtemp = 273.15_rk

    ! Stefan-Boltzmann [ W/(m**2 K**4) ]
    REAL(RK),          PARAMETER     :: stfblz = 5.67037442d-8


! Non-MKS:

    ! Molar volume of ideal gas at STP [ L/mol ] Non MKS units
    REAL(RK),          PARAMETER     :: molvol = 22.4139695_rk


! Atmospheric Constants:

    ! mean molecular weight for dry air [ g/mol ]
    ! -- 78.06% N2, 21% O2, and 0.943% A on a mole fraction basis
    !    (Source: Hobbs, 1995, pp. 69-70)
    REAL(RK),          PARAMETER     :: mwair = 28.9628_rk

    ! dry-air gas constant [ 287.07548994 J/kg-K ]
    REAL(RK),          PARAMETER     :: rdgas = 1.0e3 * rgasuniv / mwair

    ! mean molecular weight for water vapor [ g/mol ]
    REAL(RK),          PARAMETER     :: mwwat = 18.01528_rk

    ! gas constant for water vapor [ 461.52492604 J/kg-K ]
    REAL(RK),          PARAMETER     :: rwvap = 1.0e3 * rgasuniv / mwwat

    ! FSB NOTE: CPD, CVD, CPWVAP and CVWVAP are calculated assuming dry air and
    ! water vapor are classical ideal gases, i.e. vibration does not contribute
    ! to internal energy.

    ! specific heat of dry air at constant pressure [ 1004.7642148 J/kg-K ]
    REAL(RK),          PARAMETER     :: cpd = 7.0 * rdgas / 2.0

    ! specific heat of dry air at constant volume [ 717.68872485 J/kg-K ]
    REAL(RK),          PARAMETER     :: cvd = 5.0 * rdgas / 2.0

    ! specific heat for water vapor at constant pressure [ 1846.0997042 J/kg-K ]
    REAL(RK),          PARAMETER     :: cpwvap = 4.0 * rwvap

    ! specific heat for water vapor at constant volume [ 1384.5747781 J/kg-K ]
    REAL(RK),          PARAMETER     :: cvwvap = 3.0 * rwvap

    ! vapor press of water at 0 C [ Pa ]  (Source: CRC76 pp. 6-15)
    REAL(RK),          PARAMETER     :: vp0 = 611.29_rk

    ! The following values are taken from p. 641 of Stull (1988):

    ! latent heat of vaporization of water at 0 C [ J/kg ]
    REAL(RK),          PARAMETER     :: lv0 = 2.501e6

    ! Rate of change of latent heat of vaporization w.r.t. temperature [ J/kg-K ]
    REAL(RK),          PARAMETER     :: dlvdt = 2370.0

    ! latent heat of fusion of water at 0 C [ J/kg ]
    REAL(RK),          PARAMETER     :: lf0 = 3.34e5

    ! von Karman constant [ ]
    REAL(RK),          PARAMETER     :: vonk = 0.4_rk

    ! stability parameter beta for neutral conditions [ ]
    ! source:  Bonan et al. (2018)
    ! https://doi.org/10.5194/gmd-11-1467-2018
    REAL(RK),          PARAMETER     :: beta_n = 0.35_rk

    ! e-folding time to be applied to long-term past conditions for biogenic emissions (in days)
    REAL(RK), PARAMETER  :: tau_days  = 5.0_rk

    ! e-folding time to be applied to short-term past conditions for biogenic emissions (in hours)
    REAL(RK), PARAMETER  :: tau_hours = 12.0_rk

END MODULE canopy_const_mod
