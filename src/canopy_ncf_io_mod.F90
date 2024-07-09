MODULE canopy_ncf_io_mod

    public canopy_outncf_alloc, canopy_outncf_init, canopy_outncfglobal, &
        canopy_read_ncf, canopy_write_ncf, canopy_close_files
    private

!-------------------------------------------------------------------------------
! Name:     NetCDF IO
! Purpose:  Contains routines to read NetCDF output.
! Revised:  31 Aug 2004  Original version provided by U. Houston.  (S.-B. Kim)
!           21 Jul 2005  Updated to include error-handling on NetCDF functions
!                        and to conform to NetCDF standard routines for
!                        retrieving data (to eliminate type mismatches).
!                        Updated formatting.  (T. Otte)
!           02 Aug 2005  Changed order of variable declarations in some
!                        subroutines to avoid compile failure on some
!                        machines.  (T. Otte)
!           20 Jun 2006  Removed unused variables.  Changed local variables
!                        FILE to FILENAME and DATA to DATAOUT to avoid
!                        conflicts with F90 keywords.  (T. Otte)
!           19 Apr 2007  Added new routine GET_VAR_REAL2_CDF to read "real"
!                        scalars, as needed for WRFv2.2.  Added new routine
!                        GET_DIM_INT_CDF to retrieve netCDF dimensions.
!                        Changed internal error handling so that errors are
!                        passed back using a non-zero RCODE for dimension
!                        mismatches in addition to netCDF errors.  (T. Otte)
!           12 Feb 2010  Removed unused variable ID_TIME from subroutine
!                        GET_GL_ATT_INT_CDF.  (T. Otte)
!           19 Mar 2010  Removed routines GET_DIMS_CDF, GET_DIM_INT_CDF,
!                        GET_GL_ATT_INT_CDF, GET_GL_ATT_REAL_CDF,
!                        GET_GL_ATT_TEXT_CDF, and GET_DIM_ATT_INT_CDF.
!                        Removed file open and close functions from all
!                        remaining routines, and changed input argument
!                        from FILENAME to CDFID.  (T. Otte)
!           31 Aug 2011  Updated netCDF to F90.  Removed GET_TIMES_CDF.
!                        Changed F77 character declarations to F90 standard.
!                        (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!           02 Feb 2018  Added new routine GET_VAR_3D_INT_CDF.  (T. Spero)
!           22 Jun 2018  Changed module name from WRF_NETCDF to NETCDF_IO.
!                        (T. Spero)
!           19 Dec 2022  Brought module to the Canopy-App (P.C. Campbell)
!-------------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_3d_real_cdf (cdfid, var, dum3d, it, rcode)

        USE canopy_const_mod, ONLY: rk
        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        REAL(rk),          INTENT(OUT)   :: dum3d    ( : , : , : )
        INTEGER                          :: id_data
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER                          :: ny
        INTEGER                          :: nz
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var


        nx = SIZE(dum3d,1)
        ny = SIZE(dum3d,2)
        nz = SIZE(dum3d,3)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, dum3d, start=(/1,1,1,it/),  &
            count=(/nx,ny,nz,1/))

        IF ( rcode /= nf90_noerr ) then
            print*,'read error ',cdfid,var
            print*,'nx,ny,nz=',nx,ny,nz
        endif
    END SUBROUTINE get_var_3d_real_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!    SUBROUTINE get_var_3d_int_cdf (cdfid, var, idum3d, it, rcode)
!
!        USE netcdf
!
!        IMPLICIT NONE
!
!        INTEGER,           INTENT(IN)    :: cdfid
!        INTEGER                          :: id_data
!        INTEGER,           INTENT(OUT)   :: idum3d    ( : , : , : )
!        INTEGER,           INTENT(IN)    :: it
!        INTEGER                          :: nx
!        INTEGER                          :: ny
!        INTEGER                          :: nz
!        INTEGER,           INTENT(OUT)   :: rcode
!        CHARACTER(LEN=*),  INTENT(IN)    :: var
!
!        nx = SIZE(idum3d,2)
!        ny = SIZE(idum3d,1)
!        nz = SIZE(idum3d,3)
!
!        rcode = nf90_inq_varid (cdfid, var, id_data)
!        IF ( rcode /= nf90_noerr ) RETURN
!
!        rcode = nf90_get_var (cdfid, id_data, idum3d, start=(/1,1,1,it/),  &
!            count=(/nx,ny,nz,1/))
!
!    END SUBROUTINE get_var_3d_int_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_2d_real_cdf (cdfid, var, dum2d, it, rcode)

        USE canopy_const_mod, ONLY: rk
        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        REAL(rk),          INTENT(OUT)   :: dum2d    ( : , : )
        INTEGER                          :: id_data
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER                          :: ny
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        nx = SIZE(dum2d,1)
        ny = SIZE(dum2d,2)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) then
            print*,'can not find variable ',trim(var)
            RETURN
        endif
        rcode = nf90_get_var (cdfid, id_data, dum2d, start=(/1,1,it/),  &
            count=(/nx,ny,1/))
        IF ( rcode /= nf90_noerr ) print*,'read error for ',trim(var),&
            ' id_data,it,nx,ny=',id_data,it,nx,ny

    END SUBROUTINE get_var_2d_real_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_2d_int_cdf (cdfid, var, idum2d, it, rcode)

        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        INTEGER                          :: id_data
        INTEGER,           INTENT(OUT)   :: idum2d   ( : , : )
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER                          :: ny
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        nx = SIZE(idum2d,1)
        ny = SIZE(idum2d,2)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, idum2d, start=(/1,1,it/),  &
            count=(/nx,ny,1/))

    END SUBROUTINE get_var_2d_int_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_1d_real_cdf (cdfid, var, dum1d, it, rcode)

        USE canopy_const_mod, ONLY: rk
        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        REAL(rk),          INTENT(OUT)   :: dum1d    ( : )
        INTEGER                          :: id_data
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        nx = SIZE(dum1d)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, dum1d, start=(/1,it/),  &
            count=(/nx,1/))

    END SUBROUTINE get_var_1d_real_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_1d_int_cdf (cdfid, var, idum1d, it, rcode)

        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        INTEGER                          :: id_data
        INTEGER,           INTENT(OUT)   :: idum1d   ( : )
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        nx = SIZE(idum1d)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, idum1d, start=(/1,it/),  &
            count=(/nx,1/))

    END SUBROUTINE get_var_1d_int_cdf

!-------
!    SUBROUTINE get_var_1d_double_cdf (cdfid, var, dum1d, it, rcode)
!
!        USE netcdf
!
!        IMPLICIT NONE
!
!        INTEGER,           INTENT(IN)    :: cdfid
!        REAL,  INTENT(OUT)               :: dum1d    ( : )
!        INTEGER                          :: id_data
!        INTEGER,           INTENT(IN)    :: it
!        INTEGER                          :: nx
!        INTEGER,           INTENT(OUT)   :: rcode
!        CHARACTER(LEN=*),  INTENT(IN)    :: var
!        double precision, allocatable  :: dbtmp(:)
!
!        nx = SIZE(dum1d)
!        allocate(dbtmp(nx))
!
!        rcode = nf90_inq_varid (cdfid, var, id_data)
!        IF ( rcode /= nf90_noerr ) RETURN
!
!        rcode = nf90_get_var (cdfid, id_data, dbtmp, start=(/1,it/),  &
!            count=(/nx,1/))
!        dum1d(:)=sngl(dbtmp(:))
!    END SUBROUTINE get_var_1d_double_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!    SUBROUTINE get_var_real_cdf (cdfid, var, scalar, rcode)
!
!        USE canopy_const_mod, ONLY: rk
!        USE netcdf
!
!        IMPLICIT NONE
!
!        INTEGER,           INTENT(IN)    :: cdfid
!        INTEGER                          :: id_data
!        INTEGER,           INTENT(OUT)   :: rcode
!        REAL(rk),          INTENT(OUT)   :: scalar
!        CHARACTER(LEN=*),  INTENT(IN)    :: var
!
!        rcode = nf90_inq_varid (cdfid, var, id_data)
!        IF ( rcode /= nf90_noerr ) RETURN
!
!        rcode = nf90_get_var (cdfid, id_data, scalar)
!
!    END SUBROUTINE get_var_real_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


    SUBROUTINE canopy_outncf_init

        !-------------------------------------------------------------------------------
        ! Name:     Initialize output NetCDF data structures.
        ! Purpose:  Initializes output NetCDF structures.
        ! Revised:  21 Dec 2022  Initial version.  (P. C. Campbell)
        !-------------------------------------------------------------------------------
        USE canopy_canopts_mod
        use canopy_const_mod, ONLY: fillreal      !constants for canopy models
        use canopy_coord_mod   !main canopy coordinate descriptions
        use canopy_canvars_mod !main canopy variables descriptions

        IMPLICIT NONE

        !-------------------------------------------------------------------------------
        ! Time-independent 1d fields at cell centers.
        !-------------------------------------------------------------------------------

        g_level%fld = fillreal
        g_level%fldname = 'z'
        g_level%long_name = 'model interface level'
        g_level%units = 'm'
        g_level%fillvalue = fillreal
        g_level%dimnames(1) = 'modlays'
        g_level%istart(1) = 1
        g_level%iend(1) = modlays

        !-------------------------------------------------------------------------------
        ! Time-independent 2d fields at cell centers.
        !-------------------------------------------------------------------------------

        g_lat%fld = fillreal
        g_lat%fldname = 'lat'
        g_lat%long_name = 'latitude at cell centers'
        g_lat%units = 'degree_north'
        g_lat%dimnames(1) = 'nlon'
        g_lat%dimnames(2) = 'nlat'
        g_lat%istart(1) = 1
        g_lat%istart(2) = 1
        g_lat%iend(1) = nlon
        g_lat%iend(2) = nlat

        g_lon%fld = fillreal
        g_lon%fldname = 'lon'
        g_lon%long_name = 'longitude at cell centers'
        g_lon%units = 'degree_east'
        g_lon%dimnames(1) = 'nlon'
        g_lon%dimnames(2) = 'nlat'
        g_lon%istart(1) = 1
        g_lon%istart(2) = 1
        g_lon%iend(1) = nlon
        g_lon%iend(2) = nlat

        !-------------------------------------------------------------------------------
        ! Time field.
        !-------------------------------------------------------------------------------

        g_time%fld = fillreal
        g_time%fldname = 'time'
        g_time%calendar_type = 'JULIAN'
        g_time%calendar = 'JULIAN'
        g_time%cartesian_axis = "T"
        g_time%long_name = 'time'
        g_time%units = 'seconds since ' // time_start
        g_time%fillvalue = fillreal
        g_time%dimnames(1) = 'ntime'
        g_time%istart(1) = 1
!        g_time%iend(1) = ntime
        g_time%iend(1) = 1

        !-------------------------------------------------------------------------------
        ! Time-varying 2d fields at cell centers.
        !-------------------------------------------------------------------------------
        c_canheight%fld = fillreal
        c_canheight%fldname = 'canheight'
        c_canheight%long_name = 'canopy height'
        c_canheight%units = 'm'
        c_canheight%fillvalue = fillreal
        c_canheight%dimnames(1) = 'nlon'
        c_canheight%dimnames(2) = 'nlat'
        c_canheight%istart(1) = 1
        c_canheight%istart(2) = 1
        c_canheight%iend(1) = nlon
        c_canheight%iend(2) = nlat

        if (ifcanwind .or. ifcanwaf) then
            c_waf%fld = fillreal
            c_waf%fldname = 'waf'
            c_waf%long_name = 'wind adjustment factor'
            c_waf%units = '1'
            c_waf%fillvalue = fillreal
            c_waf%dimnames(1) = 'nlon'
            c_waf%dimnames(2) = 'nlat'
            c_waf%istart(1) = 1
            c_waf%istart(2) = 1
            c_waf%iend(1) = nlon
            c_waf%iend(2) = nlat

            c_flameh%fld = fillreal
            c_flameh%fldname = 'flameh'
            c_flameh%long_name = 'flame height'
            c_flameh%units = 'm'
            c_flameh%fillvalue = fillreal
            c_flameh%dimnames(1) = 'nlon'
            c_flameh%dimnames(2) = 'nlat'
            c_flameh%istart(1) = 1
            c_flameh%istart(2) = 1
            c_flameh%iend(1) = nlon
            c_flameh%iend(2) = nlat
        end if

        !-------------------------------------------------------------------------------
        ! Time-varying 3d fields at cell centers.
        !-------------------------------------------------------------------------------
        c_lad%fld = fillreal
        c_lad%fldname = 'lad'
        c_lad%long_name = 'leaf area density'
        c_lad%units = 'm2 m-3'
        c_lad%fillvalue = fillreal
        c_lad%dimnames(1) = 'nlon'
        c_lad%dimnames(2) = 'nlat'
        c_lad%dimnames(3) = 'modlays'
        c_lad%istart(1) = 1
        c_lad%istart(2) = 1
        c_lad%istart(3) = 1
        c_lad%iend(1) = nlon
        c_lad%iend(2) = nlat
        c_lad%iend(3) = modlays

        if (ifcanwind .or. ifcanwaf) then
            c_canwind%fld = fillreal
            c_canwind%fldname = 'ws'
            c_canwind%long_name = 'above/below canopy wind speed'
            c_canwind%units = 'm s-1'
            c_canwind%fillvalue = fillreal
            c_canwind%dimnames(1) = 'nlon'
            c_canwind%dimnames(2) = 'nlat'
            c_canwind%dimnames(3) = 'modlays'
            c_canwind%istart(1) = 1
            c_canwind%istart(2) = 1
            c_canwind%istart(3) = 1
            c_canwind%iend(1) = nlon
            c_canwind%iend(2) = nlat
            c_canwind%iend(3) = modlays
        end if
        if (ifcaneddy) then
            c_Kz%fld = fillreal
            c_Kz%fldname = 'kz'
            c_Kz%long_name = 'eddy diffusivity'
            c_Kz%units = 'm2 s-1'
            c_Kz%fillvalue = fillreal
            c_Kz%dimnames(1) = 'nlon'
            c_Kz%dimnames(2) = 'nlat'
            c_Kz%dimnames(3) = 'modlays'
            c_Kz%istart(1) = 1
            c_Kz%istart(2) = 1
            c_Kz%istart(3) = 1
            c_Kz%iend(1) = nlon
            c_Kz%iend(2) = nlat
            c_Kz%iend(3) = modlays
        end if
        if (ifcanphot) then
            c_rjcf%fld = fillreal
            c_rjcf%fldname = 'rjcf'
            c_rjcf%long_name = 'photolysis attenuation correction factor'
            c_rjcf%units = '1'
            c_rjcf%fillvalue = fillreal
            c_rjcf%dimnames(1) = 'nlon'
            c_rjcf%dimnames(2) = 'nlat'
            c_rjcf%dimnames(3) = 'modlays'
            c_rjcf%istart(1) = 1
            c_rjcf%istart(2) = 1
            c_rjcf%istart(3) = 1
            c_rjcf%iend(1) = nlon
            c_rjcf%iend(2) = nlat
            c_rjcf%iend(3) = modlays
        end if
        if (ifcanbio) then
            if (biospec_opt == 0 .or. biospec_opt == 1) then
                c_emi_isop%fld = fillreal
                c_emi_isop%fldname = 'emi_isop'
                c_emi_isop%long_name = 'biogenic isoprene emissions'
                c_emi_isop%units = 'kg m-3 s-1'
                c_emi_isop%fillvalue = fillreal
                c_emi_isop%dimnames(1) = 'nlon'
                c_emi_isop%dimnames(2) = 'nlat'
                c_emi_isop%dimnames(3) = 'nlays'
                c_emi_isop%istart(1) = 1
                c_emi_isop%istart(2) = 1
                c_emi_isop%istart(3) = 1
                c_emi_isop%iend(1) = nlon
                c_emi_isop%iend(2) = nlat
                c_emi_isop%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 2) then
                c_emi_myrc%fld = fillreal
                c_emi_myrc%fldname = 'emi_myrc'
                c_emi_myrc%long_name = 'biogenic myrcene emissions'
                c_emi_myrc%units = 'kg m-3 s-1'
                c_emi_myrc%fillvalue = fillreal
                c_emi_myrc%dimnames(1) = 'nlon'
                c_emi_myrc%dimnames(2) = 'nlat'
                c_emi_myrc%dimnames(3) = 'nlays'
                c_emi_myrc%istart(1) = 1
                c_emi_myrc%istart(2) = 1
                c_emi_myrc%istart(3) = 1
                c_emi_myrc%iend(1) = nlon
                c_emi_myrc%iend(2) = nlat
                c_emi_myrc%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 3) then
                c_emi_sabi%fld = fillreal
                c_emi_sabi%fldname = 'emi_sabi'
                c_emi_sabi%long_name = 'biogenic sabinene emissions'
                c_emi_sabi%units = 'kg m-3 s-1'
                c_emi_sabi%fillvalue = fillreal
                c_emi_sabi%dimnames(1) = 'nlon'
                c_emi_sabi%dimnames(2) = 'nlat'
                c_emi_sabi%dimnames(3) = 'nlays'
                c_emi_sabi%istart(1) = 1
                c_emi_sabi%istart(2) = 1
                c_emi_sabi%istart(3) = 1
                c_emi_sabi%iend(1) = nlon
                c_emi_sabi%iend(2) = nlat
                c_emi_sabi%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 4) then
                c_emi_limo%fld = fillreal
                c_emi_limo%fldname = 'emi_limo'
                c_emi_limo%long_name = 'biogenic limonene emissions'
                c_emi_limo%units = 'kg m-3 s-1'
                c_emi_limo%fillvalue = fillreal
                c_emi_limo%dimnames(1) = 'nlon'
                c_emi_limo%dimnames(2) = 'nlat'
                c_emi_limo%dimnames(3) = 'nlays'
                c_emi_limo%istart(1) = 1
                c_emi_limo%istart(2) = 1
                c_emi_limo%istart(3) = 1
                c_emi_limo%iend(1) = nlon
                c_emi_limo%iend(2) = nlat
                c_emi_limo%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 5) then
                c_emi_care%fld = fillreal
                c_emi_care%fldname = 'emi_care'
                c_emi_care%long_name = 'biogenic 3-carene emissions'
                c_emi_care%units = 'kg m-3 s-1'
                c_emi_care%fillvalue = fillreal
                c_emi_care%dimnames(1) = 'nlon'
                c_emi_care%dimnames(2) = 'nlat'
                c_emi_care%dimnames(3) = 'nlays'
                c_emi_care%istart(1) = 1
                c_emi_care%istart(2) = 1
                c_emi_care%istart(3) = 1
                c_emi_care%iend(1) = nlon
                c_emi_care%iend(2) = nlat
                c_emi_care%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 6) then
                c_emi_ocim%fld = fillreal
                c_emi_ocim%fldname = 'emi_ocim'
                c_emi_ocim%long_name = 'biogenic t-beta-ocimene emissions'
                c_emi_ocim%units = 'kg m-3 s-1'
                c_emi_ocim%fillvalue = fillreal
                c_emi_ocim%dimnames(1) = 'nlon'
                c_emi_ocim%dimnames(2) = 'nlat'
                c_emi_ocim%dimnames(3) = 'nlays'
                c_emi_ocim%istart(1) = 1
                c_emi_ocim%istart(2) = 1
                c_emi_ocim%istart(3) = 1
                c_emi_ocim%iend(1) = nlon
                c_emi_ocim%iend(2) = nlat
                c_emi_ocim%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 7) then
                c_emi_bpin%fld = fillreal
                c_emi_bpin%fldname = 'emi_bpin'
                c_emi_bpin%long_name = 'biogenic beta-pinene emissions'
                c_emi_bpin%units = 'kg m-3 s-1'
                c_emi_bpin%fillvalue = fillreal
                c_emi_bpin%dimnames(1) = 'nlon'
                c_emi_bpin%dimnames(2) = 'nlat'
                c_emi_bpin%dimnames(3) = 'nlays'
                c_emi_bpin%istart(1) = 1
                c_emi_bpin%istart(2) = 1
                c_emi_bpin%istart(3) = 1
                c_emi_bpin%iend(1) = nlon
                c_emi_bpin%iend(2) = nlat
                c_emi_bpin%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 8) then
                c_emi_apin%fld = fillreal
                c_emi_apin%fldname = 'emi_apin'
                c_emi_apin%long_name = 'biogenic alpha-pinene emissions'
                c_emi_apin%units = 'kg m-3 s-1'
                c_emi_apin%fillvalue = fillreal
                c_emi_apin%dimnames(1) = 'nlon'
                c_emi_apin%dimnames(2) = 'nlat'
                c_emi_apin%dimnames(3) = 'nlays'
                c_emi_apin%istart(1) = 1
                c_emi_apin%istart(2) = 1
                c_emi_apin%istart(3) = 1
                c_emi_apin%iend(1) = nlon
                c_emi_apin%iend(2) = nlat
                c_emi_apin%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 9) then
                c_emi_mono%fld = fillreal
                c_emi_mono%fldname = 'emi_mono'
                c_emi_mono%long_name = 'biogenic other monoterpene emissions'
                c_emi_mono%units = 'kg m-3 s-1'
                c_emi_mono%fillvalue = fillreal
                c_emi_mono%dimnames(1) = 'nlon'
                c_emi_mono%dimnames(2) = 'nlat'
                c_emi_mono%dimnames(3) = 'nlays'
                c_emi_mono%istart(1) = 1
                c_emi_mono%istart(2) = 1
                c_emi_mono%istart(3) = 1
                c_emi_mono%iend(1) = nlon
                c_emi_mono%iend(2) = nlat
                c_emi_mono%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 10) then
                c_emi_farn%fld = fillreal
                c_emi_farn%fldname = 'emi_farn'
                c_emi_farn%long_name = 'biogenic alpha-farnesene emissions'
                c_emi_farn%units = 'kg m-3 s-1'
                c_emi_farn%fillvalue = fillreal
                c_emi_farn%dimnames(1) = 'nlon'
                c_emi_farn%dimnames(2) = 'nlat'
                c_emi_farn%dimnames(3) = 'nlays'
                c_emi_farn%istart(1) = 1
                c_emi_farn%istart(2) = 1
                c_emi_farn%istart(3) = 1
                c_emi_farn%iend(1) = nlon
                c_emi_farn%iend(2) = nlat
                c_emi_farn%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 11) then
                c_emi_cary%fld = fillreal
                c_emi_cary%fldname = 'emi_cary'
                c_emi_cary%long_name = 'biogenic beta-caryophyllene emissions'
                c_emi_cary%units = 'kg m-3 s-1'
                c_emi_cary%fillvalue = fillreal
                c_emi_cary%dimnames(1) = 'nlon'
                c_emi_cary%dimnames(2) = 'nlat'
                c_emi_cary%dimnames(3) = 'nlays'
                c_emi_cary%istart(1) = 1
                c_emi_cary%istart(2) = 1
                c_emi_cary%istart(3) = 1
                c_emi_cary%iend(1) = nlon
                c_emi_cary%iend(2) = nlat
                c_emi_cary%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 12) then
                c_emi_sesq%fld = fillreal
                c_emi_sesq%fldname = 'emi_sesq'
                c_emi_sesq%long_name = 'biogenic other sesquiterpene emissions'
                c_emi_sesq%units = 'kg m-3 s-1'
                c_emi_sesq%fillvalue = fillreal
                c_emi_sesq%dimnames(1) = 'nlon'
                c_emi_sesq%dimnames(2) = 'nlat'
                c_emi_sesq%dimnames(3) = 'nlays'
                c_emi_sesq%istart(1) = 1
                c_emi_sesq%istart(2) = 1
                c_emi_sesq%istart(3) = 1
                c_emi_sesq%iend(1) = nlon
                c_emi_sesq%iend(2) = nlat
                c_emi_sesq%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 13) then
                c_emi_mbol%fld = fillreal
                c_emi_mbol%fldname = 'emi_mbol'
                c_emi_mbol%long_name = 'biogenic 232-MBO emissions'
                c_emi_mbol%units = 'kg m-3 s-1'
                c_emi_mbol%fillvalue = fillreal
                c_emi_mbol%dimnames(1) = 'nlon'
                c_emi_mbol%dimnames(2) = 'nlat'
                c_emi_mbol%dimnames(3) = 'nlays'
                c_emi_mbol%istart(1) = 1
                c_emi_mbol%istart(2) = 1
                c_emi_mbol%istart(3) = 1
                c_emi_mbol%iend(1) = nlon
                c_emi_mbol%iend(2) = nlat
                c_emi_mbol%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 14) then
                c_emi_meth%fld = fillreal
                c_emi_meth%fldname = 'emi_meth'
                c_emi_meth%long_name = 'biogenic methanol emissions'
                c_emi_meth%units = 'kg m-3 s-1'
                c_emi_meth%fillvalue = fillreal
                c_emi_meth%dimnames(1) = 'nlon'
                c_emi_meth%dimnames(2) = 'nlat'
                c_emi_meth%dimnames(3) = 'nlays'
                c_emi_meth%istart(1) = 1
                c_emi_meth%istart(2) = 1
                c_emi_meth%istart(3) = 1
                c_emi_meth%iend(1) = nlon
                c_emi_meth%iend(2) = nlat
                c_emi_meth%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 15) then
                c_emi_acet%fld = fillreal
                c_emi_acet%fldname = 'emi_acet'
                c_emi_acet%long_name = 'biogenic acetone emissions'
                c_emi_acet%units = 'kg m-3 s-1'
                c_emi_acet%fillvalue = fillreal
                c_emi_acet%dimnames(1) = 'nlon'
                c_emi_acet%dimnames(2) = 'nlat'
                c_emi_acet%dimnames(3) = 'nlays'
                c_emi_acet%istart(1) = 1
                c_emi_acet%istart(2) = 1
                c_emi_acet%istart(3) = 1
                c_emi_acet%iend(1) = nlon
                c_emi_acet%iend(2) = nlat
                c_emi_acet%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 16) then
                c_emi_co%fld = fillreal
                c_emi_co%fldname = 'emi_co'
                c_emi_co%long_name = 'biogenic carbon monoxide emissions'
                c_emi_co%units = 'kg m-3 s-1'
                c_emi_co%fillvalue = fillreal
                c_emi_co%dimnames(1) = 'nlon'
                c_emi_co%dimnames(2) = 'nlat'
                c_emi_co%dimnames(3) = 'nlays'
                c_emi_co%istart(1) = 1
                c_emi_co%istart(2) = 1
                c_emi_co%istart(3) = 1
                c_emi_co%iend(1) = nlon
                c_emi_co%iend(2) = nlat
                c_emi_co%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 17) then
                c_emi_bvoc%fld = fillreal
                c_emi_bvoc%fldname = 'emi_bvoc'
                c_emi_bvoc%long_name = 'biogenic bidi voc emissions'
                c_emi_bvoc%units = 'kg m-3 s-1'
                c_emi_bvoc%fillvalue = fillreal
                c_emi_bvoc%dimnames(1) = 'nlon'
                c_emi_bvoc%dimnames(2) = 'nlat'
                c_emi_bvoc%dimnames(3) = 'nlays'
                c_emi_bvoc%istart(1) = 1
                c_emi_bvoc%istart(2) = 1
                c_emi_bvoc%istart(3) = 1
                c_emi_bvoc%iend(1) = nlon
                c_emi_bvoc%iend(2) = nlat
                c_emi_bvoc%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 18) then
                c_emi_svoc%fld = fillreal
                c_emi_svoc%fldname = 'emi_svoc'
                c_emi_svoc%long_name = 'biogenic stress voc emissions'
                c_emi_svoc%units = 'kg m-3 s-1'
                c_emi_svoc%fillvalue = fillreal
                c_emi_svoc%dimnames(1) = 'nlon'
                c_emi_svoc%dimnames(2) = 'nlat'
                c_emi_svoc%dimnames(3) = 'nlays'
                c_emi_svoc%istart(1) = 1
                c_emi_svoc%istart(2) = 1
                c_emi_svoc%istart(3) = 1
                c_emi_svoc%iend(1) = nlon
                c_emi_svoc%iend(2) = nlat
                c_emi_svoc%iend(3) = modlays
            end if

            if (biospec_opt == 0 .or. biospec_opt == 19) then
                c_emi_ovoc%fld = fillreal
                c_emi_ovoc%fldname = 'emi_ovoc'
                c_emi_ovoc%long_name = 'biogenic other voc emissions'
                c_emi_ovoc%units = 'kg m-3 s-1'
                c_emi_ovoc%fillvalue = fillreal
                c_emi_ovoc%dimnames(1) = 'nlon'
                c_emi_ovoc%dimnames(2) = 'nlat'
                c_emi_ovoc%dimnames(3) = 'nlays'
                c_emi_ovoc%istart(1) = 1
                c_emi_ovoc%istart(2) = 1
                c_emi_ovoc%istart(3) = 1
                c_emi_ovoc%iend(1) = nlon
                c_emi_ovoc%iend(2) = nlat
                c_emi_ovoc%iend(3) = modlays
            end if

        end if

    END SUBROUTINE canopy_outncf_init


    SUBROUTINE canopy_outncf_alloc

        !-------------------------------------------------------------------------------
        ! Name:     Allocate Arrays for NetCDF output Dimensions
        ! Purpose:  Allocate arrays for NetCDF output Dimensions
        ! Revised:  23 Dec 2022  Original version.  (P. C. Campbell)
        !-------------------------------------------------------------------------------
        USE canopy_canopts_mod
        USE canopy_coord_mod
        USE canopy_canvars_mod

        IMPLICIT NONE

        INTEGER                      :: nn,set_index

        !-------------------------------------------------------------------------------
        ! Time field
        !-------------------------------------------------------------------------------

        nfld1dt = 0

        nfld1dt = nfld1dt + 1  !TIME

        if(.not.allocated(fld1dt)) ALLOCATE ( fld1dt ( nfld1dt ) )

        DO nn = 1, nfld1dt
!            ALLOCATE ( fld1dt(nn)%fld(ntime) )
            ALLOCATE ( fld1dt(nn)%fld(1) )
        ENDDO

        g_time  => fld1dt( 1)

        !-------------------------------------------------------------------------------
        ! Time-independent 1d fields at cell centers.
        !-------------------------------------------------------------------------------

        nfld1dz = 0

        nfld1dz = nfld1dz + 1   !LEVELS

        if(.not.allocated(fld1dz)) ALLOCATE ( fld1dz ( nfld1dz ) )

        DO nn = 1, nfld1dz
            ALLOCATE ( fld1dz(nn)%fld(modlays) )
        ENDDO

        g_level  => fld1dz( 1)

        !-------------------------------------------------------------------------------
        ! Time-independent 2d fields at cell centers.
        !-------------------------------------------------------------------------------

        nfld2dxy = 0

        nfld2dxy = nfld2dxy + 2   !LAT, LON

        if(.not.allocated(fld2dxy)) ALLOCATE ( fld2dxy ( nfld2dxy ) )

        DO nn = 1, nfld2dxy
            ALLOCATE ( fld2dxy(nn)%fld(nlon,nlat) )
        ENDDO

        g_lat     => fld2dxy( 1)
        g_lon     => fld2dxy( 2)

        !-------------------------------------------------------------------------------
        ! Time-varying 2d fields at cell centers.
        !-------------------------------------------------------------------------------

        nfld2dxyt = 0

        nfld2dxyt = nfld2dxyt +1  ! canopy height

        if (ifcanwind .or. ifcanwaf) then
            nfld2dxyt = nfld2dxyt + 1  !WAF
            nfld2dxyt = nfld2dxyt + 1  !FLAMEH
        end if

        if(.not.allocated(fld2dxyt)) ALLOCATE ( fld2dxyt ( nfld2dxyt ) )

        DO nn = 1, nfld2dxyt
            ALLOCATE ( fld2dxyt(nn)%fld(nlon,nlat) )
        ENDDO

        set_index = 1
        c_canheight => fld2dxyt( set_index )
        if (ifcanwind .or. ifcanwaf) then
            set_index = set_index + 1
            c_waf       => fld2dxyt( set_index )
            set_index = set_index + 1
            c_flameh    => fld2dxyt( set_index )
        end if

        !-------------------------------------------------------------------------------
        ! Time-varying 3d fields at cell centers.
        !-------------------------------------------------------------------------------

        nfld3dxyzt = 0

        nfld3dxyzt = nfld3dxyzt + 1 !LAD

        if (ifcanwind .or. ifcanwaf) then
            nfld3dxyzt = nfld3dxyzt + 1 !CANWIND
        end if

        if (ifcaneddy) then
            nfld3dxyzt = nfld3dxyzt + 1 !KZ
        end if

        if (ifcanphot) then
            nfld3dxyzt = nfld3dxyzt + 1 !RJCF
        end if

        if (ifcanbio) then
            if (biospec_opt == 0) then
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_ISOP
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_MYRC
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_SABI
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_LIMO
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_CARE
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_OCIM
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_BPIN
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_APIN
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_MONO
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_FARN
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_CARY
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_SESQ
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_MBOL
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_METH
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_ACET
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_CO
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_BVOC
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_SVOC
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_OVOC
            else
                nfld3dxyzt = nfld3dxyzt + 1 !EMI_SPEC
            end if
        end if

        if(.not.allocated(fld3dxyzt)) ALLOCATE ( fld3dxyzt ( nfld3dxyzt ) )

        DO nn = 1, nfld3dxyzt
            ALLOCATE ( fld3dxyzt(nn)%fld(nlon,nlat,modlays) )
        ENDDO

        set_index = 0

        set_index = set_index + 1
        c_lad    => fld3dxyzt( set_index )

        if (ifcanwind .or. ifcanwaf) then
            set_index = set_index + 1
            c_canwind    => fld3dxyzt( set_index )
        end if

        if (ifcaneddy) then
            set_index = set_index + 1
            c_Kz         => fld3dxyzt( set_index )
        end if

        if (ifcanphot) then
            set_index = set_index + 1
            c_rjcf       => fld3dxyzt( set_index )
        end if

        if (ifcanbio) then
            if (biospec_opt == 0 .or. biospec_opt == 1) then
                set_index = set_index + 1
                c_emi_isop   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 2) then
                set_index = set_index + 1
                c_emi_myrc   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 3) then
                set_index = set_index + 1
                c_emi_sabi   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 4) then
                set_index = set_index + 1
                c_emi_limo   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 5) then
                set_index = set_index + 1
                c_emi_care   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 6) then
                set_index = set_index + 1
                c_emi_ocim   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 7) then
                set_index = set_index + 1
                c_emi_bpin   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 8) then
                set_index = set_index + 1
                c_emi_apin   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 9) then
                set_index = set_index + 1
                c_emi_mono   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 10) then
                set_index = set_index + 1
                c_emi_farn   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 11) then
                set_index = set_index + 1
                c_emi_cary   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 12) then
                set_index = set_index + 1
                c_emi_sesq   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 13) then
                set_index = set_index + 1
                c_emi_mbol   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 14) then
                set_index = set_index + 1
                c_emi_meth   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 15) then
                set_index = set_index + 1
                c_emi_acet   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 16) then
                set_index = set_index + 1
                c_emi_co   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 17) then
                set_index = set_index + 1
                c_emi_bvoc   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 18) then
                set_index = set_index + 1
                c_emi_svoc   => fld3dxyzt( set_index )
            end if
            if (biospec_opt == 0 .or. biospec_opt == 19) then
                set_index = set_index + 1
                c_emi_ovoc   => fld3dxyzt( set_index )
            end if
        end if

    END SUBROUTINE canopy_outncf_alloc

    SUBROUTINE canopy_outncfglobal (cdfid_in, fl)

        !-------------------------------------------------------------------------------
        ! Name:     Output netCDF Global Attributes
        ! Purpose:  Write netCDF global attributes.
        ! Revised:  23 Dec 2022  Original version.  (P. C. Campbell)
        !-------------------------------------------------------------------------------

        USE canopy_coord_mod
        USE canopy_canvars_mod
        USE netcdf

        IMPLICIT NONE

        INTEGER,            INTENT(IN)  :: cdfid_in
        !  CHARACTER(LEN=32)               :: cstr
        CHARACTER(LEN=256), INTENT(IN)  :: fl
        CHARACTER(LEN=32),  PARAMETER   :: pname      = 'CANOPY_OUTNCFGLOBAL'
        INTEGER                         :: rcode
        CHARACTER(LEN=32)               :: var

        !-------------------------------------------------------------------------------
        ! Error, warning, and informational messages.
        !-------------------------------------------------------------------------------

        CHARACTER(LEN=256), PARAMETER :: f9300 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR CREATING ATTRIBUTE FOR', a, &
        & /, 1x, '***   IN FILE ', a, &
        & /, 1x, '***   ', a, &
        & /, 1x, 70('*'))"

        !-------------------------------------------------------------------------------
        ! Define global attributes.
        !-------------------------------------------------------------------------------

        var = "PROGNAME"
        rcode = nf90_put_att (cdfid_in, nf90_global, var, progname)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

        var = "VERSION"
        rcode = nf90_put_att (cdfid_in, nf90_global, var, ver)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

        var = "CODE_DATE"
        rcode = nf90_put_att (cdfid_in, nf90_global, var, vdate)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

        !  var = "INPUT_MODEL"
        !  IF ( met_model == 2 ) THEN
        !    cstr = "WRF ARW " // TRIM(met_release)
        !  ELSE IF ( met_model == 3 ) THEN
        !    cstr = "FV3 " // TRIM(met_release)
        !  ELSE
        !    cstr = " "
        !  ENDIF
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, cstr)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF

        var = "im"
        rcode = nf90_put_att (cdfid_in, nf90_global, var, nlon)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

        var = "jm"
        rcode = nf90_put_att (cdfid_in, nf90_global, var, nlat)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

        var = "levels"
        rcode = nf90_put_att (cdfid_in, nf90_global, var, modlays)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

        var = "times"
        rcode = nf90_put_att (cdfid_in, nf90_global, var, ntime)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

        !  var = "NTHIK"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, nthik)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "GDTYP"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, gdtyp_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "P_ALP"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, p_alp_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "P_BET"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, p_bet_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "P_GAM"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, p_gam_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "XCENT"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, xcent_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "YCENT"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, ycent_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "XORIG"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, xorig_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "YORIG"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, yorig_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "XCELL"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, xcell_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "YCELL"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, ycell_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "VGTYP"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, vgtyp_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "VGTOP"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, vgtop_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "VGTOP"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, vgtop_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "VGLVLS"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, vglvs_gd(:))
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "CEN_LAT"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_cen_lat)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "CEN_LON"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_cen_lon)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "TRUELAT1"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_tru1)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "TRUELAT2"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_tru2)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MOAD_CEN_LAT"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_proj_clat)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_REF_LAT"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_ref_lat)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "STAND_LON"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_proj_clon)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "DX"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_resoln)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "DY"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_resoln)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "PTOP"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_ptop)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_CUMULUS"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_cumulus)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_SHAL_CU"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_shal_cu)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_MICROPHYS"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_expl_moist)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_LW_RAD"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_lw_rad)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_SW_RAD"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_sw_rad)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_PBL"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_pbl)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_SFC_LAY"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_sfc_lay)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_LSM"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_soil_lsm)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_URBAN_PHYS"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_urban_phys)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "LAND_USE_SOURCE"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_lu_src)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_FDDA_3DAN"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_3dan)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  IF ( met_fdda_3dan > 0 ) THEN  ! 3d nudging (any variety)
        !
        !    var = "MET_FDDA_3DAN_WIND"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gv3d)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    var = "MET_FDDA_3DAN_TEMP"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gt3d)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    var = "MET_FDDA_3DAN_MOIS"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gq3d)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    IF ( met_fdda_3dan == 2 ) THEN  ! spectral nudging only
        !      var = "MET_FDDA_3DAN_GEOP"
        !      rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gph3d)
        !      IF ( rcode /= nf90_noerr ) THEN
        !        WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                        TRIM(nf90_strerror(rcode))
        !        CALL graceful_stop (pname)
        !      ENDIF
        !    ENDIF
        !
        !  ENDIF  ! 3d nudging
        !
        !  var = "MET_FDDA_SFAN"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_sfan)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  IF ( met_fdda_sfan > 0 ) THEN  ! surface nudging
        !
        !    var = "MET_FDDA_SFAN_WIND"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gvsfc)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    var = "MET_FDDA_SFAN_TEMP"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gtsfc)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    var = "MET_FDDA_SFAN_MOIS"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gqsfc)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !  ENDIF  ! surface nudging
        !
        !  var = "MET_FDDA_OBS"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_obs)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  IF ( met_fdda_obs > 0 ) THEN  ! observation nudging
        !
        !    var = "MET_FDDA_OBS_WIND"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_giv)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    var = "MET_FDDA_OBS_TEMP"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_git)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    var = "MET_FDDA_OBS_MOIS"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_giq)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !  ENDIF  ! surface nudging
        !
        !  var = "MET_HYBRID"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_hybrid)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "EARTH_RADIUS"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, eradm)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF

    END SUBROUTINE canopy_outncfglobal


    SUBROUTINE canopy_write_ncf (OUTPREFX)

        !-------------------------------------------------------------------------------
        ! Name:     Output netCDF File
        ! Purpose:  Create a netCDF file of Canopy-App Output
        ! Revised:  23 Dec 2022  Original version.  (P.C Campbell Spero)
        !-------------------------------------------------------------------------------

        USE canopy_canopts_mod
        USE canopy_canmet_mod
        USE canopy_coord_mod
        USE canopy_canvars_mod
        USE canopy_files_mod
        USE netcdf

        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT( IN )  :: OUTPREFX
        CHARACTER(LEN=256)              :: fl
        CHARACTER(LEN=16),  PARAMETER   :: pname      = 'CANOPY_WRITE_NCF'
        INTEGER                         :: rcode
        INTEGER                         :: dim_nx
        INTEGER                         :: dim_ny
        INTEGER                         :: dim_nz
        INTEGER                         :: dim_time
        INTEGER                         :: dim_timestr
!        LOGICAL,            SAVE        :: first      = .TRUE.
        INTEGER,  SAVE,     ALLOCATABLE :: id_fld     ( : )
        INTEGER,  SAVE                  :: it         = 0
        INTEGER,            PARAMETER   :: len_time   = 19
        INTEGER                         :: n
        INTEGER                         :: nn
        INTEGER                         :: ntot
        INTEGER                         :: nvars
        CHARACTER(LEN=32)               :: var

        !-------------------------------------------------------------------------------
        ! Error, warning, and informational messages.
        !-------------------------------------------------------------------------------

        CHARACTER(LEN=256), PARAMETER :: f9100 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR CREATING DIMENSION FOR ', a, &
        & /, 1x, '***   IN FILE ', a, &
        & /, 1x, '***   ', a, &
        & /, 1x, 70('*'))"

        CHARACTER(LEN=256), PARAMETER :: f9200 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR DEFINING VARIABLE ', a, &
        & /, 1x, '***   IN FILE ', a, &
        & /, 1x, '***   ', a, &
        & /, 1x, 70('*'))"

        CHARACTER(LEN=256), PARAMETER :: f9300 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR CREATING ATTRIBUTE FOR', a, &
        & /, 1x, '***   IN FILE ', a, &
        & /, 1x, '***   ', a, &
        & /, 1x, 70('*'))"

        CHARACTER(LEN=256), PARAMETER :: f9350 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR ENDING DEFINITIONS ', &
        & /, 1x, '***   IN FILE ', a, &
        & /, 1x, '***   ', a, &
        & /, 1x, 70('*'))"

        CHARACTER(LEN=256), PARAMETER :: f9400 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR WRITING VARIABLE ', a, &
        & /, 1x, '***   TO FILE ', a, &
        & /, 1x, '***   ', a, &
        & /, 1x, 70('*'))"

        CHARACTER(LEN=256), PARAMETER :: f9500 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR CREATING NETCDF FILE', &
        & /, 1x, '***   FILE = ', a, &
        & /, 1x, '***   ', a, &
        & /, 1x, 70('*'))"


        if (infmt_opt .eq. 0) then !Input format is 2D and output is ONLY written to 2D

            !-------------------------------------------------------------------------------
            ! Allocate necessary variables.
            !-------------------------------------------------------------------------------

            it = it + 1

            nvars = nfld1dt + nfld1dz + nfld2dxy + nfld2dxyt + nfld3dxyzt

            IF ( .NOT. ALLOCATED ( id_fld ) ) ALLOCATE ( id_fld ( nvars ) )

            !-------------------------------------------------------------------------------
            ! If first time calling this routine, set up the netCDF output file.
            !-------------------------------------------------------------------------------

!            IF ( first ) THEN

            !-----------------------------------------------------------------------------
            ! Create netCDF file.
            !-----------------------------------------------------------------------------
            fl = TRIM(OUTPREFX)//trim('.nc')

            rcode = nf90_create (fl, nf90_hdf5, cdfid_m)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9500) TRIM(pname), TRIM(fl), TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF

            !-----------------------------------------------------------------------------
            ! Set up dimensions.
            !-----------------------------------------------------------------------------
            var = "time"
            rcode = nf90_def_dim (cdfid_m, TRIM(var), nf90_unlimited, dim_time)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9100) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            var = "timestr"
            rcode = nf90_def_dim (cdfid_m, TRIM(var), len_time, dim_timestr)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9100) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF

            var = "grid_xt"
            rcode = nf90_def_dim (cdfid_m, TRIM(var), nlon, dim_nx)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9100) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF

            var = "grid_yt"
            rcode = nf90_def_dim (cdfid_m, TRIM(var), nlat, dim_ny)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9100) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF

            var = "level"
            rcode = nf90_def_dim (cdfid_m, TRIM(var), modlays, dim_nz)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9100) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF

            !-----------------------------------------------------------------------------
            ! Define variables that will populate the file.
            !-----------------------------------------------------------------------------

            !-------------------------------------------------------------------------------
            ! Time field.
            !-------------------------------------------------------------------------------
            DO n = 1, nfld1dt
                var = TRIM(fld1dt(n)%fldname)
                rcode = nf90_def_var (cdfid_m, TRIM(var), nf90_float,  &
                    (/ dim_time /), id_fld(n))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9200) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit(2)
                ENDIF
            ENDDO
            ntot = nfld1dt

            !-------------------------------------------------------------------------------
            ! Time-independent 1d fields at cell centers.
            !-------------------------------------------------------------------------------

            DO n = 1, nfld1dz
                nn = ntot + n
                var = TRIM(fld1dz(n)%fldname)
                rcode = nf90_def_var (cdfid_m, TRIM(var), nf90_float,  &
                    (/ dim_nz /), id_fld(nn))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9200) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit(2)
                ENDIF
            ENDDO
            ntot = ntot + nfld1dz

            !-------------------------------------------------------------------------------
            ! Time-independent 2d fields at cell centers.
            !-------------------------------------------------------------------------------

            DO n = 1, nfld2dxy
                nn = ntot + n
                var = TRIM(fld2dxy(n)%fldname)
                rcode = nf90_def_var (cdfid_m, TRIM(var), nf90_float,  &
                    (/ dim_nx, dim_ny /), id_fld(nn))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9200) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit(2)
                ENDIF
            ENDDO
            ntot = ntot + nfld2dxy

            !-------------------------------------------------------------------------------
            ! Time-varying 2d fields at cell centers.
            !-------------------------------------------------------------------------------

            DO n = 1, nfld2dxyt
                nn = ntot + n
                var = TRIM(fld2dxyt(n)%fldname)
                rcode = nf90_def_var (cdfid_m, TRIM(var), nf90_float,  &
                    (/ dim_nx, dim_ny, dim_time /), id_fld(nn))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9200) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit(2)
                ENDIF
            ENDDO
            ntot = ntot + nfld2dxyt

            !-------------------------------------------------------------------------------
            ! Time-varying 3d fields at cell centers.
            !-------------------------------------------------------------------------------

            DO n = 1, nfld3dxyzt
                nn = ntot + n
                var = TRIM(fld3dxyzt(n)%fldname)
                rcode = nf90_def_var (cdfid_m, TRIM(var), nf90_float,  &
                    (/ dim_nx, dim_ny, dim_nz, dim_time /), id_fld(nn))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9200) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit(2)
                ENDIF
            ENDDO
            ntot = ntot + nfld3dxyzt

            !-----------------------------------------------------------------------------
            ! Define global attributes.
            !-----------------------------------------------------------------------------

            CALL canopy_outncfglobal (cdfid_m, fl)

            !-----------------------------------------------------------------------------
            ! Define attributes for the variables.
            !-----------------------------------------------------------------------------

            !-------------------------------------------------------------------------------
            ! Time field.
            !-------------------------------------------------------------------------------

            DO n = 1, nfld1dt
                var = TRIM(fld1dt(n)%fldname)
                rcode = nf90_put_att (cdfid_m, id_fld(n), 'calendar',  &
                    TRIM(fld1dt(n)%calendar))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
                rcode = nf90_put_att (cdfid_m, id_fld(n), 'calendar_type',  &
                    TRIM(fld1dt(n)%calendar_type))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
                rcode = nf90_put_att (cdfid_m, id_fld(n), 'cartesian_axis',  &
                    TRIM(fld1dt(n)%cartesian_axis))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
                rcode = nf90_put_att (cdfid_m, id_fld(n), 'long_name',  &
                    TRIM(fld1dt(n)%long_name))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
                rcode = nf90_put_att (cdfid_m, id_fld(n), 'units', TRIM(fld1dt(n)%units))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
            ENDDO
            ntot = nfld1dt

            !-------------------------------------------------------------------------------
            ! Time-independent 1d fields at cell centers.
            !-------------------------------------------------------------------------------

            DO n = 1, nfld1dz
                nn = ntot + n
                var = TRIM(fld1dz(n)%fldname)
                rcode = nf90_put_att (cdfid_m, id_fld(nn), 'long_name',  &
                    TRIM(fld1dz(n)%long_name))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
                rcode = nf90_put_att (cdfid_m, id_fld(nn), 'units', TRIM(fld1dz(n)%units))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
            ENDDO
            ntot = ntot + nfld1dz

            !-------------------------------------------------------------------------------
            ! Time-independent 2d fields at cell centers.
            !-------------------------------------------------------------------------------

            DO n = 1, nfld2dxy
                nn = ntot + n
                var = TRIM(fld2dxy(n)%fldname)
                rcode = nf90_put_att (cdfid_m, id_fld(nn), 'long_name',  &
                    TRIM(fld2dxy(n)%long_name))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
                rcode = nf90_put_att (cdfid_m, id_fld(nn), 'units', TRIM(fld2dxy(n)%units))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
            ENDDO
            ntot = ntot + nfld2dxy

            !-------------------------------------------------------------------------------
            ! Time-varying 2d fields at cell centers.
            !-------------------------------------------------------------------------------

            DO n = 1, nfld2dxyt
                nn = ntot + n
                var = TRIM(fld2dxyt(n)%fldname)
                rcode = nf90_put_att (cdfid_m, id_fld(nn), 'long_name',  &
                    TRIM(fld2dxyt(n)%long_name))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
                rcode = nf90_put_att (cdfid_m, id_fld(nn), 'units',  &
                    TRIM(fld2dxyt(n)%units))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
                rcode = nf90_put_att (cdfid_m, id_fld(nn), '_FillValue',  &
                    fld2dxyt(n)%fillvalue)
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
            ENDDO
            ntot = ntot + nfld2dxyt

            !-------------------------------------------------------------------------------
            ! Time-varying 3d fields at cell centers.
            !-------------------------------------------------------------------------------

            DO n = 1, nfld3dxyzt
                nn = ntot + n
                var = TRIM(fld3dxyzt(n)%fldname)
                rcode = nf90_put_att (cdfid_m, id_fld(nn), 'long_name',  &
                    TRIM(fld3dxyzt(n)%long_name))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
                rcode = nf90_put_att (cdfid_m, id_fld(nn), 'units',  &
                    TRIM(fld3dxyzt(n)%units))
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
                rcode = nf90_put_att (cdfid_m, id_fld(nn), '_FillValue',  &
                    fld3dxyzt(n)%fillvalue)
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
            ENDDO
            ntot = ntot + nfld3dxyzt

            !-------------------------------------------------------------------------------
            ! Take file out of "define mode" so that variables may be written to it.
            !-------------------------------------------------------------------------------

            rcode = nf90_enddef (cdfid_m)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9350) TRIM(pname), TRIM(fl), TRIM(nf90_strerror(rcode))
                CALL exit (2)
            ENDIF

!            ENDIF  ! first = .TRUE.

            !-------------------------------------------------------------------------------
            ! Assign pointer variables to the respective arrays from canopy-app.
            !-------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------
            ! Time-independent 1d fields at cell centers.
            !-------------------------------------------------------------------------------

            g_level%fld = zk

            !-------------------------------------------------------------------------------
            ! Time-independent 2d fields at cell centers.
            !-------------------------------------------------------------------------------

            g_lat%fld = variables_2d%lat
            g_lon%fld = variables_2d%lon

            !-------------------------------------------------------------------------------
            ! Time field.
            !-------------------------------------------------------------------------------

            g_time%fld = (it-1)*time_intvl

            !-------------------------------------------------------------------------------
            ! Time-varying 2d fields at cell centers.
            !-------------------------------------------------------------------------------
            c_canheight%fld = variables_2d%ch
            if (ifcanwind .or. ifcanwaf) then
                c_waf%fld = waf_2d
                c_flameh%fld = flameh_2d
            end if

            !-------------------------------------------------------------------------------
            ! Time-varying 3d fields at cell centers.
            !-------------------------------------------------------------------------------
            c_lad%fld  = lad_3d
            if (ifcanwind .or. ifcanwaf) then
                c_canwind%fld  = canWIND_3d
            end if
            if (ifcaneddy) then
                c_Kz%fld       = Kz_3d
            end if
            if (ifcanphot) then
                c_rjcf%fld     = rjcf_3d
            end if
            if (ifcanbio) then
                if (biospec_opt == 0 .or. biospec_opt == 1) then
                    c_emi_isop%fld = emi_isop_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 2) then
                    c_emi_myrc%fld = emi_myrc_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 3) then
                    c_emi_sabi%fld = emi_sabi_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 4) then
                    c_emi_limo%fld = emi_limo_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 5) then
                    c_emi_care%fld = emi_care_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 6) then
                    c_emi_ocim%fld = emi_ocim_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 7) then
                    c_emi_bpin%fld = emi_bpin_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 8) then
                    c_emi_apin%fld = emi_apin_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 9) then
                    c_emi_mono%fld = emi_mono_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 10) then
                    c_emi_farn%fld = emi_farn_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 11) then
                    c_emi_cary%fld = emi_cary_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 12) then
                    c_emi_sesq%fld = emi_sesq_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 13) then
                    c_emi_mbol%fld = emi_mbol_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 14) then
                    c_emi_meth%fld = emi_meth_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 15) then
                    c_emi_acet%fld = emi_acet_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 16) then
                    c_emi_co%fld   = emi_co_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 17) then
                    c_emi_bvoc%fld = emi_bvoc_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 18) then
                    c_emi_svoc%fld = emi_svoc_3d
                end if
                if (biospec_opt == 0 .or. biospec_opt == 19) then
                    c_emi_ovoc%fld = emi_ovoc_3d
                end if
            end if

            !-------------------------------------------------------------------------------
            ! Write variables.
            !-------------------------------------------------------------------------------

            write(*,*)  'Writing NetCDF Output'
            write(*,*)  '-------------------------------'

            !-------------------------------------------------------------------------------
            ! Time field.
            !-------------------------------------------------------------------------------
            write(*,*)  'Writing Time field'
            write(*,*)  '-------------------------------'

            DO n = 1, nfld1dt
                var = TRIM(fld1dt(n)%fldname)
                rcode = nf90_put_var (cdfid_m, id_fld(n), fld1dt(n)%fld,  &
                    start = (/ 1 /),  &
                    count = (/ 1 /) )
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9400) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
            ENDDO
            ntot = nfld1dt

!            IF ( first ) THEN  ! write time-independent fields

            !-------------------------------------------------------------------------------
            ! Time-independent 1d fields at cell centers.
            !-------------------------------------------------------------------------------
            write(*,*)  'Writing Time-independent 1d fields'
            write(*,*)  '-------------------------------'
!            ENDIF  ! first

            DO n = 1, nfld1dz
                nn = ntot + n
                var = TRIM(fld1dz(n)%fldname)
!                IF ( first ) THEN  ! write time-independent fields
                rcode = nf90_put_var (cdfid_m, id_fld(nn), fld1dz(n)%fld)
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9400) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
!                ENDIF  ! first
            ENDDO

            ntot = ntot + nfld1dz

            !-------------------------------------------------------------------------------
            ! Time-independent 2d fields at cell centers.
            !-------------------------------------------------------------------------------
!            IF ( first ) THEN  ! write time-independent fields
            write(*,*)  'Writing Time-independent 2d fields'
            write(*,*)  '-------------------------------'
!            ENDIF  ! first
            DO n = 1, nfld2dxy
                nn = ntot + n
                var = TRIM(fld2dxy(n)%fldname)
!                IF ( first ) THEN  ! write time-independent fields
                rcode = nf90_put_var (cdfid_m, id_fld(nn), fld2dxy(n)%fld)
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9400) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
!                ENDIF  ! first
            ENDDO

            ntot = ntot + nfld2dxy

            !-------------------------------------------------------------------------------
            ! Time-varying 2d fields at cell centers.
            !-------------------------------------------------------------------------------
            write(*,*)  'Writing Time-varying 2d fields'
            write(*,*)  '-------------------------------'
            DO n = 1, nfld2dxyt
                nn = ntot + n
                var = TRIM(fld2dxyt(n)%fldname)
                rcode = nf90_put_var (cdfid_m, id_fld(nn), fld2dxyt(n)%fld,  &
                    start = (/ 1, 1, 1 /),  &
                    count = (/ fld2dxyt(n)%iend(1), &
                    fld2dxyt(n)%iend(2), 1 /) )
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9400) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
            ENDDO
            ntot = ntot + nfld2dxyt

            !-------------------------------------------------------------------------------
            ! Time-varying 3d fields at cell centers.
            !-------------------------------------------------------------------------------
            write(*,*)  'Writing Time-varying 3d fields'
            write(*,*)  '-------------------------------'

            DO n = 1, nfld3dxyzt
                nn = ntot + n
                var = TRIM(fld3dxyzt(n)%fldname)
                rcode = nf90_put_var (cdfid_m, id_fld(nn), fld3dxyzt(n)%fld,  &
                    start = (/ 1, 1, 1, 1 /),  &
                    count = (/ fld3dxyzt(n)%iend(1), &
                    fld3dxyzt(n)%iend(2), &
                    fld3dxyzt(n)%iend(3), 1 /) )
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (6,f9400) TRIM(pname), TRIM(var), TRIM(fl),  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit (2)
                ENDIF
            ENDDO
            ntot = ntot + nfld3dxyzt

        end if !infmt_opt

!        first = .FALSE.
    END SUBROUTINE canopy_write_ncf


    SUBROUTINE canopy_read_ncf(infile)

        !-------------------------------------------------------------------------------
        ! Name:     Read Canopy Met/Sfc Inputs from NetCDF
        ! Purpose:  Read Canopy Met/Sfc Inputs from NetCDF
        ! Revised:  19 Dec 2022  Original version.  (P.C. Campbell)
        !-------------------------------------------------------------------------------
        USE canopy_canopts_mod !main canopy option descriptions
        use canopy_const_mod, ONLY: rk
        USE canopy_coord_mod
        USE canopy_canmet_mod
        USE netcdf

        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT( IN )   :: infile
        CHARACTER(LEN=*), PARAMETER      :: pname      = 'canopy_read_ncf'
        CHARACTER(LEN=24)                :: rddate
        CHARACTER(LEN=32)                :: date_init
        INTEGER                          :: cdfid,rcode,varid,it
        double precision                 :: rdtime
        REAL(rk)                         :: variables_3d_real(nlon,nlat,var3d_set)
        REAL(rk)                         :: variables_2d_real(nlon,nlat)
        INTEGER                          :: variables_2d_int(nlon,nlat)
        REAL(rk)                         :: variables_1d_lev_real(var3d_set)
!        REAL(rk)                         :: variables_real(nlon*nlat)
!        INTEGER                          :: variables_int(nlon*nlat)

        !-------------------------------------------------------------------------------
        ! Error, warning, and informational messages.
        !-------------------------------------------------------------------------------

        !    CHARACTER(LEN=256), PARAMETER :: f9400 = "(/, 1x, 70('*'), &
        !    & /, 1x, '*** SUBROUTINE: ', a, &
        !    & /, 1x, '***   ERROR RETRIEVING VARIABLE FROM NCF INPUT FILE', &
        !    & /, 1x, '***   VARIABLE = ', a, &
        !    & /, 1x, '***   RCODE = ', a, &
        !    & /, 1x, 70('*'))"

        CHARACTER(LEN=256), PARAMETER :: f9410 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR RETRIEVING NCF ID FROM NCF INPUT FILE', &
        & /, 1x, '***   VARIABLE = ', a, &
        & /, 1x, '***   NCF: ', a, &
        & /, 1x, 70('*'))"

        !    CHARACTER(LEN=256), PARAMETER :: f9420 = "(/, 1x, 70('*'), &
        !    & /, 1x, '*** SUBROUTINE: ', a, &
        !    & /, 1x, '***   ERROR INQUIRING ABOUT VAR IN NCF INPUT FILE', &
        !    & /, 1x, '***   VARIABLE = ', a, &
        !    & /, 1x, '***   NCF: ', a, &
        !    & /, 1x, 70('*'))"

        !    CHARACTER(LEN=256), PARAMETER :: f9430 = "(/, 1x, 70('*'), &
        !    & /, 1x, '*** SUBROUTINE: ', a, &
        !    & /, 1x, '***   ERROR RETRIEVING DIMS FROM NCF INPUT FILE', &
        !    & /, 1x, '***   VARIABLE = ', a, &
        !    & /, 1x, '***   NCF: ', a, &
        !    & /, 1x, 70('*'))"

        CHARACTER(LEN=256), PARAMETER :: f9440 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR OPENING NCF INPUT FILE', &
        & /, 1x, '***   FILE = ', a, &
        & /, 1x, '***   NCF: ', a, &
        & /, 1x, 70('*'))"

        CHARACTER(LEN=256), PARAMETER :: f9450 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR RETRIEVING VAR ATT IN NCF INPUT FILE', &
        & /, 1x, '***   VARIABLE = ', a, &
        & /, 1x, '***   NCF: ', a, &
        & /, 1x, 70('*'))"

        !-------------------------------------------------------------------------------
        ! Open NCF input file
        !-------------------------------------------------------------------------------

        rcode = nf90_open (infile, nf90_nowrite, cdfid)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (*,f9440) TRIM(pname), TRIM(infile)
            call exit(2)
        ENDIF

        !-------------------------------------------------------------------------------
        ! Extract and check input lat/lon dimensions and time information.
        !-------------------------------------------------------------------------------
        !time information
        rcode = nf90_inq_varid (cdfid, 'time', varid)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (*,f9410) TRIM(pname), 'time',  &
                TRIM(nf90_strerror(rcode))
            CALL exit(2)
        ENDIF

        rcode = nf90_get_att (cdfid, varid, 'units', date_init)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (*,f9450) TRIM(pname), 'SIMULATION_START_DATE',  &
                TRIM(nf90_strerror(rcode))
            CALL exit(2)
        ENDIF

        rddate = date_init(13:31) // '.0000'
        rddate(11:11)='-'

        rcode=nf90_get_var(cdfid,varid,rdtime)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (*,f9410) TRIM(pname), 'time',  &
                TRIM(nf90_strerror(rcode))
            CALL exit(2)
        ENDIF

        !grid 1D or 2D  lat/lon variables
!        it = int(rdtime) + 1
        it = 1  !Assume input files are full integer time steps

        if (infmt_opt .eq. 0) then !Input format is 2D

            CALL get_var_2d_real_cdf (cdfid, 'lat', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'lon',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%lat = variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%lat=reshape(variables_2d%lat,[size(variables_2d%lat)])
            CALL get_var_2d_real_cdf (cdfid, 'lon', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'lon',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%lon=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%lon=reshape(variables_2d%lon,[size(variables_2d%lon)])
            !Canopy input met/sfc variables
            !Clumping index
            CALL get_var_2d_real_cdf (cdfid, 'clu', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'clu',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%clu=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%clu=reshape(variables_2d%clu,[size(variables_2d%clu)])
            !Cosine of solar zenith angle
            CALL get_var_2d_real_cdf (cdfid, 'csz', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'csz',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%csz=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%csz=reshape(variables_2d%csz,[size(variables_2d%csz)])
            !Forest Fraction
            CALL get_var_2d_real_cdf (cdfid, 'canfrac', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'canfrac',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%canfrac=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%canfrac=reshape(variables_2d%canfrac,[size(variables_2d%canfrac)])
            !Forest canopy height
            CALL get_var_2d_real_cdf (cdfid, 'ch', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'ch',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%ch=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%ch=reshape(variables_2d%ch,[size(variables_2d%ch)])
            !Fire Radiative Power
            CALL get_var_2d_real_cdf (cdfid, 'frp', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'frp',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%frp=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%frp=reshape(variables_2d%frp,[size(variables_2d%frp)])
            !Reference height above canopy
            CALL get_var_2d_real_cdf (cdfid, 'href', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'href',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%href=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%href=reshape(variables_2d%href,[size(variables_2d%href)])
            !Leaf Area Index
            CALL get_var_2d_real_cdf (cdfid, 'lai', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'lai',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%lai=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%lai=reshape(variables_2d%lai,[size(variables_2d%lai)])
            !Monin-Obukhov Length
            CALL get_var_2d_real_cdf (cdfid, 'mol', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'mol',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%mol=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%mol=reshape(variables_2d%mol,[size(variables_2d%mol)])
            !Friction velocity
            CALL get_var_2d_real_cdf (cdfid, 'fricv', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'fricv',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%fricv=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%fricv=reshape(variables_2d%fricv,[size(variables_2d%fricv)])
            !Reference U Wind Speed (at HREF)
            CALL get_var_2d_real_cdf (cdfid, 'ugrd10m', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'ugrd10m',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%ugrd10m=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%ugrd10m=reshape(variables_2d%ugrd10m,[size(variables_2d%ugrd10m)])
            !Reference V Wind Speed (at HREF)
            CALL get_var_2d_real_cdf (cdfid, 'vgrd10m', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'vgrd10m',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%vgrd10m=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%vgrd10m=reshape(variables_2d%vgrd10m,[size(variables_2d%vgrd10m)])
            !Surface (veg+soil) Roughness Length
            CALL get_var_2d_real_cdf (cdfid, 'sfcr', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'sfcr',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%sfcr=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%sfcr=reshape(variables_2d%sfcr,[size(variables_2d%sfcr)])
            !Vegetation Type
            CALL get_var_2d_int_cdf (cdfid, 'vtype', variables_2d_int, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'vtype',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%vtype=variables_2d_int
            !Also reshape to 1D array for 1D calculation and output
!            variables%vtype=reshape(variables_2d%vtype,[size(variables_2d%vtype)])
            !Soil Type
            CALL get_var_2d_int_cdf (cdfid, 'sotyp', variables_2d_int, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'sotyp',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%sotyp=variables_2d_int
            !Also reshape to 1D array for 1D calculation and output
!            variables%sotyp=reshape(variables_2d%sotyp,[size(variables_2d%sotyp)])
            !Surface pressure
            CALL get_var_2d_real_cdf (cdfid, 'pressfc', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'pressfc',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%pressfc=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%pressfc=reshape(variables_2d%pressfc,[size(variables_2d%pressfc)])
            !instantaneous surface downward shortwave flux
            CALL get_var_2d_real_cdf (cdfid, 'dswrf', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'dswrf',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%dswrf=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%dswrf=reshape(variables_2d%dswrf,[size(variables_2d%dswrf)])
            !instantaneous surface sensible heat net flux
            CALL get_var_2d_real_cdf (cdfid, 'shtfl', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'shtfl',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%shtfl=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%shtfl=reshape(variables_2d%shtfl,[size(variables_2d%shtfl)])
            !Surface temperature
            CALL get_var_2d_real_cdf (cdfid, 'tmpsfc', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'tmpsfc',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%tmpsfc=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%tmpsfc=reshape(variables_2d%tmpsfc,[size(variables_2d%tmpsfc)])
            !2-meter temperature
            CALL get_var_2d_real_cdf (cdfid, 'tmp2m', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'tmp2m',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%tmp2m=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%tmp2m=reshape(variables_2d%tmp2m,[size(variables_2d%tmp2m)])
            !2-meter specific humidity
            CALL get_var_2d_real_cdf (cdfid, 'spfh2m', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'spfh2m',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%spfh2m=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%spfh2m=reshape(variables_2d%spfh2m,[size(variables_2d%spfh2m)])
            !Height of planetary boundary layer
            CALL get_var_2d_real_cdf (cdfid, 'hpbl', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'hpbl',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%hpbl=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%hpbl=reshape(variables_2d%hpbl,[size(variables_2d%hpbl)])
            !Mass precipitation rate
            CALL get_var_2d_real_cdf (cdfid, 'prate_ave', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'prate_ave',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%prate_ave=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%prate_ave=reshape(variables_2d%prate_ave,[size(variables_2d%prate_ave)])
            !Volumetric soil moisture level 1
            CALL get_var_2d_real_cdf (cdfid, 'soilw1', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'soilw1',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%soilw1=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%soilw1=reshape(variables_2d%soilw1,[size(variables_2d%soilw1)])
            !Volumetric soil moisture level 2
            CALL get_var_2d_real_cdf (cdfid, 'soilw2', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'soilw2',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%soilw2=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%soilw2=reshape(variables_2d%soilw2,[size(variables_2d%soilw2)])
            !Volumetric soil moisture level 3
            CALL get_var_2d_real_cdf (cdfid, 'soilw3', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'soilw3',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%soilw3=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%soilw3=reshape(variables_2d%soilw3,[size(variables_2d%soilw3)])
            !Volumetric soil moisture level 4
            CALL get_var_2d_real_cdf (cdfid, 'soilw4', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'soilw4',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%soilw4=variables_2d_real
            !Also reshape to 1D array for 1D calculation and output
!            variables%soilw4=reshape(variables_2d%soilw4,[size(variables_2d%soilw4)])
            !Wilting point
            CALL get_var_2d_real_cdf (cdfid, 'wilt', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'wilt',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%wilt=variables_2d_real

            !Also reshape to 1D array for 1D calculation and output
!            variables%wilt=reshape(variables_2d%wilt,[size(variables_2d%wilt)])
            !Ozone W126
            CALL get_var_2d_real_cdf (cdfid, 'ozone_w126', variables_2d_real, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'ozone_w126',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            variables_2d%ozone_w126=variables_2d_real

            !Also reshape to 1D array for 1D calculation and output
!            variables%ozone_w126=reshape(variables_2d%ozone_w126,[size(variables_2d%ozone_w126)])

            !3D Input Level Profile
            if (var3d_opt .eq. 1) then
                CALL get_var_1d_real_cdf (cdfid, 'lev', variables_1d_lev_real, it, rcode)
                IF ( rcode /= nf90_noerr ) THEN
                    WRITE (*,f9410) TRIM(pname), 'lev',  &
                        TRIM(nf90_strerror(rcode))
                    CALL exit(2)
                ENDIF
                variables_1d%lev=variables_1d_lev_real
            end if
            !GEDI 3D PAVD Profile
            if (pavd_opt .eq. 1 ) then
                if (var3d_opt .ne. 1) then !check to make sure variables_3d are allocated
                    write(*,*)  'Wrong choice of VAR3D_OPT ', var3d_opt, ' in namelist...exiting'
                    call exit(2)
                else !read PAVD
                    !Forest Plant Area Volume Density (PAVD) Profile
                    CALL get_var_3d_real_cdf (cdfid, 'pavd', variables_3d_real, it, rcode)
                    IF ( rcode /= nf90_noerr ) THEN
                        WRITE (*,f9410) TRIM(pname), 'pavd',  &
                            TRIM(nf90_strerror(rcode))
                        CALL exit(2)
                    ENDIF
                    variables_3d%pavd=variables_3d_real
                end if
            end if



        else if (infmt_opt .eq. 1) then !Input format is 1D

            CALL get_var_1d_real_cdf (cdfid, 'lat', variables%lat, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'lat',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF

            CALL get_var_1d_real_cdf (cdfid, 'lon', variables%lon, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'lon',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF

            !Canopy input met/sfc variables
            !Clumping index
            CALL get_var_1d_real_cdf (cdfid, 'clu', variables%clu, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'clu',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Cosine of solar zenith angle
            CALL get_var_1d_real_cdf (cdfid, 'csz', variables%csz, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'csz',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Forest Fraction
            CALL get_var_1d_real_cdf (cdfid, 'canfrac', variables%canfrac, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'canfrac',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Forest canopy height
            CALL get_var_1d_real_cdf (cdfid, 'ch', variables%ch, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'ch',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Fire Radiative Power
            CALL get_var_1d_real_cdf (cdfid, 'frp', variables%frp, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'frp',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Reference height above canopy
            CALL get_var_1d_real_cdf (cdfid, 'href', variables%href, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'href',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Leaf Area Index
            CALL get_var_1d_real_cdf (cdfid, 'lai', variables%lai, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'lai',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Monin-Obukhov Length
            CALL get_var_1d_real_cdf (cdfid, 'mol', variables%mol, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'mol',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Friction velocity
            CALL get_var_1d_real_cdf (cdfid, 'fricv', variables%fricv, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'fricv',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Reference U Wind Speed (at HREF)
            CALL get_var_1d_real_cdf (cdfid, 'ugrd10m', variables%ugrd10m, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'ugrd10m',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Reference V Wind Speed (at HREF)
            CALL get_var_1d_real_cdf (cdfid, 'vgrd10m', variables%vgrd10m, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'vgrd10m',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Surface (veg+soil) Roughness Length
            CALL get_var_1d_real_cdf (cdfid, 'sfcr', variables%sfcr, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'sfcr',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Vegetation Type
            CALL get_var_1d_int_cdf (cdfid, 'vtype', variables%vtype, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'vtype',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Soil Type
            CALL get_var_1d_int_cdf (cdfid, 'sotyp', variables%sotyp, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'sotyp',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Surface pressure
            CALL get_var_1d_real_cdf (cdfid, 'pressfc', variables%pressfc, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'pressfc',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !instantaneous surface downward shortwave flux
            CALL get_var_1d_real_cdf (cdfid, 'dswrf', variables%dswrf, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'dswrf',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !instantaneous surface sensible heat net flux
            CALL get_var_1d_real_cdf (cdfid, 'shtfl', variables%shtfl, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'shtfl',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Surface temperature
            CALL get_var_1d_real_cdf (cdfid, 'tmpsfc', variables%tmpsfc, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'tmpsfc',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !2-meter temperature
            CALL get_var_1d_real_cdf (cdfid, 'tmp2m', variables%tmp2m, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'tmp2m',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !2-meter specific humidity
            CALL get_var_1d_real_cdf (cdfid, 'spfh2m', variables%spfh2m, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'spfh2m',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Height of planetary boundary layer
            CALL get_var_1d_real_cdf (cdfid, 'hpbl', variables%hpbl, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'hpbl',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Mass precipitation rate
            CALL get_var_1d_real_cdf (cdfid, 'prate_ave', variables%prate_ave, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'prate_ave',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
        else
            write(*,*)  'Wrong INFMT_OPT choice of ', infmt_opt, ' in namelist...exiting'
            call exit(2)
        end if !Input Format (1D or 2D)

    END SUBROUTINE canopy_read_ncf

    SUBROUTINE canopy_close_files(OUTPREFX)

!-------------------------------------------------------------------------------
! Name:     Close Files
! Purpose:  Close I/O API files.
! Revised:  10 Sep 2001  Original version.  (T. Otte)
!           09 Jan 2002  Changed "stop" statements to calls to "m3exit" for
!                        graceful shut-down of I/O API files.  (T. Otte)
!           29 Aug 2011  Replaced module IODECL3 with I/O API module M3UTILIO.
!                        Replaced F77 character declarations with F90 standard.
!                        Improved error handling.  (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!           19 Dec 2018  Added runtime option to choose output format.
!                        (T. Spero)
!           04 Aug 2023  Using for canopy-app closing output files (P. Campbell)
!-------------------------------------------------------------------------------

        use canopy_files_mod
        use netcdf

        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT( IN )    :: OUTPREFX
        CHARACTER(LEN=256)                :: fl
        CHARACTER(LEN=16),  PARAMETER     :: pname      = 'CLOSE_FILES'
        INTEGER                           :: rcode

!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

        CHARACTER(LEN=256), PARAMETER :: f9100 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR CLOSING NETCDF FILE', &
        & /, 1x, '***   FILE = ', a, &
        & /, 1x, '***   ', a, &
        & /, 1x, 70('*'))"

!-------------------------------------------------------------------------------
! Gracefully close output files.
!-------------------------------------------------------------------------------

        fl = TRIM(OUTPREFX)//trim('.nc')

        rcode = nf90_close (cdfid_m)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9100) TRIM(pname), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

    END SUBROUTINE canopy_close_files

END MODULE canopy_ncf_io_mod
