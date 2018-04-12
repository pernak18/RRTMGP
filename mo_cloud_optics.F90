! Module: mo_cloud_optics
!
! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe, Jennifer Delamere, and Rick Pernak?
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
! BSD 3-clause license, 
! see http://opensource.org/licenses/BSD-3-Clause
!
!
! Description:
! This is the interface for routines that receive cloud physical 
! properties and return cloud optical properties by band using
! Pade formulations. 

module mo_cloud_optics

  ! force declaration of all variables ("good practice")
  implicit none

  use mo_rte_kind, only: wp
  use mo_optical_props, only: ty_optical_props_1scl, &
    ty_optical_props_2scl, ty_optical_props_nscl
  use mo_cloud_optical_props, only: ty_cloud_optical_props_arry, &
    ty_cloud_optical_props_1scl, ty_cloud_optical_props_2str, &
    ty_cloud_optical_props_nstr
  use netcdf

  interface read_field
    module procedure read_scalar, read_1d_field, read_2d_field, &
      read_3d_field, read_4d_field
  end interface

  ! ---------------------------------------------------------------
  ! Attributes for entire class

  ! dimensions
  integer :: ncol, nlay, nlev, nband, nrough, nsizereg, &
    nsizeice, nsizeliq, ncoeff_ext, ncoeff_ssa, ncoeff_asy

  ! cloud physical properties
  integer :: roughness
  real(wp), dimension(:,:), allocatable :: &
    cld_fraction, cld_liq_water_path, cld_ice_water_path, &
    cld_tot_water_path, r_eff_liq, r_eff_ice

  ! metadata
  character(len=2) :: spectral_domain
  character(len=4) :: approx, method
  character(len=*) :: cld_coeff_file
  logical :: do_lw, do_pade

  !https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/295774
  ! declaring this way means i don't have to call allocate() as long
  ! as i set an option at compile time?
  ! Pade coefficients/LUT data
  real(wp), dimension(:,:,:), allocatable :: &
    pade_extliq, pade_ssaliq, pade_asyliq
  real(wp), dimension(:,:,:,:), allocatable :: &
    pade_extice, pade_ssaice, pade_asyice
  real(wp), dimension(:,:), allocatable :: &
    lut_extliq, lut_ssaliq, lut_asyliq
  real(wp), dimension(:,:,:), allocatable :: &
    lut_extice, lut_ssaice, lut_asyice

  ! masks for cloud presence
  logical, dimension(:,:), allocatable :: &
    ice_cloud_present, liq_cloud_present

  ! cloud optical properties
  real(wp), dimension(:,:,:), allocatable :: extliq, ssaliq, asyliq
  real(wp), dimension(:,:,:,:), allocatable :: extice, ssaic, asyice

  ! Pade size (effective radius) regime boundaries
  integer, dimension(4) :: &
    r_eff_bounds_liq_ext, r_eff_bounds_ice_ext, &
    r_eff_bounds_liq_ssa, r_eff_bounds_ice_ssa, &
    r_eff_bounds_liq_asy, r_eff_bounds_ice_asy
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  ! Methods for entire class
  public :: error_check, make_cloud_mask, &
    load_cloud_optics, calc_optical_properties, &
    stop_on_err, get_dim_length, create_var, dim_exists, var_exists, &
    read_field, write_field
  ! ---------------------------------------------------------------

  type, extends(ty_optical_props_1scl), public :: &
    ty_cloud_optics_band_1scl
    ! Attributes in extended (sub)class
    real(wp), dimension(:,:,:), allocatable :: tau
  end type ty_cloud_optics_band

  type, extends(ty_optical_props_2str), public :: &
    ty_cloud_optics_band_2str
    ! Attributes in extended (sub)class
    real(wp), dimension(:,:,:), allocatable :: tau, ssa, asy
  end type ty_cloud_optics_band

  type, extends(ty_optical_props_nstr), public :: &
    ty_cloud_optics_band_nstr
    ! Attributes in extended (sub)class
    real(wp), dimension(:,:,:), allocatable :: tau, ssa, asy
    real(wp), dimension(:,:,:,:), allocatable :: phi
  end type ty_cloud_optics_band

  contains

  function init(this, gas_spec, cld_fraction, &
    play, plev, tlay, tsfc, &
    cld_liq_water_path, cld_ice_water_path, ice_rough, &
    r_eff_liq, r_eff_ice, method, is_lw, is_pade, delta_scale)

    ! initialization of the cloud_optics class
    ! ---------------------------------------------------------------
    !   this -- object of class ty_cloud_optics_*_band
    !   gas_spec -- object of class ty_gas_optics_specification, 
    !     which contains state information at each layer and level
    !   cld_fraction -- (nCol x nLay) float array of percentage 
    !     of column occupied by cloud
    !   play -- (nCol x nLay) float array of layer pressures [Pa, mb]
    !   plev -- (nCol x nLev) float array of level pressures [Pa, mb]
    !   tlay -- (nCol x nLay) float array of layer temperatures [K]
    !   tsfc -- nCol float vector of surface skin temperatures [K]
    !   cld_???_water_path -- (nCol x nLay) float arrays of ice and 
    !     liquid water path
    !   ice_rough -- Ice surface roughness category needed for 
    !     Yang (2013) ice optics parameterization 
    !     (1 = none, 2 = medium, 3 = high)
    !   r_eff_??? -- (nCol x nLay) float arrays of ice and liquid 
    !     effective radii (or size) [microns]
    !   method -- string, method (1-scalar, 2-stream, n-stream) used 
    !     for solving radiative transfer equation with scattering 
    !     ("1scl", "2str", or "nstr" only)
    !   is_lw -- boolean, are the specifications for the LW domain?
    !     (otherwise SW)
    !   is_pade -- boolean, is the Pade approximation being used in
    !     the optical property calculations? (otherwise LookUp Table
    !     -- LUT -- data are used)
    !   delta_scale -- boolean
    ! ---------------------------------------------------------------

    ! ---------------------------------------------------------------

    ! force declaration of all variables ("good practice")
    implicit none

    ! ---------------------------------------------------------------
    ! User input
    ! what should i be doing with gas_spec? i never use it
    class(ty_gas_optics_specification), intent(in) :: gas_spec

    real(wp), dimension(:,:), intent(in) :: play, plev, tlay
    real(wp), dimension(:), intent(in) :: tsfc
    integer, intent(in) :: icergh
    logical, intent(in) :: is_lw, do_pade, delta_scale

    real(wp), intent(in), dimension(:,:), allocatable :: &
      cld_fraction, cld_liq_water_path, cld_ice_water_path, &
      r_eff_liq, r_eff_ice
    ! ---------------------------------------------------------------

    ! ---------------------------------------------------------------
    ! Output is an object of either the 1scl, 2str, or nstr class
    ! can i use a conditional in the declaration for this?
    character(len=4), intent(in):: method
    character(len=128): err_msg
    if (method .eq. "1scl") then
      class(ty_cloud_optics_band_1scl), intent(out) :: this
    elseif (method .eq. "2str") then
      class(ty_cloud_optics_band_2str), intent(out) :: this
    elseif (method .eq. "nstr") then
    else
      class(ty_cloud_optics_band_nstr), intent(out) :: this
      err_msg = "Please specify '1scl', '2str', or 'nstr' for method"
      call stop_on_err(err_msg)
    endif
    ! ---------------------------------------------------------------

    ! ---------------------------------------------------------------
    ! Attribute assignments
    ! dimensions
    this%ncol = size(play, 1)
    this%nlay = size(play, 2)
    this%nlev = size(plev, 2)

    ! "metadata" properties
    this%method = method
    this%do_lw = is_lw
    this%do_pade = is_pade

    ! can these boolean-dependent assignments be more compact?
    if (is_lw) then
      this%spectral_domain = 'lw'
    else
      this%spectral_domain = 'sw'
    endif

    if (is_pade) then
      this%approx = 'pade'
    else
      this%approx = 'lut'
    endif

    this%delta_scale = delta_scale

    this%cld_coeff_file = 'rrtmgp-' // &
      this%spectral_domain // '-inputs-cloud-optics-' // &
      trim(this%approx) //'.nc'

    ! physical property arrays
    ! need to do some validation on these guys in this%error_check!!
    this%cld_fraction = cld_fraction
    this%cld_liq_water_path = cld_liq_water_path
    this%cld_ice_water_path = cld_ice_water_path
    this%r_eff_liq = r_reff_liq
    this%r_eff_ice = r_eff_ice
    this%roughness = ice_rough
    do icol = 1, this%ncol
      do ilyr = 1, this%nlayers
        this%cld_tot_water_path(icol, ilyr) = &
          cld_liq_water_path(icol, ilyr) + &
          cld_ice_water_path(icol, ilyr)
      enddo ! icol
    enddo ! ilyr

    ! load LUT or Pade Coefficients
    this%load_cloud_optics()

    this%error_check()

    ! are clouds present?
    this%make_cloud_mask()

    ! compute optical properties
    this%calc_optical_properties()
    ! ---------------------------------------------------------------
  end function init

  function load_cloud_optics(this)
    ! from an input netCDF, load cloud optics data into 
    ! class attributes (should work regardless of the spectral domain 
    ! and Pade-LUT formalism)

    ! local variables will just be used as shortcuts
    integer, private :: ncol, nlay, nlev, nband, nrough, nsizereg, &
      ncoeff_ext, ncoeff_ssa_asy
    character(len=2), private :: lw_sw
    character(len=4), private :: pade_lut
    character(len=*), private :: extLiqStr, ssaLiqStr, asyLiqStr
    character(len=*), private :: extIceStr, ssaIceStr, asyIceStr

    lw_sw = this%spectral_domain
    pade_lut = this%approx

    ! Pade/LUT data extraction
    if (nf90_open(&
       trim(this%cld_coeff_file), NF90_WRITE, ncid) &
       /= NF90_NOERR) &
       call stop_on_err("load_cloud_optics(): can't open file " // &
       trim(cld_coeff_file))

    ! LUT and Pade common dimensions, with shortcuts
    this%nband = get_dim_length(ncid, 'nband_' // lw_sw)
    this%nrough = get_dim_length(ncid, 'nrghice')
    nBand = this%nband
    nRough = this%ncoeff_nrough

    extLiqStr = pade_lut // '_extliq_' // lw_sw
    ssaLiqStr = pade_lut // '_ssaliq_' // lw_sw
    asyLiqStr = pade_lut // '_asyliq_' // lw_sw
    extIceStr = pade_lut // '_extice_' // lw_sw
    ssaIceStr = pade_lut // '_ssaice_' // lw_sw
    asyIceStr = pade_lut // '_asyice_' // lw_sw

    if (this%do_pade) then
      ! Pade-exclusive dimensions and shortcuts
      this%nsizereg = get_dim_length(ncid, 'nsizereg')
      this%ncoeff_ext = get_dim_length(ncid, 'ncoeff_ext')
      this%ncoeff_ssa = get_dim_length(ncid, 'ncoeff_ssa_g')
      this%ncoeff_asy = get_dim_length(ncid, 'ncoeff_ssa_g')
      nSizeReg = this%nsizereg
      nCoeffExt = this%ncoeff_ext
      nCoeffSSA = this%ncoeff_ssa
      nCoeffAsy = this%ncoeff_asy

      ! Pade coefficient input arrays
      ! i don't think i need to do the allocation anymore since 
      ! i'm declaring these guys with "allocatable"

      ! liquid
      this%pade_extliq  = this%read_field(&
        ncid, extLiqStr, nCoeffExt, nSizeReg, nBand)
      this%pade_ssaliq  = this%read_field(&
        ncid, ssaLiqStr, nCoeffSSA, nSizeReg, nBand)
      this%pade_asyliq  = this%read_field(&
        ncid, asyLiqStr, nCoeffAsy, nSizeReg, nBand)

      ! ice
      this%pade_extice  = this%read_field(&
        ncid, extIceStr, nRough, nCoeffExt, nSizeReg, nBand)
      this%pade_ssaice  = this%read_field(&
        ncid, ssaIceStr, nRough, nCoeffSSA, nSizeReg, nBand)
      this%pade_asyice  = this%read_field(&
        ncid, asyIceStr, nRough, nCoeffAsy, nSizeReg, nBand)

      ! Particle size regimes for Pade formulations
      ! these bounds were determined by fitting the Pade polynomials 
      ! to the lookup table data (the fits had to be done in groups of 
      ! of effective radius)
      ! SW and LW are mostly the same (defaults are LW)
      this%r_eff_bounds_liq_ext = [2, 10, 35, 60]
      this%r_eff_bounds_liq_ssa = [2, 10, 35, 60]
      this%r_eff_bounds_liq_asy = [2, 9, 20, 60]
      this%r_eff_bounds_ice_ext = [10, 20, 30, 180]
      this%r_eff_bounds_ice_ssa = [10, 20, 30, 180]
      this%r_eff_bounds_ice_asy = [10, 20, 30, 180]

      if (this%do_lw) then
        this%r_eff_bounds_liq_asy(2) = 8
        this%r_eff_bounds_ice_asy(2) = 50
      endif ! LW
    else ! LUT
      ! LUT-exclusive dimensions and shortcuts
      this%nsizeice = get_dim_length(ncid, 'nsize_ice')
      this%nsizeliq = get_dim_length(ncid, 'nsize_liq')
      nSizeIce = this%nsizeice
      nSizeLiq = this%nsizeliq

      ! LUT data arrays
      ! i don't think i need to do the allocation anymore since 
      ! i'm declaring these guys with "allocatable"

      ! liquid
      this%lut_extliq  = &
        this%read_field(ncid, extLiqStr, nSizeLiq, nBand)
      this%lut_ssaliq  = &
        this%read_field(ncid, ssaLiqStr, nSizeLiq, nBand)
      this%lut_asyliq  = &
        this%read_field(ncid, asyLiqStr, nSizeLiq, nBand)

      ! ice
      this%lut_extice  = this%read_field(&
        ncid, extIceStr, nRough, nBand, nSizeIce)
      this%lut_ssaice  = this%read_field(&
        ncid, ssaIceStr, nRough, nBand, nSizeIce)
      this%lut_asyice  = this%read_field(&
        ncid, asyIceStr, nRough, nBand, nSizeIce)
    endif ! Pade/LUT

    ncid = nf90_close(ncid)
  end function load_cloud_optics

  function error_check(this)
    ! check the data for any invalid values

    character(len=128) :: error_msg

    error_msg = ''

    icergh = this%roughness
    if (icergh < 1 .or. icergh > 3) then
      error_msg = 'cloud optics: ' \\
        'cloud ice surface roughness flag is out of bounds'
      return
    endif

  end function error_check

  function make_cloud_mask(this)
    ! are values in cloud optics object valid?
    ! are clouds present?
    ! populate logical arrays based on whether physical properties 
    ! pass a series of conditionals

    logical, private, allocatable, dimension(:,:) :: &
      lCldBool, iCldBool

    ! cloud water path (liquid + ice), effective sizes
    real(wp) :: cwp, rei, rel

    ! effective size bounds
    integer, dimension(4) :: liqbounds, icebounds

    ! minimum value for cloud quantities
    real(wp), parameter :: cldmin = 1.e-20     

    icebounds = this%r_eff_bounds_ice_ext
    liqbounds = this%r_eff_bounds_liq_asy
    do icol = 1, this%ncol
      do ilyr = 1, this%nlayers
        cwp = this%cld_tot_water_path(icol, ilyr)
        rel = this%r_eff_liq(icol, ilyr)
        rei = this%r_eff_ice(icol, ilyr)

        ! liquid cloud logical array conditionals
        lCldBool(icol, ilyr) = merge(.true., .false., &
          (this%cldfrac(icol,ilyr) .gt. cldmin .and. cwp .gt. cldmin))

        lCldBool(icol, ilyr) = merge(.true., .false., &
          (rel .gt. 0.0_wp &
           this%cld_liq_water_path(icol, ilyr) .gt. 0.0_wp) )

        ! the rel conditional isn't exactly what it was in the 
        ! original code (2.5 <= rel <= 20), but i'm using the r_eff
        ! bounds here and not "magic numbers"
        lCldBool(icol, ilyr) = merge(.true., .false., &
          (rel .ge. liqbounds[1] .and. rel .le. liqbounds[3]) )

        ! ice cloud logical array conditionals
        iCldBool(icol, ilyr) = merge(.true., .false., &
          (this%cldfrac(icol,ilyr) .gt. cldmin .and. cwp .gt. cldmin))
        iCldBool(icol, ilyr) = merge(.true., .false., &
          (rei .gt. 0.0_wp &
           this%cld_ice_water_path(icol, ilyr) .gt. 0.0_wp) )
        iCldBool(icol, ilyr) = merge(.true., .false., &
          (rei .ge. icebounds[1] .and. rei .le. icebounds[4]) )
      enddo ! icol
    enddo ! ilyr

    this%ice_cloud_present = iCldBool
    this%liq_cloud_present = lCldBool

  end function make_cloud_mask()

  function calc_optical_properties(this)
    ! from cloud optics object, calculate optical properties:
    !       tau/OD/extinction/ext, 
    !       omega/w/single scatter albedo/ssa, 
    !       g/asymmetry/asy
    ! given the specifications in its attributes

    ! ------- Local -------

    ! cloud effective sizes (microns)
    real(wp) :: radliq, radice

    integer :: icol, ilyr, ibnd, irade, iradw, iradg

    ! for LUT coefficients and indexing
    integer, icergh, irad
    float arg, slope, var1, var2

    ! extinction coefficient, single scattering albedo, and 
    ! asymmetry arrays for both phases
    real(wp), dimension(:,:), allocatable :: &
      extliq, extice, ssaliq, ssaice, asyliq, asyice

    ! liquid and ice cloud extinction optical depth, no delta scaling
    real(wp) :: tauliq, tauice

    ! ice and liquid scattering terms, (liquid + ice?) asymmetry
    real(wp) :: scatice, scatliq, asycld

    ! Forward scattering terms
    real(wp) :: forwliq, forwice

    ! delta scaling extinction OD and single scattering albedo
    real(wp) :: tauliq_del, tauice_del, ssaliq_del, ssaice_del

    icergh = this%roughness

    do icol = 1, this%ncol
      do ilyr = 1, this%nlayers
        ! liquid
        if (this%liq_cloud_present(icol, ilyr)) then
          radliq = this%r_eff_liq(icol, ilyr)

          do ibnd = 1, this%nband
            if (this%do_pade) then
              ! determine index of effective size domain in which 
              ! radliq exists
              irade = this%get_irad(radliq, .true., 'ext')
              iradw = this%get_irad(radliq, .true., 'ssa')
              iradg = this%get_irad(radliq, .true., 'asy')

              ! compute ext, ssa, and asy using Pade approximation
              ! think of how to get rid of these magic numbers (number
              ! of numerator and denominator terms -- 3 and 2)
              extliq(ilyr, ibnd) = this%pade(&
                this%pade_extliq(:, irade, ibnd), radliq, 3, 3)
              ssaliq(ilyr, ibnd) = 1-this%pade(&
                this%pade_ssaliq(:, iradw, ibnd), radliq, 3, 2)
              asyliq(ilyr, ibnd) = this%pade(&
                this%pade_asyliq(:, iradg, ibnd), radliq, 3, 2)
            else
              ! not entirely sure about these arithmetic, but i think
              ! this is pretty much what was in the original code
              arg = radliq - 1.5_wp
              irad = int(arg)
              slope = arg - real(irad)

              ! compute ext, ssa, and asy by linearly interpolating
              ! to LUT data
              var1 = this%lut_extliq(ibnd, irad)
              var2 = this%lut_extliq(ibnd, irad+1)
              extliq(ilyr, ibnd) = slope * (var2-var1) + var1

              var1 = this%lut_ssaliq(ibnd, irad)
              var2 = this%lut_ssaliq(ibnd, irad+1)
              ssaliq(ilyr, ibnd) = slope * (var2-var1) + var1

              var1 = this%lut_asyliq(ibnd, irad)
              var2 = this%lut_asyliq(ibnd, irad+1)
              asyliq(ilyr, ibnd) = slope * (var2-var1) + var1
            endif ! Pade liquid
          enddo ! band
        else
          ! fill values if no cloud present
          extliq(ilyr,ibnd) = 0.0_wp
          ssaliq(ilyr,ibnd) = 0.0_wp
          asyliq(ilyr,ibnd) = 0.0_wp
        endif ! liquid cloud

        ! Ice
        if (this%liq_cloud_present(icol, ilyr)) then
          radice = this%r_eff_ice(icol,ilyr)

          do ibnd = 1, nband
            if (this%do_pade) then
              ! determine index of effective size domain in which 
              ! radliq exists
              irade = this%get_irad(radice, .false., 'ext')
              iradw = this%get_irad(radice, .false., 'ssa')
              iradg = this%get_irad(radice, .false., 'asy')

              ! Derive optical properties for selected ice roughness
              extice(ilyr, ibnd) = this%pade(&
                this%pade_extice(icergh,:,irade, ibnd), radice, 3, 3)
              ssaice(ilyr, ibnd) = 1-this%pade(&
                this%pade_ssaice(icergh,:,iradw, ibnd), radice, 3, 2)
              asyice(ilyr, ibnd) = this%pade(&
                this%pade_asyice(icergh,:,iradg, ibnd), radice, 3, 2)
            else
              arg = radice / 10.0_wp
              irad = int(arg)
              slope = arg - real(irad)

              ! compute ext, ssa, and asy by linearly interpolating
              ! to LUT data
              var1 = this%lut_extice(icergh, ibnd, irad)
              var2 = this%lut_extice(icergh, ibnd, irad+1)
              extliq(ilyr, ibnd) = slope * (var2-var1) + var1

              var1 = this%lut_ssaice(icergh, ibnd, irad)
              var2 = this%lut_ssaice(icergh, ibnd, irad+1)
              ssaliq(ilyr, ibnd) = slope * (var2-var1) + var1

              var1 = this%lut_asyice(icergh, ibnd, irad)
              var2 = this%lut_asyice(icergh, ibnd, irad+1)
              asyliq(ilyr, ibnd) = slope * (var2-var1) + var1
            endif ! Pade ice
          enddo ! end band loop
        else
          ! fill values if no cloud present
          extliq(ilyr, ibnd) = 0.0_wp
          ssaliq(ilyr, ibnd) = 0.0_wp
          asyliq(ilyr, ibnd) = 0.0_wp
        endif ! liquid cloud

        ! Combine liquid and ice contributions into total cloud 
        ! optical properties
        ! by now, the calculations are independent of the path we 
        ! used to get here (Pade or LUT)
        do ibnd = 1, this%nbnd
          tauliq = this%cld_liq_water_path(icol, ilyr) * &
            extliq(ilyr, ibnd)
          tauice = this%cld_ice_water_path(icol, ilyr) * &
            extice(ilyr, ibnd)

          ! still doing this provision?
          !if (tauice == 0.0_wp .and. tauliq == 0.0_wp) then
          !  error_msg = &
          !    'cloud optics: cloud optical depth is zero'
          !  return
          !endif

          if (this%delta_scale) then
            ! LW and SW with no delta scaling looked the same to me
            ! so i put them in the same block
            scatliq = ssaliq(ilyr, ibnd) * tauliq
            scatice = ssaice(ilyr, ibnd) * tauice

            select type(this)
              type is (ty_cloud_optics_band_1scl)
                this%tau(icol, ilyr, ibnd) = tauice + tauliq
              type is (ty_cloud_optics_band_2str)
                this%tau(icol, ilyr, ibnd) = tauice + tauliq
                this%ssa(icol, ilyr, ibnd) = (scatice + scatliq) / &
                  this%tau(icol, ilyr, ibnd)
                this%asy(icol, ilyr, ibnd) = &
                  (scatice * asyice(ilyr, ibnd) + &
                  (scatliq * asyliq(ilyr, ibnd))) / &
                  (scatice + scatliq)
              type is (ty_cloud_optics_band_nstr)
                this%tau(icol, ilyr, ibnd) = tauice + tauliq
                this%ssa(icol, ilyr, ibnd) = (scatice + scatliq) / &
                  this%tau(icol, ilyr, ibnd)
                asycld = (scatice * asyice(ilyr, ibnd) + &
                  (scatliq * asyliq(ilyr, ibnd))) / &
                  (scatice + scatliq)
                this%phi(1, icol, ilyr, ibnd) = 1.0_wp
                this%phi(2, icol, ilyr, ibnd) = asycld
                this%phi(3, icol, ilyr, ibnd) = asycld * asycld
            end select
          else
            ! SHORTWAVE, delta-scaling
            ! no quality assurance yet on whether i got the formulas
            ! correct (should probably find a reference for these)
            forwliq = asyliq(ilyr, ibnd) * asyliq(ilyr, ibnd)
            forwice = asyice(ilyr, ibnd) * asyice(ilyr, ibnd)
            tauliq_del = (1.0_wp - forwliq * ssaliq(ilyr, ibnd)) * &
              tauliq
            tauice_del = (1.0_wp - forwice * ssaice(ilyr, ibnd)) * &
              tauice

            select type(this)
              type is (ty_cloud_optics_band_1scl)
                this%tau(icol, ilyr, ibnd) = tauliq_del + tauice_del
              type is (ty_cloud_optics_band_2str)
                this%tau(icol, ilyr, ibnd) = tauliq_del + tauice_del
                ssaliq_del = ssaliq(ilyr, ibnd) * &
                  (1._wp - forwliq) / &
                  (1._wp - forwliq * ssaliq(ilyr, ibnd))
                ssaice_del = ssaice(ilyr, ibnd) * &
                  (1._wp - forwice) / &
                  (1._wp - forwice * ssaice(ilyr, ibnd))
                scatliq = ssaliq_del * tauliq_del
                scatice = ssaice_del * tauice_del
                this%ssa(icol, ilyr, ibnd) = (scatliq + scatice) / &
                  this%tau(icol,ilyr,ibnd)
                this%asy(icol, ilyr, ibnd) = &
                  (scatice * (asyice(ilyr, ibnd) - forwice) / &
                  (1._wp - forwice) + &
                  scatliq * (asyliq(ilyr, ibnd) - forwliq) / &
                  (1._wp - forwliq)) / (scatice + scatliq)
              type is (ty_cloud_optics_band_nstr)
                this%tau(icol,ilyr,ibnd) = tauliq_del + tauice_del
                ssaliq_del = ssaliq(ilyr, ibnd) * &
                  (1._wp - forwliq) / &
                  (1._wp - forwliq * ssaliq(ilyr, ibnd))
                ssaice_del = ssaice(ilyr, ibnd) * &
                  (1._wp - forwice) / &
                  (1._wp - forwice * ssaice(ilyr, ibnd))
                scatliq = ssaliq_del * tauliq_del
                scatice = ssaice_del * tauice_del
                this%ssa(icol,ilyr,ibnd) = (scatliq + scatice) / &
                  this%tau(icol,ilyr,ibnd)
                asycld = &
                  (scatliq * (asyliq(ilyr, ibnd) - forwliq) / &
                  (1._wp - forwliq) + &
                  scatice * (asyice(ilyr,ibnd) - forwice) / &
                  (1._wp - forwice)) / (scatliq + scatice)
                this%phi(1,icol,ilyr,ibnd) = 1.0_wp
                this%phi(2,icol,ilyr,ibnd) = asycld
                this%phi(3,icol,ilyr,ibnd) = 0.0_wp
            end select
          endif ! LW/SW
        enddo ! band loop
      enddo ! End layer loop
    enddo ! End column loop
  end function calc_optical_properties

  function get_irad(this, rad, is_liquid, param) result(get_irad)
    real(wp), intent(in) :: rad             ! particle radius
    logical, intent(in) :: is_liquid
    character(len=3), intent(in) :: param   ! ext/ssa/asy
    integer :: get_irad                     ! irad index

    ! Local variables
    integer :: irad, nrad
    real(wp), dimension(4) :: sizreg

    if (phase .eq. 'liq') nrad = 2 else nrad = 3

    if (is_liquid) then
      if (param .eq. 'ext') sizreg = this%r_eff_bounds_liq_ext(:)
      if (param .eq. 'ssa') sizreg = this%r_eff_bounds_liq_ssa(:)
      if (param .eq. 'asy') sizreg = this%r_eff_bounds_liq_ext(:)
    else
      if (param .eq. 'ext') sizreg = this%r_eff_bounds_ice_ext(:)
      if (param .eq. 'ssa') sizreg = this%r_eff_bounds_ice_ssa(:)
      if (param .eq. 'asy') sizreg = this%r_eff_bounds_ice_ext(:)
      ! roughness conditional? r_eff_bounds_ice_* are identical...
    endif ! phase

    do irad = 1, nrad
      if (rad .gt. sizreg(irad) .and. rad .le. sizreg(irad+1)) &
        get_irad = irad
    enddo

  end function get_irad

  function pade(this, coeff, rEff, n_num, n_den) result(out_pade)
    ! Pade formulation: Extinction Coefficient 
    ! (Hogan and Bozzo, ECMWF, TM787, 2016)

    integer, intent(in) :: n_num, n_den
    integer :: nCoeff = n_num + n_den
    real(wp), intent(in) :: coeff(nCoeff)

    integer :: idx
    real(wp) :: numCoeff(n_num), denCoeff(n_den)
    real(wp) :: num = 0, den = 1, coeffSum = 0, rEffSum = 0
    real(wp) :: out_pade

    do idx = 1, n_num
      coeffSum = coeff(idx) * rEffSum
      num = num + coeffSum
      rEffSum = rEffSum + rEff
    enddo

    coeffSum = 0
    rEffSum = rEff
    do idx = n_num, nCoeff
      coeffSum = coeff(idx) * rEffSum
      den = den + coeffSum
      rEffSum = rEffSum + rEff
    enddo

    ! need to validate my calculation!
    out_pade = num / den

  end function pade()

  subroutine stop_on_err(this, msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg

    if(msg /= "") then
      write(error_unit, *) msg
      stop
    end if
  end subroutine

  function get_dim_length(this, ncid, dimname)
    !
    ! Get the length of a dimension from an open netCDF file
    ! This is unfortunately a two-step process
    ! Pernak: I have not modified this function 
    ! (except for aesthetic changes and adding "this" argument)
    !
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer :: get_dim_length

    integer :: dimid

    if(nf90_inq_dimid(ncid, trim(dimname), dimid) == NF90_NOERR) then
      if(nf90_inquire_dimension(&
        ncid, dimid, len=get_dim_length) /= NF90_NOERR) &
        get_dim_length = 0
    else
      get_dim_length = 0
    end if

  end function get_dim_length

  function get_data_size(this, ncid, varName, n)
    !
    ! Returns the extents of a netcdf variable on disk
    ! Pernak: I have not modified this function 
    ! (except for aesthetic changes and adding "this" argument)
    !
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: n
    integer                      :: get_data_size(n)

    integer :: i
    integer :: varid, ndims, dimids(n)

    get_data_size(n) = -1
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("get_data_size: can't find variable " // &
      trim(varName))

    if(&
      nf90_inquire_variable(ncid,varid,ndims=ndims) /= NF90_NOERR) &
      call stop_on_err(&
        "get_data_size: can't get information for variable " // &
        trim(varName))

    if(ndims /= n) &
      call stop_on_err("get_data_size:  variable " // &
        trim(varName) // " has the wrong number of dimensions" )

    if(nf90_inquire_variable(ncid, varid, dimids = dimids) /= &
      NF90_NOERR) &
      call stop_on_err(&
        "get_data_size: can't read dimension ids for variable " // &
        trim(varName))

    do i = 1, n
      if(nf90_inquire_dimension(ncid, dimids(i), &
        len=get_data_size(i)) /= NF90_NOERR) &
        call stop_on_err(&
          "get_data_size: can't get dimension len for variable " // &
          trim(varName))
    end do

  end function get_data_size

  function read_scalar(this, ncid, varName)
    ! Pernak: I have not modified this function 
    ! (except for aesthetic changes and adding "this" argument)

    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    real(wp)                     :: read_scalar

    integer :: varid

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // &
      trim(varName))

    if(nf90_get_var(ncid, varid, read_scalar)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // &
      trim(varName))

  end function read_scalar

  function read_1d_field(this, ncid, varName, nx)
    ! Pernak: I have not modified this function 
    ! (except for aesthetic changes and adding "this" argument)

    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx
    real(wp), dimension(nx)      :: read_1d_field

    integer :: varid

    if(any(get_data_size(ncid, varName, 1) /= [nx])) &
      call stop_on_err("read_field: variable " // trim(varName) // &
        " size is inconsistent.")

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // &
      trim(varName))

    if(nf90_get_var(ncid, varid, read_1d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // &
      trim(varName))

  end function read_1d_field

  function read_2d_field(this, ncid, varName, nx, ny)
    ! Pernak: I have not modified this function 
    ! (except for aesthetic changes and adding "this" argument)

    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx, ny
    real(wp), dimension(nx, ny)  :: read_2d_field

    integer :: varid
    if(any(get_data_size(ncid, varName, 2) /= [nx, ny])) &
      call stop_on_err("read_field: variable " // trim(varName) // &
        " size is inconsistent.")

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // &
      trim(varName))

    if(nf90_get_var(ncid, varid, read_2d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // &
      trim(varName))

  end function read_2d_field

  function read_3d_field(this, ncid, varName, nx, ny, nz)
    ! Pernak: I have not modified this function 
    ! (except for aesthetic changes and adding "this" argument)

    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx, ny, nz
    real(wp), dimension(nx, ny, nz)  :: read_3d_field

    integer :: varid

    if(any(get_data_size(ncid, varName, 3) /= [nx, ny, nz])) &
      call stop_on_err("read_field: variable " // trim(varName) // &
      " size is inconsistent.")

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // &
      trim(varName))

    if(nf90_get_var(ncid, varid, read_3d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // &
      trim(varName))

  end function read_3d_field

  function read_4d_field(this, ncid, varName, nw, nx, ny, nz)
    ! Pernak: I have not modified this function 
    ! (except for aesthetic changes and adding "this" argument)

    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nw, nx, ny, nz
    real(wp), dimension(nw, nx, ny, nz)  :: read_4d_field
    integer :: varid

    if(any(get_data_size(ncid, varName, 4) /= [nw, nx, ny, nz])) &
      call stop_on_err("read_field: variable " // trim(varName) // &
      " size is inconsistent." )

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // &
      trim(varName))

    if(nf90_get_var(ncid, varid, read_4d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // &
      trim(varName))

  end function read_4d_field

end module mo_cloud_optics

