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

  interface write_field 
    module procedure write_1d_int_field, write_2d_int_field, &
      write_1d_field, write_2d_field, write_3d_field, write_4d_field
  end interface 

  ! ---------------------------------------------------------------
  ! Attributes for entire class

  ! dimensions
  integer :: ncol, nlay, nlev, nband, nrough, nsizereg, &
    nsizeice, nsizeliq, ncoeff_ext, ncoeff_ssa, ncoeff_asy

  ! optical properties
  integer :: icergh
  real(wp), dimension(:,:), allocatable :: &
    cld_fraction, cld_liq_water_path, cld_ice_water_path, &
    r_eff_liq, r_eff_ice
  character(len=2) :: spectral_domain
  character(len=4) :: approx
  character(len=*) :: cld_coeff_file
  logical :: do_lw, do_pade

  ! Pade coefficients/LUT data
  ! not sure how to handle the two different allocations...LUT 
  ! tables have 1 less dimension for both phases
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
  public :: load_cloud_optics, calc_optical_properties, &
    stop_on_err, get_dim_length, create_var, dim_exists, var_exists, &
    read_field, write_field
  ! ---------------------------------------------------------------

  type, extends(ty_spectral_disc), public :: ty_cloud_optics_band_1scl
    ! Attributes in extended (sub)class
    real(wp), dimension(:,:,:), allocatable :: optical_depth
  end type ty_cloud_optics_band

  type, extends(ty_spectral_disc), public :: ty_cloud_optics_band_2str
    ! Attributes in extended (sub)class
    real(wp), dimension(:,:,:), allocatable :: &
      optical_depth, ssalbedo, asymmetry
  end type ty_cloud_optics_band

  type, extends(ty_spectral_disc), public :: ty_cloud_optics_band_nstr
    ! Attributes in extended (sub)class
    real(wp), dimension(:,:,:), allocatable :: optical_depth, ssalbedo
    real(wp), dimension(:,:,:,:), allocatable :: phase
  end type ty_cloud_optics_band

  contains

  function init(this, gas_spec, cld_fraction, &
    play, plev, tlay, tsfc, &
    cld_liq_water_path, cld_ice_water_path, ice_rough, &
    r_eff_liq, r_eff_ice, method, is_lw, is_pade)

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
    ! ---------------------------------------------------------------

    ! ---------------------------------------------------------------

    ! force declaration of all variables ("good practice")
    implicit none

    ! ---------------------------------------------------------------
    ! User input
    class(ty_gas_optics_specification), intent(in) :: gas_spec

    real(wp), dimension(:,:), intent(in) :: play, plev, tlay
    real(wp), dimension(:), intent(in) :: tsfc
    integer, intent(in) :: icergh

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
      class(ty_cloud_optics_band_nstr), intent(out) :: this
    else
      err_msg = "Please specify '1scl', '2str', or 'nstr' for method"
      call stop_on_err(err_msg)
    endif
    ! ---------------------------------------------------------------

    ! ---------------------------------------------------------------
    ! Intermediates (neither input nor output)
    character(len=*) :: cld_coeff_file, domain_str, approx_str, &
      extliq_str, ssaliq_str, asyliq_str, &
      extice_str, ssaice_str, asyice_str

    integer, dimension(4) :: &
      r_eff_bounds_liq_ext, r_eff_bounds_ice_ext, &
      r_eff_bounds_liq_ssa, r_eff_bounds_ice_ssa, &
      r_eff_bounds_liq_asy, r_eff_bounds_ice_asy 
    ! ---------------------------------------------------------------

    ! ---------------------------------------------------------------
    ! Attribute assignments
    ! dimensions
    this%ncol = size(play, 1)
    this%nlay = size(play, 2)
    this%nlev = size(plev, 2)

    ! "metadata" properties
    this%do_lw = is_lw
    this%do_pade = is_pade

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

    this%cld_coeff_file = 'rrtmgp-' // &
      this%spectral_domain // '-inputs-cloud-optics-' // &
      this%approx //'.nc'

    ! physical property arrays
    allocate(this%cld_fraction(this%ncol, this%nlay))
    allocate(this%cld_liq_water_path(&
      this%ncol, this%nlay))
    allocate(cloud_spec%cld_ice_water_path(&
      this%ncol, this%nlay))
    allocate(this%r_eff_liq(this%ncol, this%nlay))
    allocate(this%r_eff_ice(this%ncol, this%nlay))

    this%cld_fraction = cld_fraction
    this%cld_liq_water_path = cld_liq_water_path
    this%cld_ice_water_path = cld_ice_water_path
    this%r_eff_liq = r_reff_liq
    this%r_eff_ice = r_eff_ice

    ! load LUT or Pade Coefficients
    load_cloud_optics(this)

    ! compute optical properties
    calc_optical_properties(this)
    ! ---------------------------------------------------------------
  end function init

  function load_cloud_optics(this)
    ! from an input netCDF, load cloud optics data into 
    ! class attributes

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
      allocate(this%extliq(nCoeffExt, nSizeReg, nBand))
      allocate(this%ssaliq(nCoeffSSA, nSizeReg, nBand))
      allocate(this%asyliq(nCoeffAsy, nSizeReg, nBand))
      allocate(this%extice(nRough, nCoeffExt, nSizeReg, nBand))
      allocate(this%ssaice(nRough, nCoeffSSA, nSizeReg, nBand))
      allocate(this%asyice(nRough, nCoeffAsy, nSizeReg, nBand))

      ! liquid
      this%extliq  = read_field(&
        ncid, extLiqStr, nCoeffExt, nSizeReg, nBand)
      this%ssaliq  = read_field(&
        ncid, ssaLiqStr, nCoeffSSA, nSizeReg, nBand)
      this%asyliq  = read_field(&
        ncid, asyLiqStr, nCoeffAsy, nSizeReg, nBand)

      ! ice
      this%extice  = read_field(&
        ncid, extIceStr, nRough, nCoeffExt, nSizeReg, nBand)
      this%ssaice  = read_field(&
        ncid, ssaIceStr, nRough, nCoeffSSA, nSizeReg, nBand)
      this%asyice  = read_field(&
        ncid, asyIceStr, nRough, nCoeffAsy, nSizeReg, nBand)

      ! Particle size regimes for Pade formulations
      ! these bounds were determined by fitting the Pade polynomials 
      ! to the lookup table data (the fits had to be done in groups of 
      ! of effective radius)
      ! SW and LW are mostly the same (defaults are LW)
      this%r_eff_bounds_liq_ext = (/2,10,35,60/)
      this%r_eff_bounds_liq_ssa = (/2,10,35,60/)
      this%r_eff_bounds_liq_asy = (/2,9,20,60/)
      this%r_eff_bounds_ice_ext = (/10,20,30,180/)
      this%r_eff_bounds_ice_ssa = (/10,20,30,180/)
      this%r_eff_bounds_ice_asy = (/10,20,30,180/)

      if (this%do_lw) then
        r_eff_bounds_liq_asy(2) = 8
        r_eff_bounds_ice_asy(2) = 50
      endif ! LW
    else ! LUT
      ! LUT-exclusive dimensions and shortcuts
      this%nsizeice = get_dim_length(ncid, 'nsize_ice')
      this%nsizeliq = get_dim_length(ncid, 'nsize_liq')
      nSizeIce = this%nsizeice
      nSizeLiq = this%nsizeliq

      ! LUT data input arrays
      allocate(this%extliq(nSizeLiq, nBand))
      allocate(this%ssaliq(nSizeLiq, nBand))
      allocate(this%asyliq(nSizeLiq, nBand))
      allocate(this%extice(nRough, nBand, nSizeIce))
      allocate(this%ssaice(nRough, nBand, nSizeIce))
      allocate(this%asyice(nRough, nBand, nSizeIce))

      ! liquid
      this%extliq  = read_field(ncid, extLiqStr, nSizeLiq, nBand)
      this%ssaliq  = read_field(ncid, ssaLiqStr, nSizeLiq, nBand)
      this%asyliq  = read_field(ncid, asyLiqStr, nSizeLiq, nBand)

      ! ice
      this%extice  = read_field(&
        ncid, extIceStr, nRough, nBand, nSizeIce)
      this%ssaice  = read_field(&
        ncid, ssaIceStr, nRough, nBand, nSizeIce)
      this%asyice  = read_field(&
        ncid, asyIceStr, nRough, nBand, nSizeIce)
    endif ! Pade/LUT

    ncid = nf90_close(ncid) 
  end function load_cloud_optics

  function calc_optical_properties(this)
    ! from cloud optics object, calculate optical properties 
    ! (tau,OD/extinction/ext, omega/single scatter albedo/ssa, 
    ! g/asymmetry/asy) given the specifications in its attributes
  end function calc_optical_properties
end module mo_cloud_optics

