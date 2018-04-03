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
  integer :: icergh
  integer :: ncol, nlay, nlev, nband, nrough, nsizereg, &
    ncoeff_ext, ncoeff_ssa_asy
  real(wp), dimension(:,:), allocatable :: &
    cld_fraction, cld_liq_water_path, cld_ice_water_path, &
    r_eff_liq, r_eff_ice
  ! ---------------------------------------------------------------

  ! ---------------------------------------------------------------
  ! Methods for entire class
  public :: load_cloud_optics, calc_optical_properties, &
            stop_on_err, get_dim_length, read_field, create_var, &
            dim_exists, var_exists, write_field
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
    real(wp), dimension(:,:,:,:), allocatable :: p
  end type ty_cloud_optics_band

  function init(cld_spec, gas_spec, cld_fraction, &
    play, plev, tlay, tsfc, &
    cld_liq_water_path, cld_ice_water_path, ice_rough, &
    r_eff_liq, r_eff_ice, &
    pade_liq_tau, pade_liq_g, pade_ice_tau, pade_ice_g)

    ! initialization of the cloud_optics class
    ! ---------------------------------------------------------------
    !   cld_spec -- object of class ty_cloud_optics_*_band
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
    !   pade_???_tau, pade_???_g -- Pade coefficients size limits for 
    !     ice and liquid for tau and g regimes (integer vectors; need 
    !     to be the same number of elements?)
    ! ---------------------------------------------------------------

    ! ---------------------------------------------------------------

    ! force declaration of all variables ("good practice")
    implicit none

    ! User input
    ! can i just use a conditional in the declarations for cld_spec?
    ! assuming, of course, that we have some sort of keyword to 
    ! specify whether we want a 1-scalar, 2-stream, or n-stream class
    class(ty_cloud_optics_band), intent(out) :: cld_spec
    class(ty_gas_optics_specification), intent(in) :: gas_spec

    real(wp), dimension(:,:), intent(in) :: play, plev, tlay
    real(wp), dimension(:), intent(in) :: tsfc
    integer, intent(in) :: icergh

    real(wp), intent(in), dimension(:,:), allocatable :: &
      cld_fraction, cld_liq_water_path, cld_ice_water_path, &
      r_eff_liq, r_eff_ice
    ! ---------------------------------------------------------------

    ! ---------------------------------------------------------------
    ! Attribute assignments
    ! dimensions
    cld_spec%ncol = size(play, 1)
    cld_spec%nlay = size(play, 2)
    cld_spec%nlev = cld_spec%nlay + 1
    cld_spec%nband = ! from netCDF
    cld_spec%nrough = ! from netCDF
    cld_spec%nsizereg = ! from netCDF
    cld_spec%ncoeff_ext = ! from netCDF
    cld_spec%ncoeff_ssa_asy = ! from netCDF

    ! arrays
    allocate(cloud_spec%cld_fraction(cld_spec%ncol, cld_spec%nlay))
    allocate(cloud_spec%cld_liq_water_path(&
      cld_spec%ncol, cld_spec%nlay))
    allocate(cloud_spec%cld_ice_water_path(&
      cld_spec%ncol, cld_spec%nlay))
    allocate(cloud_spec%r_eff_liq(cld_spec%ncol, cld_spec%nlay))
    allocate(cloud_spec%r_eff_ice(cld_spec%ncol, cld_spec%nlay))

    cld_spec%cld_fraction = cld_fraction
    cld_spec%cld_liq_water_path = cld_liq_water_path
    cld_spec%cld_ice_water_path = cld_ice_water_path
    cld_spec%r_eff_liq = r_reff_liq
    cld_spec%r_eff_ice = r_eff_ice
    ! ---------------------------------------------------------------

  end function init

  function load_cloud_optics(cld_spec)
    ! from an input netCDF, load cloud optics data into 
    ! class attributes
  end function load_cloud_optics

  ! ---------------------------------------------------------------
  ! Not sure if this guy is necessary anymore...but it was first
  ! requested by Robert along with the load function

  function calc_optical_properties(cld_spec)
    ! from cloud optics object, calculate optical properties 
    ! given the specifications in its attributes
  end function calc_optical_properties
  ! ---------------------------------------------------------------
end module mo_cloud_optics

