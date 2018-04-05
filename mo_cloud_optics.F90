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

  ! cloud physical properties
  integer :: roughness
  real(wp), dimension(:,:), allocatable :: &
    cld_fraction, cld_liq_water_path, cld_ice_water_path, &
    r_eff_liq, r_eff_ice

  ! metadata
  character(len=2) :: spectral_domain
  character(len=4) :: approx
  character(len=*) :: cld_coeff_file
  logical :: do_lw, do_pade

  ! Pade coefficients/LUT data
  ! not sure how to handle the two different allocations...LUT 
  ! tables have 1 less dimension for both phases
  real(wp), dimension(:,:,:), allocatable :: &
    co_extliq, co_ssaliq, co_asyliq
  real(wp), dimension(:,:,:,:), allocatable :: &
    co_extice, co_ssaic, co_asyice

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
    allocate(this%cld_liq_water_path(this%ncol, this%nlay))
    allocate(this%cld_ice_water_path(this%ncol, this%nlay))
    allocate(this%r_eff_liq(this%ncol, this%nlay))
    allocate(this%r_eff_ice(this%ncol, this%nlay))
    this%roughness = ice_rough

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
      allocate(this%co_extliq(nCoeffExt, nSizeReg, nBand))
      allocate(this%co_ssaliq(nCoeffSSA, nSizeReg, nBand))
      allocate(this%co_asyliq(nCoeffAsy, nSizeReg, nBand))
      allocate(this%co_extice(nRough, nCoeffExt, nSizeReg, nBand))
      allocate(this%co_ssaice(nRough, nCoeffSSA, nSizeReg, nBand))
      allocate(this%co_asyice(nRough, nCoeffAsy, nSizeReg, nBand))

      ! liquid
      this%co_extliq  = read_field(&
        ncid, extLiqStr, nCoeffExt, nSizeReg, nBand)
      this%co_ssaliq  = read_field(&
        ncid, ssaLiqStr, nCoeffSSA, nSizeReg, nBand)
      this%co_asyliq  = read_field(&
        ncid, asyLiqStr, nCoeffAsy, nSizeReg, nBand)

      ! ice
      this%co_extice  = read_field(&
        ncid, extIceStr, nRough, nCoeffExt, nSizeReg, nBand)
      this%co_ssaice  = read_field(&
        ncid, ssaIceStr, nRough, nCoeffSSA, nSizeReg, nBand)
      this%co_asyice  = read_field(&
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
      allocate(this%co_extliq(nSizeLiq, nBand))
      allocate(this%co_ssaliq(nSizeLiq, nBand))
      allocate(this%co_asyliq(nSizeLiq, nBand))
      allocate(this%co_extice(nRough, nBand, nSizeIce))
      allocate(this%co_ssaice(nRough, nBand, nSizeIce))
      allocate(this%co_asyice(nRough, nBand, nSizeIce))

      ! liquid
      this%co_extliq  = read_field(ncid, extLiqStr, nSizeLiq, nBand)
      this%co_ssaliq  = read_field(ncid, ssaLiqStr, nSizeLiq, nBand)
      this%co_asyliq  = read_field(ncid, asyLiqStr, nSizeLiq, nBand)

      ! ice
      this%co_extice  = read_field(&
        ncid, extIceStr, nRough, nBand, nSizeIce)
      this%co_ssaice  = read_field(&
        ncid, ssaIceStr, nRough, nBand, nSizeIce)
      this%co_asyice  = read_field(&
        ncid, asyIceStr, nRough, nBand, nSizeIce)
    endif ! Pade/LUT

    ncid = nf90_close(ncid)
  end function load_cloud_optics

  function calc_optical_properties(this)
    ! from cloud optics object, calculate optical properties:
    !       tau/OD/extinction/ext, 
    !       omega/single scatter albedo/ssa, 
    !       g/asymmetry/asy
    ! given the specifications in its attributes

    ! ------- Input -------

    ! cloud specification data
    class(ty_optics_1scl), intent(inout) :: this

    ! ------- Output -------
    ! Dimensions: (ncol x nlayers x nbndlw)
    class(ty_cloud_optical_props_arry), intent(inout) :: &
      cloud_optical_props

    ! ------- Local -------

    character(len=128) :: error_msg

    ! Local pade coefficients for extinction and (ssa, g)
    real(wp) :: p_e(6), p(5)

    ! cloud water path (liquid + ice), cloud liquid droplet 
    ! radius (microns), and cloud ice effective size (microns)
    real(wp) :: cwp, radliq, radice, factor, fint

    ! minimum value for cloud quantities
    real(wp), parameter :: cldmin = 1.e-20     

    integer :: index, icol, ilyr, ibnd, irad, irade, irads, iradg

    ! ice surface roughness (1 = none, 2 = medium, 3 = high)
    integer :: icergh                          

    ! liquid extinction coefficient, single scattering albedo, and 
    ! asymmetry parameter
    real(wp) :: extliq(nlayers, nbnd), ssaliq(nlayers, nbnd), &
                asyliq(nlayers, nbnd)

    ! ice extinction coefficients, single scattering albedo, and 
    ! asymmetry parameter
    real(wp) :: extice(nlayers, nbnd), ssaice(nlayers, nbnd), &
                asyice(nlayers, nbnd)

    ! liquid and ice cloud extinction optical depth, no delta scaling
    real(wp) :: tauliq, tauice

    ! ice and liquid scattering terms, (liquid + ice?) asymmetry
    real(wp) :: scatice, scatliq, asycld

    ! Forward scattering terms
    real(wp) :: forwliq, forwice

    ! delta scaling extinction OD and single scattering albedo
    real(wp) :: tauliq_del, tauice_del, ssaliq_del, ssaice_del

    ! ------- Definitions -------

    ! ------- Error checking -------
    error_msg = ''
    icergh = this%icergh
    if (icergh < 1 .or. icergh > 3) then
      error_msg = 'cloud optics: ' \\
        'cloud ice surface roughness flag is out of bounds'
      return
    endif

    ! RP comments: 
    ! As per all other RRTMGP examples, the validity of input data 
    ! should be checked the outside computational loop. That means 
    ! checking that sizes are in range, values are positive, etc. 
    !
    ! One could also compute a logical "cloud present" mask in this 
    ! loop and use it below.  
    !   
    ! No magic numbers. When checking validity of e.g. size refer to 
    ! the size regime arrays. 

    ! RP comments
    ! Some (many?) of these values will be replaced. Better to use 
    ! merge() in loops below. 

    ! Initialize
    extliq(:,:) = 0.0_wp
    ssaliq(:,:) = 1.0_wp
    asyliq(:,:) = 0.0_wp
    extice(:,:) = 0.0_wp
    ssaice(:,:) = 1.0_wp
    asyice(:,:) = 0.0_wp

    select type(cloud_optical_props)
      type is (ty_cloud_optical_props_1scl)
      cloud_optical_props%taucld(:,:,:) = 0.0_wp
      type is (ty_cloud_optical_props_2str)
        cloud_optical_props%taucld(:,:,:) = 0.0_wp
        cloud_optical_props%ssacld(:,:,:) = 1.0_wp
        cloud_optical_props%asycld(:,:,:) = 0.0_wp
      type is (ty_cloud_optical_props_nstr)
        cloud_optical_props%taucld(:,:,:) = 0.0_wp
        cloud_optical_props%ssacld(:,:,:) = 1.0_wp
        cloud_optical_props%pcld(:,:,:,:) = 0.0_wp
    end select

    ! Main column loop
    do icol = 1, this%ncol

    ! Cloud optical properties from Pade coefficients

      do ilyr = 1, this%nlayers
        cwp = this%ciwp(icol,ilyr) + this%clwp(icol,ilyr)
        if (this%cldfrac(icol,ilyr) .gt. cldmin .and. &
            cwp .gt. cldmin) then

          ! Derive optical properties from Pade functions 
          ! Original Pade formulations for EXT, SSA, G

          ! Liquid
          radliq = this%rel(icol,ilyr)
          if (radliq .gt. 0.0_wp .and. &
              this%clwp(icol,ilyr) .gt. 0.0_wp) then 

            ! For liquid OP, particle size is limited to 2.5-20.0 
            ! microns (LW only?)
            if (radliq .lt. 2.5_wp .or. radliq .gt. 20.0_wp) then 
              error_msg = 'cloud optics: ' \\
                'liquid effective radius is out of bounds'
              return
            endif

            ! Define coefficient particle size regime for 
            ! current size: extinction, ssa
            irade = this%get_irad(&
              radliq, 'liq', this%is_lw, 'ext')
            irads = this%get_irad(&
              radliq, 'liq', this%is_lw, 'ssa')
            iradg = this%get_irad(&
              radliq, 'liq', this%is_lw, 'asy')

            do ibnd = 1, this%nband
              extliq(ilyr,ibnd) = this%pade_ext(&
                this%pade_extliq(ibnd,irade,:), radliq)
              ssaliq(ilyr,ibnd) = this%pade_ssa(&
                this%pade_ssaliq(ibnd,irads,:), radliq)
              asyliq(ilyr,ibnd) = this%pade_asy(&
                this%pade_asyliq(ibnd,iradg,:), radliq)
            enddo
          endif ! radliq provision

          ! Ice
          radice = this%rei(icol,ilyr)
          if (radice .gt. 0.0_wp .and. &
              this%ciwp(icol,ilyr) .gt. 0.0_wp) then 

            ! For Yang (2013) ice OP, particle size is limited to 
            ! 10.0-180.0 microns (LW only?)
            if (radice .lt. 10.0_wp .or. radice .gt. 180.0_wp) then
              error_msg = &
                'cloud optics: ice effective radius is out of bounds'
              return
            endif

            ! Define coefficient particle size regime for current 
            ! size: extinction, ssa
            irade = this%get_irad(&
              radice, 'ice', this&is_lw, 'ext')
            irads = this%get_irad(&
              radice, 'ice', this&is_lw, 'ssa')
            iradg = this%get_irad(&
              radice, 'ice', this&is_lw, 'asy')

            do ibnd = 1, nband
              ! Derive optical properties for selected ice roughness
              extice(ilyr,ibnd) = this%pade_ext(&
                this%pade_extice(ibnd,irade,:,icergh), radice)
              ssaice(ilyr,ibnd) = this%pade_ssa(&
                this%pade_ssaice(ibnd,irads,:,icergh), radice)
              asyice(ilyr,ibnd) = this%pade_asy(&
                this%pade_asyice(ibnd,iradg,:,icergh), radice)
            enddo ! end band loop
          endif ! radice
         endif ! cldfrac and cwp provision
      enddo ! End layer loop
    enddo ! End column loop

    ! Combine liquid and ice contributions into total cloud 
    ! optical properties
    do icol = 1, this%ncol
      do ilyr = 1, this%nlayers
        cwp = ciwp(icol,ilyr) + clwp(icol,ilyr)
        if (this%cldfrac(icol,ilyr) .gt. cldmin .and. &
            cwp .gt. cldmin) then 

          do ibnd = 1, this%nbnd
            tauice = ciwp(icol,ilyr) * extice(ilyr,ibnd)
            tauliq = clwp(icol,ilyr) * extliq(ilyr,ibnd)

            if (tauice == 0.0_wp .and. tauliq == 0.0_wp) then
              error_msg = &
                'cloud optics: cloud optical depth is zero'
              return
            endif

            if this%is_lw .eq. 1 then
              ! LONGWAVE NO DSCALE
              scatice = ssaice(ilyr,ibnd) * tauice
              scatliq = ssaliq(ilyr,ibnd) * tauliq

              select type(cloud_optical_props)
                type is (ty_cloud_optical_props_1scl)
                  cloud_optical_props%taucld(icol,ilyr,ibnd) = &
                    tauice + tauliq
                type is (ty_cloud_optical_props_2str)
                  cloud_optical_props%taucld(icol,ilyr,ibnd) = &
                    tauice + tauliq
                  cloud_optical_props%ssacld(icol,ilyr,ibnd) = &
                    (scatice + scatliq) / &
                    cloud_optical_props%taucld(icol,ilyr,ibnd)
                  cloud_optical_props%asycld(icol,ilyr,ibnd) = &
                    (scatice * asyice(ilyr,ibnd) + &
                    (scatliq * asyliq(ilyr,ibnd))) / &
                    (scatice + scatliq)
                type is (ty_cloud_optical_props_nstr)
                  if (clwp(icol,ilyr) .gt. 0.0_wp) then 
                    error_msg = 'cloud optics: n-stream option ' \\
                      'not available in longwave for liquid clouds'
                    return
                  else
                    cloud_optical_props%taucld(icol,ilyr,ibnd) = &
                      tauice + tauliq
                    cloud_optical_props%ssacld(icol,ilyr,ibnd) = &
                      (scatice + scatliq) / & 
                      cloud_optical_props%taucld(icol,ilyr,ibnd)
                    asycld = &
                      (scatice * asyice(ilyr,ibnd) + &
                      (scatliq * asyliq(ilyr,ibnd))) / &
                      (scatice + scatliq)
                    cloud_optical_props%pcld(1,icol,ilyr,ibnd) = &
                      1.0_wp
                    cloud_optical_props%pcld(2,icol,ilyr,ibnd) = &
                      asycld
                    cloud_optical_props%pcld(3,icol,ilyr,ibnd) = &
                      asycld**2
                  endif ! clwp provision
              end select
            else
              ! SHORTWAVE
              forwliq = asyliq(ilyr,ibnd)**2
              forwice = asyice(ilyr,ibnd)**2
              tauliq_del = (1.0_wp - forwliq * ssaliq(ilyr,ibnd)) * &
                tauliq
              tauice_del = (1.0_wp - forwice * ssaice(ilyr,ibnd)) * &
                tauice

              select type(cloud_optical_props)
                type is (ty_cloud_optical_props_1scl)
                  cloud_optical_props%taucld(icol,ilyr,ibnd) = &
                    tauliq_del + tauice_del
                type is (ty_cloud_optical_props_2str)
                  cloud_optical_props%taucld(icol,ilyr,ibnd) = &
                    tauliq_del + tauice_del
                  ssaliq_del = &
                    ssaliq(ilyr,ibnd) * (1._wp - forwliq) / &
                    (1._wp - forwliq * ssaliq(ilyr,ibnd))
                  ssaice_del = &
                    ssaice(ilyr,ibnd) * (1._wp - forwice) / &
                    (1._wp - forwice * ssaice(ilyr,ibnd))
                  scatliq = ssaliq_del * tauliq_del
                  scatice = ssaice_del * tauice_del
                  cloud_optical_props%ssacld(icol,ilyr,ibnd) = &
                    (scatliq + scatice) / &
                    cloud_optical_props%taucld(icol,ilyr,ibnd)
                  cloud_optical_props%asycld(icol,ilyr,ibnd) = &
                    (scatice * (asyice(ilyr,ibnd) - forwice) / &
                    (1._wp - forwice) + &
                    scatliq * (asyliq(ilyr,ibnd) - forwliq) / &
                    (1._wp - forwliq)) / &
                    (scatice + scatliq)
                type is (ty_cloud_optical_props_nstr)
                  cloud_optical_props%taucld(icol,ilyr,ibnd) = &
                    tauliq_del + tauice_del
                  ssaliq_del = ssaliq(ilyr,ibnd) * &
                               (1._wp - forwliq) / &
                               (1._wp - forwliq * ssaliq(ilyr,ibnd))
                  ssaice_del = ssaice(ilyr,ibnd) * &
                               (1._wp - forwice) / &
                               (1._wp - forwice * ssaice(ilyr,ibnd))
                  scatliq = ssaliq_del * tauliq_del
                  scatice = ssaice_del * tauice_del
                  cloud_optical_props%ssacld(icol,ilyr,ibnd) = &
                    (scatliq + scatice) / &
                     cloud_optical_props%taucld(icol,ilyr,ibnd)
                  asycld = &
                    (scatliq * (asyliq(ilyr,ibnd) - forwliq) / &
                    (1._wp - forwliq) + &
                    scatice * (asyice(ilyr,ibnd) - forwice) / &
                    (1._wp - forwice)) / &
                    (scatliq + scatice)
                  cloud_optical_props%pcld(1,icol,ilyr,ibnd) = 1.0_wp
                  cloud_optical_props%pcld(2,icol,ilyr,ibnd) = asycld
                  cloud_optical_props%pcld(3,icol,ilyr,ibnd) = 0.0_wp
                end select
              else
                select type(cloud_optical_props)
                  type is (ty_cloud_optical_props_1scl)
                    cloud_optical_props%taucld(icol,ilyr,ibnd) = &
                      tauice + tauliq
                  type is (ty_cloud_optical_props_2str)
                    cloud_optical_props%taucld(icol,ilyr,ibnd) = &
                      tauice + tauliq
                    ! Pernak: this is the only substantial difference 
                    ! with LW (no delta-scaling) -- can't it just go 
                    ! outside the select block? how would 2-str even
                    ! work if scatice and scatliq are defined here?
                    scatice = ssaice(ilyr,ibnd) * tauice
                    scatliq = ssaliq(ilyr,ibnd) * tauliq

                    cloud_optical_props%ssacld(icol,ilyr,ibnd) = &
                      (scatice + scatliq) / &
                       cloud_optical_props%taucld(icol,ilyr,ibnd)
                    cloud_optical_props%asycld(icol,ilyr,ibnd) = &
                      (scatice * asyice(ilyr,ibnd) + &
                      (scatliq * asyliq(ilyr,ibnd))) / &
                      (scatice + scatliq)
                  type is (ty_cloud_optical_props_nstr)
                    cloud_optical_props%taucld(icol,ilyr,ibnd) = &
                      tauice + tauliq
                    cloud_optical_props%ssacld(icol,ilyr,ibnd) = &
                      (scatice + scatliq) / &
                      cloud_optical_props%taucld(icol,ilyr,ibnd)
                    asycld = &
                      (scatice * asyice(ilyr,ibnd) + &
                      (scatliq * asyliq(ilyr,ibnd))) / &
                      (scatice + scatliq)
                    cloud_optical_props%pcld(1,icol,ilyr,ibnd) = &
                      1.0_wp
                    cloud_optical_props%pcld(2,icol,ilyr,ibnd) = &
                      asycld
                    cloud_optical_props%pcld(3,icol,ilyr,ibnd) = 
                      asycld**2
                end select
              endif ! SW delta scaling
            endif ! LW/SW
          enddo ! band loop
        endif ! cldmin provision
      enddo ! layer loop
    enddo ! column loop
  end function calc_optical_properties
end module mo_cloud_optics

