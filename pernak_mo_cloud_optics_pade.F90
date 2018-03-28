! Module: mo_cloud_optics_pade
!
! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
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

module mo_cloud_optics_pade
  use mo_rte_kind,            only: wp
  ! Pernak: can we just use nscl? I don't see how 1scl, 2scl, or 
  ! nscl are differentiated
  use mo_optical_props,       only: ty_optical_props_1scl, &
                                    ty_optical_props_2scl
                                    ty_optical_props_nscl
  use mo_cloud_optical_props, only: ty_cloud_optical_props_arry, &
                                    ty_cloud_optical_props_1scl, &
                                    ty_cloud_optical_props_2str, &
                                    ty_cloud_optical_props_nstr
  use netcdf

  implicit none
  private
  public :: is_lw
  public :: stop_on_err, get_dim_length, read_field, create_var, &
            dim_exists, var_exists, write_field

  interface read_field
    module procedure read_scalar, read_1d_field, read_2d_field, &
                     read_3d_field, read_4d_field
  end interface

  interface write_field 
    module procedure write_1d_int_field, write_2d_int_field, &
                     write_1d_field, write_2d_field, &
                     write_3d_field, write_4d_field
  end interface 

!--------------------------------------------------------------------
!
! RP comments: 
!  It's true that ty_gas_optics extends ty_spectral_disc, but this is 
!  unique because the gas optics actually sets the spectral
!  discretization. Instead we should extend ty_optical_props_1scl,
!  _2str, and _nstr from mo_optical_props. I don't see a way to avoid
!  extending each with the same data components.
!
!  validation() will need to be extended to ensure that the cloud 
!  optics data are present,  that cloud physical and optical 
!  properties arrays are allocated, the same size as the tau/omega0/g
!  arrays, ... 
!
!  This code can be used generally -- it's just optical properties by 
!  band from data -- so maybe call it "ty_optics_1scl" etc.? 
!
!  Pernak: for now, i only have a ty_optics_1scl type extension, and 
!  the others would just be a copy and paste with changes to the type
!  name. until i understand this better, we will not have the 
!  flexibility to extend to 2- and n-stream
!--------------------------------------------------------------------

  type, extends(ty_optical_props_1scl), public :: ty_optics_1scl

    ! User input
    ! Ice surface roughness category - needed for Yang (2013) ice 
    ! optics parameterization (1 = none, 2 = medium, 3 = high)
    integer :: icergh

    ! Cloud physical properties (dimensions, cloud fraction, cloud ice
    ! water path, cloud liquid water path, cloud ice effective radius
    ! (microns), cloud liquid effective radius (microns) )
    integer :: ncol, nlay
    real(wp), dimension(:,:), allocatable :: cldfrac
    real(wp), dimension(:,:), allocatable :: ciwp
    real(wp), dimension(:,:), allocatable :: clwp
    real(wp), dimension(:,:), allocatable :: rei
    real(wp), dimension(:,:), allocatable :: rel

    !-----------------------------------------------------------------
    ! RP comments: 
    !  Cloud optics types shouldn't have hardcoded LW/SW distinctions. 
    !  Better to have explicit wavelength or wavenumber limits on 
    !  bands provided at initialization. 
    !-----------------------------------------------------------------

    ! Model input: LUT
    ! Pade coefficient dimensions
    integer :: nband, nrghice, nsizereg, ncoeff_ext, ncoeff_ssa_g

    ! RP comments: 
    !  Magic numbers are undesirable
    !  Size regimes should be parameterized here: number of regimes, 
    !  size limits.
    !  Pernak response: not entirely sure how to do this -- via
    !    allocation?

    ! Particle size regimes for Pade formulations
    !integer, dimension(2,4) :: pade_sizreg_liq, pade_sizreg_ice
    !integer, dimension(:,:) :: pade_sizreg_liq, pade_sizreg_ice
    integer, dimension(:), allocatable :: pade_sizreg_liq_g, &
                                          pade_sizreg_liq_tau
    integer, dimension(:), allocatable :: pade_sizreg_ice_g
                                          pade_sizreg_ice_tau

    ! Pade cloud coefficients
    real(wp), dimension(:,:,:  ), allocatable :: pade_extliq
    real(wp), dimension(:,:,:  ), allocatable :: pade_ssaliq
    real(wp), dimension(:,:,:  ), allocatable :: pade_asyliq
    real(wp), dimension(:,:,:,:), allocatable :: pade_extice
    real(wp), dimension(:,:,:,:), allocatable :: pade_ssaice
    real(wp), dimension(:,:,:,:), allocatable :: pade_asyice

    contains

    generic, public :: init_cldopt  => init_cldopt
    generic, public :: cloud_optics => cloud_optics

    ! Internal procedures
    procedure, private :: init_cldopt
    procedure, private :: cloud_optics
    procedure, private :: get_irad
    procedure, private :: pade_ext
    procedure, private :: pade_ssa
    procedure, private :: pade_asy

  end type ty_optics_1scl

  ! mo_cloud_optics_pade methods
  contains

  !-------------------------------------------------------------------
  !
  ! Cloud optics initialization 
  !
  !-------------------------------------------------------------------

  ! Initialize object based on data read from netCDF file however the 
  ! user desires. Pade coefficient data used for cloud 
  ! property conversion
  !-------------------------------------------------------------------

  function init_cldopt(cloud_spec, cld_coeff_file, is_lw) &
    result(error_msg)

    ! Pade cloud optics
    class(ty_optics_1scl), intent(inout) :: cloud_spec

    ! Cloud coefficient optical property input file - Pade
    character(len=*), intent(in) :: cld_coeff_file

    logical, intent(in) :: is_lw
    character(len = 128) error_msg

    ! Local variables
    integer :: ncid, nband, nrghice, nsizereg, &
      ncoeff_ext, ncoeff_ssa_g
    integer, dimension(:), allocatable :: pade_sizreg_liq_g, &
                                          pade_sizreg_liq_tau
    integer, dimension(:), allocatable :: pade_sizreg_ice_g
                                          pade_sizreg_ice_tau

    ! Particle size regimes for Pade formulations
    cloud_spec%pade_sizreg_liq_g(:) = (/2,10,35,60/)
    cloud_spec%pade_sizreg_ice_g(:) = (/10,20,30,180/)

    if (is_lw) then
      cloud_spec%pade_sizreg_liq_tau(:) = (/2,8,20,60/)
      cloud_spec%pade_sizreg_ice_tau(:) = (/10,20,50,180/)
    else
      cloud_spec%pade_sizreg_liq_tau(:) = (/2,9,20,60/)
      cloud_spec%pade_sizreg_ice_tau(:) = (/10,20,30,180/)
    endif

    ! RP comments: 
    !  Types should be initialzed with data directly. 
    !  No I/O in code we will distribute. 
    !  See the initialization of gas_optics
    !  Pernak: doing this here with some other attribute assignments
    !cloud_spec%ncol = ...
    !cloud_spec%nlayers = ...
    !cloud_spec%nband = nband ! defined later... 
    !cloud_spec%is_lw = is_lw
    !cloud_spec%cldfrac = ... ! cloud fraction (ncol x nlayers)
    !cloud_spec%clwp = ... ! cloud liquid water path (ncol x nlayers)
    !cloud_spec%ciwp = ... ! cloud ice water path (ncol x nlayers)
    !cloud_spec%rel = ... ! cloud liquid particle effective 
                          ! radius (microns), (ncol x nlayers)
    !cloud_spec%rei = ... ! cloud ice particle effective 
                          ! radius (microns), (ncol x nlayers)

    error_msg = ""

    ! Open cloud optical property coefficient file
    if (nf90_open(trim(cld_coeff_file), NF90_WRITE, ncid) &
        /= NF90_NOERR) & call stop_on_err(&
        "init_cldopt: can't open file " // trim(cld_coeff_file))

    ! read in array dimensions
    nband        = get_dim_length(ncid, 'nband') 
    nrghice      = get_dim_length(ncid, 'nrghice') 
    nsizereg     = get_dim_length(ncid, 'nsizereg') 
    ncoeff_ext   = get_dim_length(ncid, 'ncoeff_ext') 
    ncoeff_ssa_g = get_dim_length(ncid, 'ncoeff_ssa_g') 

    ! create new object attributes and populate array dimensions
    cloud_spec%nband        = nband
    cloud_spec%nrghice      = nrghice
    cloud_spec%nsizereg     = nsizereg
    cloud_spec%ncoeff_ext   = ncoeff_ext
    cloud_spec%ncoeff_ssa_g = ncoeff_ssa_g

    ! Allocate cloud property Pade coefficient input arrays
    allocate(cloud_spec%pade_extliq(nband, nsizereg, ncoeff_ext))
    allocate(cloud_spec%pade_ssaliq(nband, nsizereg, ncoeff_ssa_g))
    allocate(cloud_spec%pade_asyliq(nband, nsizereg, ncoeff_ssa_g))
    allocate(cloud_spec%pade_extice(nband, nsizereg, ncoeff_ext, &
      nrghice))
    allocate(cloud_spec%pade_ssaice(nband, nsizereg, ncoeff_ssa_g, &
      nrghice))
    allocate(cloud_spec%pade_asyice(nband, nsizereg, ncoeff_ssa_g, &
      nrghice))

    ! more object attributes: parameter arrays
    cloud_spec%pade_extliq  = read_field(ncid, 'pade_extliq', &
      nband, nsizereg, ncoeff_ext) 
    cloud_spec%pade_ssaliq  = read_field(ncid, 'pade_ssaliq', &
      nband, nsizereg, ncoeff_ssa_g) 
    cloud_spec%pade_asyliq  = read_field(ncid, 'pade_asyliq', &
      nband, nsizereg, ncoeff_ssa_g) 
    cloud_spec%pade_extice  = read_field(ncid, 'pade_extice', &
      nband, nsizereg, ncoeff_ext, nrghice) 
    cloud_spec%pade_ssaice  = read_field(ncid, 'pade_ssaice', &
      nband, nsizereg, ncoeff_ssa_g, nrghice) 
    cloud_spec%pade_asyice  = read_field(ncid, 'pade_asyice', &
      nband, nsizereg, ncoeff_ssa_g, nrghice) 

    ncid = nf90_close(ncid) 

  end function init_cldopt

  ! -----------------------------------------------------------------
  !
  ! Derive cloud optical properties from provided cloud 
  ! physical pproperties
  !
  ! -----------------------------------------------------------------

  function cloud_optics(cloud_spec, cloud_optical_props) &
    result(error_msg)
    ! ---------------------------------------------------------------
    ! Purpose:  Compute the cloud optical properties for each 
    ! cloudy layer.

    ! ------- Input -------

    ! cloud specification data
    class(ty_optics_1scl), intent(inout) :: cloud_spec

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
    icergh = cloud_spec%icergh
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
    do icol = 1, cloud_spec%ncol

    ! Cloud optical properties from Pade coefficients

      do ilyr = 1, cloud_spec%nlayers
        cwp = cloud_spec%ciwp(icol,ilyr) + cloud_spec%clwp(icol,ilyr)
        if (cloud_spec%cldfrac(icol,ilyr) .gt. cldmin .and. &
            cwp .gt. cldmin) then

          ! Derive optical properties from Pade functions 
          ! Original Pade formulations for EXT, SSA, G

          ! Liquid
          radliq = cloud_spec%rel(icol,ilyr)
          if (radliq .gt. 0.0_wp .and. &
              cloud_spec%clwp(icol,ilyr) .gt. 0.0_wp) then 

            ! For liquid OP, particle size is limited to 2.5-20.0 
            ! microns (LW only?)
            if (radliq .lt. 2.5_wp .or. radliq .gt. 20.0_wp) then 
              error_msg = 'cloud optics: ' \\
                'liquid effective radius is out of bounds'
              return
            endif

            ! Define coefficient particle size regime for 
            ! current size: extinction, ssa
            irade = cloud_spec%get_irad(&
              radliq, 'liq', cloud_spec%is_lw, 'ext')
            irads = cloud_spec%get_irad(&
              radliq, 'liq', cloud_spec%is_lw, 'ssa')
            iradg = cloud_spec%get_irad(&
              radliq, 'liq', cloud_spec%is_lw, 'asy')

            do ibnd = 1, cloud_spec%nband
              extliq(ilyr,ibnd) = cloud_spec%pade_ext(&
                cloud_spec%pade_extliq(ibnd,irade,:), radliq)
              ssaliq(ilyr,ibnd) = cloud_spec%pade_ssa(&
                cloud_spec%pade_ssaliq(ibnd,irads,:), radliq)
              asyliq(ilyr,ibnd) = cloud_spec%pade_asy(&
                cloud_spec%pade_asyliq(ibnd,iradg,:), radliq)
            enddo
          endif ! radliq provision

          ! Ice
          radice = cloud_spec%rei(icol,ilyr)
          if (radice .gt. 0.0_wp .and. &
              cloud_spec%ciwp(icol,ilyr) .gt. 0.0_wp) then 

            ! For Yang (2013) ice OP, particle size is limited to 
            ! 10.0-180.0 microns (LW only?)
            if (radice .lt. 10.0_wp .or. radice .gt. 180.0_wp) then
              error_msg = &
                'cloud optics: ice effective radius is out of bounds'
              return
            endif

            ! Define coefficient particle size regime for current 
            ! size: extinction, ssa
            irade = cloud_spec%get_irad(&
              radice, 'ice', cloud_spec&is_lw, 'ext')
            irads = cloud_spec%get_irad(&
              radice, 'ice', cloud_spec&is_lw, 'ssa')
            iradg = cloud_spec%get_irad(&
              radice, 'ice', cloud_spec&is_lw, 'asy')

            do ibnd = 1, nband
              ! Derive optical properties for selected ice roughness
              extice(ilyr,ibnd) = cloud_spec%pade_ext(&
                cloud_spec%pade_extice(ibnd,irade,:,icergh), radice)
              ssaice(ilyr,ibnd) = cloud_spec%pade_ssa(&
                cloud_spec%pade_ssaice(ibnd,irads,:,icergh), radice)
              asyice(ilyr,ibnd) = cloud_spec%pade_asy(&
                cloud_spec%pade_asyice(ibnd,iradg,:,icergh), radice)
            enddo ! end band loop
          endif ! radice
         endif ! cldfrac and cwp provision
      enddo ! End layer loop
    enddo ! End column loop

    ! Combine liquid and ice contributions into total cloud 
    ! optical properties
    do icol = 1, cloud_spec%ncol
      do ilyr = 1, cloud_spec%nlayers
        cwp = ciwp(icol,ilyr) + clwp(icol,ilyr)
        if (cloud_spec%cldfrac(icol,ilyr) .gt. cldmin .and. &
            cwp .gt. cldmin) then 

          do ibnd = 1, cloud_spec%nbnd
            tauice = ciwp(icol,ilyr) * extice(ilyr,ibnd)
            tauliq = clwp(icol,ilyr) * extliq(ilyr,ibnd)

            if (tauice == 0.0_wp .and. tauliq == 0.0_wp) then
              error_msg = &
                'cloud optics: cloud optical depth is zero'
              return
            endif

            if cloud_Spec%is_lw .eq. 1 then
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
  end function cloud_optics

  !-------------------------------------------------------------------
  !
  ! Ancillary functions
  !
  !-------------------------------------------------------------------

  ! the load function Robert requested
  !  We will need a load function to add the cloud optics data 
  !  (coefficients or lookup tables); only one function is needed and 
  !  all the types can point to it. 
  function load_cloud_optics()
  end function load_cloud_optics

  ! the compute function Robert requested:
  !  We will need a compute function that fills in the optical 
  !  properties arrays (tau/omega0/g or whatever)
  !
  !  compute() should include a variable of ty_spectral_disc as an 
  !  input argument. If the cloud optics variable doesn't know how 
  !  to provide data on that spectral grid an error string should 
  !  be returned. It would also be possible to have a single type 
  !  that can provide both LW and SW values. 
  function calc_optical_props()
  end function calc_optical_props

  function get_irad(cloud_spec,rad,phase,param)
    class(ty_cloud_optics_pade), intent(inout) :: cloud_spec

    ! particle radius
    real(wp), intent(in) :: rad

    ! phase (liq/ice), ext/ssa/asy
    character(len=3), intent(in) :: phase, param

    integer :: get_irad

    ! Local variables
    integer :: irad
    real(wp), dimension(4) :: sizreg

    if (phase .eq. 'liq') then
      if (param .eq. 'ext' .or. param .eq. 'ssa') then
        sizreg = cloud_spec%pade_sizreg_liq_g(:)
      else (param .eq. 'asy') 
        sizreg = cloud_spec%pade_sizreg_liq_tau(:)
      endif

      do irad = 1, 2
        if (rad .gt. sizreg(irad) .and. rad .le. sizreg(irad+1)) &
          get_irad = irad
      enddo
    endif ! liquid

    if (phase .eq. 'ice') then
      sizreg = cloud_spec%pade_sizreg_ice_g(:)
      if (cloud_spec%icergh .eq. 1 .and. cloud_spec%is_lw .and. &
         (param .eq. 'ssa' .or. param .eq. 'asy')) &
        sizreg = cloud_spec%pade_sizreg_ice_tau(:)

      do irad = 1, 3 
        if (rad .gt. sizreg(irad) .and. rad .le. sizreg(irad+1)) &
          get_irad = irad
      enddo
    endif ! ice
  end function get_irad

  ! RP comments: 
  ! Much better to replace repeated code with a single routine to 
  ! evaluate Pade approximants since the formula is entirely general 
  ! https://en.wikipedia.org/wiki/Padé_approximant
  !
  !   function pade(n_num, num, n_den, denom)
  !     integer, intent(in) :: n_num, n_den
  !     real(wp), intent(in) :: num(n_num), den(n_den) 
  !
  !  Direct exponentiation is very expensive. Please evaluate 
  !  polynomals with Horner 
  !(p(1) + p(2)*reff + p(3)*reff**2  = p(1) + reff*(p(2) + reff * p3))
  !  or better as a loop over the number of terms 

  ! If ssa is fit as coalbedo (1 - ssa) that can be done where this 
  ! routine is called. 

!---------------------------------------------------------------------------
  function pade_ext(cloud_spec,p,reff)
     class(ty_cloud_optics_pade), intent(inout) :: cloud_spec
     real(wp), intent(in) :: p(6)            ! extinction Pade coefficients
     real(wp), intent(in) :: reff            ! particle radius (microns)
     real(wp) :: pade_ext
! Pade formulation: Extinction Coefficient (Hogan and Bozzo, ECMWF, TM787, 2016)
     pade_ext = (p(1) + p(2)*reff + p(3)*reff**2) / &
                (1.0_wp + p(4)*reff + p(5)*reff**2 + p(6)*reff**3)
  end function pade_ext

!---------------------------------------------------------------------------
  function pade_ssa(cloud_spec,p,reff)
     class(ty_cloud_optics_pade), intent(inout) :: cloud_spec
     real(wp), intent(in) :: p(5)            ! ssa or g Pade coefficients
     real(wp), intent(in) :: reff            ! particle radius (microns)
     real(wp) :: pade_ssa
! Pade formulation: Single Scattering Albedo (Hogan and Bozzo, ECMWF, TM787, 2016)
     pade_ssa = 1.0_wp - (p(1)+p(2)*reff+p(3)*reff**2) / &
                      (1.0_wp+p(4)*reff+p(5)*reff**2)
  end function pade_ssa

!---------------------------------------------------------------------------
  function pade_asy(cloud_spec,p,reff)
     class(ty_cloud_optics_pade), intent(inout) :: cloud_spec
     real(wp), intent(in) :: p(5)            ! ssa or g Pade coefficients
     real(wp), intent(in) :: reff            ! particle radius (microns)
     real(wp) :: pade_asy
! Pade formulation: Asymmetry Parameter (Hogan and Bozzo, ECMWF, TM787, 2016)
     pade_asy = ((p(1)+p(2)*reff+p(3)*reff**2) / &
                (1.0_wp+p(4)*reff+p(5)*reff**2))
  end function pade_asy

  function is_lw(fileName) 
    !----------------------------------------------------------------
    !
    ! Does this file contain variables needed to do LW calculations? 
    ! I don't think this is necessary anymore now that is_lw is a 
    ! an attribute in this Pade cloud optics class
    !----------------------------------------------------------------

    character(len=*), intent(in   ) :: fileName
    logical                         :: is_lw
    
    integer :: ncid, dimid, status

    is_lw = .false.
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) & 
      call stop_on_err("is_lw: can't find file " // trim(fileName))

    if(nf90_inq_dimid(ncid, 'nband_lw', dimid) == NF90_NOERR) &
      is_lw = .true.

    ncid = nf90_close(ncid) 
  end function is_lw

  subroutine stop_on_err(msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg

    if(msg /= "") then
      write(error_unit, *) msg
      stop
    end if
  end subroutine

  function get_dim_length(ncid, dimname)
    !
    ! Get the length of a dimension from an open netCDF file
    !  This is unfortunately a two-step process
    !

    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer :: get_dim_length, dimid

    ! so....get_dim_length can only be undefined or 0?
    if(nf90_inq_dimid(ncid, trim(dimname), dimid) == NF90_NOERR) then
      if(nf90_inquire_dimension(ncid, dimid, len=get_dim_length) /= &
         NF90_NOERR) get_dim_length = 0
    else
      get_dim_length = 0
    end if
  end function get_dim_length

  function get_data_size(ncid, varName, n)
    !
    ! Returns the extents of a netcdf variable on disk
    !

    integer,          intent(in) :: ncid, n
    character(len=*), intent(in) :: varName
    integer                      :: get_data_size(n), i, &
                                    varid, ndims, dimids(n)

    get_data_size(n) = -1

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("get_data_size: can't find variable " // &
      trim(varName))

    if(nf90_inquire_variable(ncid, varid, ndims = ndims) /= &
       NF90_NOERR) &
      call stop_on_err(&
      "get_data_size: can't get information for variable " // &
      trim(varName))

    if(ndims /= n) &
      call stop_on_err("get_data_size:  variable " // &
      trim(varName) // &
      " has the wrong number of dimensions" )
    
    if(nf90_inquire_variable(ncid, varid, dimids = dimids) /= &
       NF90_NOERR) &
      call stop_on_err(&
      "get_data_size: can't read dimension ids for variable " // &
      trim(varName))

    do i = 1, n
      if(nf90_inquire_dimension(&
         ncid, dimids(i), len = get_data_size(i)) /= NF90_NOERR) &
        call stop_on_err(&
        "get_data_size: can't get dim lengths for variable " // &
        trim(varName))
    end do
  end function get_data_size

  function read_scalar(ncid, varName)
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

  function read_1d_field(ncid, varName, nx)
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

  function read_2d_field(ncid, varName, nx, ny)
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

  function read_3d_field(ncid, varName, nx, ny, nz)
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

  function read_4d_field(ncid, varName, nw, nx, ny, nz)
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

  function var_exists(ncid, varName) 
    !
    ! Does this variable exist (have a valid var_id) in the open 
    ! netCDF file? 
    ! 

    integer,          intent(in) :: ncid 
    character(len=*), intent(in) :: varName
    logical                      :: var_exists
    integer                      :: varId 

    var_exists = &
      nf90_inq_varid(ncid, trim(varName), varid) == NF90_NOERR
  end function var_exists 

  function dim_exists(ncid, dimName) 
    !
    ! Does this dimension exist (have a valid dim_id) in the open 
    ! netCDF file? 
    ! 

    integer,          intent(in) :: ncid 
    character(len=*), intent(in) :: dimName
    logical                      :: dim_exists
    integer                      :: dimid

    dim_exists = &
      nf90_inq_dimid(ncid, trim(dimName), dimid) == NF90_NOERR
  end function dim_exists 

  subroutine create_dim(ncid, dimName, dimLength) 
    !
    ! Check to see if a dimiable with this name exists in the file
    !   If so, check against current size
    !   If not, create with specified dimensions  
    ! 
    integer,          intent(in) :: ncid 
    character(len=*), intent(in) :: dimName
    integer,          intent(in) :: dimLength
    integer                      :: i, dimid 

    if(dim_exists(ncid, dimName)) then 
      if (dimLength /= get_dim_length(ncid, trim(dimName))) & 
          call stop_on_err("dim " // trim(dimName) // &
          " is present but incorrectly sized.") 
    else 
      if(nf90_redef(ncid) /= NF90_NOERR) & 
        call stop_on_err(&
        "create_dim: can't put file into redefine mode") 

      if(&
        nf90_def_dim(ncid, dimName, dimLength, dimid) /= NF90_NOERR) & 
        call stop_on_err("create_dim: can't define dimension " // &
        trim(dimName))

      if(nf90_enddef(ncid) /= NF90_NOERR) & 
        call stop_on_err("create_dim: can't end redefinition??")
    end if 
  end subroutine create_dim 

  subroutine create_var(ncid, varName, dimNames, dimLengths, dataType)
    !
    ! Check to see if a variable with this name exists in the file
    !   If so, check against current size
    !   If not, create with specified dimensions  
    ! datatype: NF90_DOUBLE, NF90_FLOAT, NF90_INT, etc.
    ! 

    integer,          intent(in) :: ncid 
    character(len=*), intent(in) :: varName
    character(len=*), intent(in) :: dimNames(:) 
    integer,          intent(in) :: dimLengths(:)
    integer, optional, intent(in) :: dataType
    
    integer :: i, varid, xtype, dimIds(size(dimNames))
    
    if(var_exists(ncid, varName)) then 
      do i = 1, size(dimNames) 
        if (dimLengths(i) /= &
            get_dim_length(ncid, trim(dimNames(i)))) & 
          call stop_on_err("Variable " // trim(varName) // &
          " is present but incorrectly sized.") 
      end do 
    else 
      do i = 1, size(dimNames) 
        if(nf90_inq_dimid(&
           ncid, trim(dimNames(i)), dimIds(i)) /= NF90_NOERR) & 
          call stop_on_err("create_var: Can't get id for dim " // &
          trim(dimnames(i)))
      end do

      if(nf90_redef(ncid) /= NF90_NOERR) & 
        call stop_on_err(&
          "create_var: can't put file into redefine mode") 

      xtype = NF90_DOUBLE

      if(present(dataType)) xtype = dataType

      if(nf90_def_var(ncid, varName, xtype, dimIds, varid) /= &
        NF90_NOERR) call stop_on_err(&
        "create_var: can't define variable " // trim(varName))

      if(nf90_enddef(ncid) /= NF90_NOERR) & 
        call stop_on_err("create_dim: can't end redefinition??")
    end if 
  end subroutine create_var 

  !------------------------------------------------------------------
  ! Writing functions 
  !------------------------------------------------------------------
  function write_1d_int_field(ncid, varName, var) result(err_msg) 
    integer,                intent(in) :: ncid 
    character(len=*),       intent(in) :: varName
    integer, dimension(:),  intent(in) :: var
    character(len=128)                 :: err_msg 
     
    integer :: varid 
    
    err_msg = ""
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) then
      err_msg = "write_field: can't find variable " // trim(varName)
      return
    end if 
    if(nf90_put_var(ncid, varid, var)  /= NF90_NOERR) &
      err_msg = "write_field: can't write variable " // trim(varName)
  end function write_1d_int_field

  function write_2d_int_field(ncid, varName, var) result(err_msg) 
    integer,                  intent(in) :: ncid 
    character(len=*),         intent(in) :: varName
    integer, dimension(:,:),  intent(in) :: var
    character(len=128)                   :: err_msg 
     
    integer :: varid 
    
    err_msg = ""
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) then
      err_msg = "write_field: can't find variable " // trim(varName)
      return
    end if 
    if(nf90_put_var(ncid, varid, var)  /= NF90_NOERR) &
      err_msg = "write_field: can't write variable " // trim(varName)
  end function write_2d_int_field

  function write_1d_field(ncid, varName, var) result(err_msg) 
    integer,                intent(in) :: ncid 
    character(len=*),       intent(in) :: varName
    real(wp), dimension(:), intent(in) :: var
    character(len=128)                 :: err_msg 
     
    integer :: varid 
    
    err_msg = ""
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) then
      err_msg = "write_field: can't find variable " // trim(varName)
      return
    end if 
    if(nf90_put_var(ncid, varid, var)  /= NF90_NOERR) &
      err_msg = "write_field: can't write variable " // trim(varName)
  end function write_1d_field

  function write_2d_field(ncid, varName, var) result(err_msg) 
    integer,                  intent(in) :: ncid 
    character(len=*),         intent(in) :: varName
    real(wp), dimension(:,:), intent(in) :: var
    character(len=128)                   :: err_msg 
     
    integer :: varid 
    
    err_msg = ""
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) then
      err_msg = "write_field: can't find variable " // trim(varName)
      return
    end if 
    if(nf90_put_var(ncid, varid, var)  /= NF90_NOERR) &
      err_msg = "write_field: can't write variable " // trim(varName)
  end function write_2d_field

  function write_3d_field(ncid, varName, var) result(err_msg) 
    integer,                    intent(in) :: ncid 
    character(len=*),           intent(in) :: varName
    real(wp), dimension(:,:,:), intent(in) :: var
    character(len=128)                     :: err_msg 
     
    integer :: varid 
    
    err_msg = ""
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) then
      err_msg = "write_field: can't find variable " // trim(varName)
      return
    end if 
    if(nf90_put_var(ncid, varid, var)  /= NF90_NOERR) &
      err_msg = "write_field: can't write variable " // trim(varName)  
  end function write_3d_field

  function write_4d_field(ncid, varName, var) result(err_msg) 
    integer,                    intent(in) :: ncid 
    character(len=*),           intent(in) :: varName
    real(wp), dimension(:,:,:,:), intent(in) :: var
    character(len=128)                     :: err_msg 
     
    integer :: varid 
    
    err_msg = ""
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) then
      err_msg = "write_field: can't find variable " // trim(varName)
      return
    end if 
    if(nf90_put_var(ncid, varid, var)  /= NF90_NOERR) &
      err_msg = "write_field: can't write variable " // trim(varName)  
  end function write_4d_field

end module mo_cloud_optics_pade
