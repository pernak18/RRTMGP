! Module: mo_cloud_optics

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
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description:
! This is the interface for routines that receive cloud physical properties
! and return cloud optical properties by band using LUT input data. 
!

module mo_cloud_optics_lut
  use mo_rte_kind,            only: wp
  use mo_spectral_disc,       only: ty_spectral_disc
  use mo_cloud_optical_props, only: ty_cloud_optical_props_arry, &
                                    ty_cloud_optical_props_1scl, &
                                    ty_cloud_optical_props_2str, &
                                    ty_cloud_optical_props_nstr
  use netcdf

  implicit none
  private
  public :: is_lw
  public :: stop_on_err, get_dim_length, read_field, create_var, dim_exists, var_exists, write_field

  interface read_field
    module procedure read_scalar, read_1d_field, read_2d_field, read_3d_field, read_4d_field
  end interface
  interface write_field 
    module procedure write_1d_int_field, write_2d_int_field, &
                     write_1d_field, write_2d_field, write_3d_field, write_4d_field
  end interface 

  ! -----------------------------------------------------------------------------------
  type, extends(ty_spectral_disc), public :: ty_cloud_optics_lut
!    private

    ! User input
    ! Ice surface roughness category - needed for Yang (2013) ice optics parameterization
    integer :: icergh                                   ! (1 = none, 2 = medium, 3 = high)

    ! Delta-scaling for cloud optical properties (shortwave only)
    integer :: idelscl                                  ! (0 = no delta-scaling, 1 = with delta-scaling)

    ! Cloud physical property dimensions
    integer :: ncol, nlay
    ! Cloud physical properties                         ! (ncol,nlay)
    real(wp), dimension(:,:), allocatable :: cldfrac    ! cloud fraction
    real(wp), dimension(:,:), allocatable :: ciwp       ! cloud ice water path
    real(wp), dimension(:,:), allocatable :: clwp       ! cloud liquid water path
    real(wp), dimension(:,:), allocatable :: rei        ! cloud ice particle effective size (microns)
    real(wp), dimension(:,:), allocatable :: rel        ! cloud liquid particle effective radius (microns)

    ! Model input: LUT
    ! Lookup table cloud coefficient dimensions
    integer :: nband_lw,        &
               nband_sw,        &
               nrghice
    integer :: nsize_liq,       & 
               nsize_ice
    ! Lookup table cloud coefficients
    real(wp), dimension(:,:    ), allocatable :: lut_extliq     ! (nsize_liq, nband_lw)
    real(wp), dimension(:,:    ), allocatable :: lut_ssaliq     ! (nsize_liq, nband_lw)
    real(wp), dimension(:,:    ), allocatable :: lut_asyliq     ! (nsize_liq, nband_lw)
    real(wp), dimension(:,:,:  ), allocatable :: lut_extice     ! (nsize_ice, nband_lw, nrghice)
    real(wp), dimension(:,:,:  ), allocatable :: lut_ssaice     ! (nsize_ice, nband_lw, nrghice)
    real(wp), dimension(:,:,:  ), allocatable :: lut_asyice     ! (nsize_ice, nband_lw, nrghice)

! ------------------------------------------------------------------------------------------
  contains

    generic, public :: init_cldopt  => init_cldopt_lw, init_cldopt_sw
    generic, public :: cloud_optics => cloud_optics_lw, cloud_optics_sw
    ! Internal procedures
    procedure, private :: init_cldopt_lw
    procedure, private :: init_cldopt_sw
    procedure, private :: cloud_optics_lw
    procedure, private :: cloud_optics_sw

  end type ty_cloud_optics_lut

  contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Cloud optics initialization functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  ! Initialize object based on data read from netCDF file however the user desires. 
  !   Lookup table data used for longwave cloud property conversion
  ! 
  function init_cldopt_lw(cloud_spec, cld_coeff_file, is_lw) result(error_msg)

    ! LUT
    class(ty_cloud_optics_lut), intent(inout) :: cloud_spec

    ! Cloud coefficient optical property input file - LUT
    character(len=*), intent(in) :: cld_coeff_file

    logical, intent(in)  :: is_lw

    character(len = 128) error_msg

    ! -----------------
    ! Local variables
    integer :: ncid
    integer :: nband_lw
    integer :: nrghice
    integer :: nsize_liq
    integer :: nsize_ice

    ! -----------------
    error_msg = ""

    ! Open cloud optical property coefficient file
    if(nf90_open(trim(cld_coeff_file), NF90_WRITE, ncid) /= NF90_NOERR) & 
       call stop_on_err("init_cldopt_lw(): can't open file " // trim(cld_coeff_file))

    cloud_spec%nband_lw  = get_dim_length(ncid,'nband_lw') 
    cloud_spec%nrghice   = get_dim_length(ncid,'nrghice') 
    cloud_spec%nsize_liq = get_dim_length(ncid,'nsize_liq') 
    cloud_spec%nsize_ice = get_dim_length(ncid,'nsize_ice') 

    nband_lw  = cloud_spec%nband_lw
    nrghice   = cloud_spec%nrghice
    nsize_liq = cloud_spec%nsize_liq
    nsize_ice = cloud_spec%nsize_ice

    ! Allocate cloud property lookup table input arrays
    allocate(cloud_spec%lut_extliq(nsize_liq, nband_lw))
    allocate(cloud_spec%lut_ssaliq(nsize_liq, nband_lw))
    allocate(cloud_spec%lut_asyliq(nsize_liq, nband_lw))
    allocate(cloud_spec%lut_extice(nsize_ice, nband_lw, nrghice))
    allocate(cloud_spec%lut_ssaice(nsize_ice, nband_lw, nrghice))
    allocate(cloud_spec%lut_asyice(nsize_ice, nband_lw, nrghice))

    cloud_spec%lut_extliq   = read_field(ncid, 'lut_extliq_lw',  nsize_liq, nband_lw) 
    cloud_spec%lut_ssaliq   = read_field(ncid, 'lut_ssaliq_lw',  nsize_liq, nband_lw) 
    cloud_spec%lut_asyliq   = read_field(ncid, 'lut_asyliq_lw',  nsize_liq, nband_lw) 
    cloud_spec%lut_extice   = read_field(ncid, 'lut_extice_lw',  nsize_ice, nband_lw, nrghice) 
    cloud_spec%lut_ssaice   = read_field(ncid, 'lut_ssaice_lw',  nsize_ice, nband_lw, nrghice) 
    cloud_spec%lut_asyice   = read_field(ncid, 'lut_asyice_lw',  nsize_ice, nband_lw, nrghice) 

    ncid = nf90_close(ncid) 

  end function init_cldopt_lw

  !--------------------------------------------------------------------------------------------------------------------
  ! Initialize object based on data read from netCDF file however the user desires. 
  !   Lookup table data used for shortwave cloud property conversion
  ! 
  function init_cldopt_sw(cloud_spec, cld_coeff_file) result(error_msg)

    ! LUT
    class(ty_cloud_optics_lut), intent(inout) :: cloud_spec

    ! Cloud coefficient optical property input file - LUT
    character(len=*), intent(in) :: cld_coeff_file

    character(len = 128) error_msg

    ! -----------------
    ! Local variables
    integer :: ncid
    integer :: nband_sw
    integer :: nrghice
    integer :: nsize_liq
    integer :: nsize_ice

    ! -----------------
    error_msg = ""

    ! Open cloud optical property coefficient file
    if(nf90_open(trim(cld_coeff_file), NF90_WRITE, ncid) /= NF90_NOERR) & 
       call stop_on_err("init_cldopt_sw(): can't open file " // trim(cld_coeff_file))

    cloud_spec%nband_sw  = get_dim_length(ncid,'nband_sw') 
    cloud_spec%nrghice   = get_dim_length(ncid,'nrghice') 
    cloud_spec%nsize_liq = get_dim_length(ncid,'nsize_liq') 
    cloud_spec%nsize_ice = get_dim_length(ncid,'nsize_ice') 

    nband_sw  = cloud_spec%nband_sw
    nrghice   = cloud_spec%nrghice
    nsize_liq = cloud_spec%nsize_liq
    nsize_ice = cloud_spec%nsize_ice

    ! Allocate cloud property lookup table input arrays
    allocate(cloud_spec%lut_extliq(nsize_liq, nband_sw))
    allocate(cloud_spec%lut_ssaliq(nsize_liq, nband_sw))
    allocate(cloud_spec%lut_asyliq(nsize_liq, nband_sw))
    allocate(cloud_spec%lut_extice(nsize_ice, nband_sw, nrghice))
    allocate(cloud_spec%lut_ssaice(nsize_ice, nband_sw, nrghice))
    allocate(cloud_spec%lut_asyice(nsize_ice, nband_sw, nrghice))

    cloud_spec%lut_extliq   = read_field(ncid, 'lut_extliq_sw',  nsize_liq, nband_sw) 
    cloud_spec%lut_ssaliq   = read_field(ncid, 'lut_ssaliq_sw',  nsize_liq, nband_sw) 
    cloud_spec%lut_asyliq   = read_field(ncid, 'lut_asyliq_sw',  nsize_liq, nband_sw) 
    cloud_spec%lut_extice   = read_field(ncid, 'lut_extice_sw',  nsize_ice, nband_sw, nrghice) 
    cloud_spec%lut_ssaice   = read_field(ncid, 'lut_ssaice_sw',  nsize_ice, nband_sw, nrghice) 
    cloud_spec%lut_asyice   = read_field(ncid, 'lut_asyice_sw',  nsize_ice, nband_sw, nrghice) 

    ncid = nf90_close(ncid) 

  end function init_cldopt_sw

  ! ------------------------------------------------------------------------------
  !
  ! Derive LW cloud optical properties from provided cloud physical properties
  !
  ! ------------------------------------------------------------------------------
  function cloud_optics_lw( &
  ! Input
                   cloud_spec, &
                   ncol, nlayers, nbndlw, is_lw, &
                   cldfrac, clwp, ciwp, rel, rei, &
  ! Output
                   cloud_optical_props) &
                   result(error_msg)
  ! ------------------------------------------------------------------------------

  ! Purpose:  Compute the cloud optical properties for each cloudy layer.

  ! ------- Input -------

      class(ty_cloud_optics_lut), intent(in) :: cloud_spec
                                                 ! cloud specification data
      integer, intent(in) :: ncol                ! total number of columns
      integer, intent(in) :: nlayers             ! total number of layers
      integer, intent(in) :: nbndlw              ! number of LW bands

      logical, intent(in) :: is_lw

      real(wp), intent(in) :: cldfrac(:,:)       ! cloud fraction
                                                 !    Dimensions: (ncol,nlayers)
      real(wp), intent(in) :: ciwp(:,:)          ! cloud ice water path
                                                 !    Dimensions: (ncol,nlayers)
      real(wp), intent(in) :: clwp(:,:)          ! cloud liquid water path
                                                 !    Dimensions: (ncol,nlayers)
      real(wp), intent(in) :: rei(:,:)           ! cloud ice particle effective size (microns)
                                                 !    Dimensions: (ncol,nlayers)
      real(wp), intent(in) :: rel(:,:)           ! cloud liquid particle effective radius (microns)
                                                 !    Dimensions: (ncol,nlayers)

! ------- Output -------
      class(ty_cloud_optical_props_arry), intent(inout) :: cloud_optical_props
                                                 ! Dimensions: (ncol,nlayers,nbndlw)

! ------- Local -------

      character(len=128)    :: error_msg

      real(wp) :: cwp                            ! cloud water path (liquid + ice)
      real(wp), parameter :: cldmin = 1.e-20     ! minimum value for cloud quantities
      real(wp) :: radliq                         ! cloud liquid droplet radius (microns)
      real(wp) :: radice                         ! cloud ice effective size (microns)
      real(wp) :: factor                         ! 
      real(wp) :: fint                           ! 

      integer :: index                           !
      integer :: icol, ilyr, ibnd                !
      integer :: irad, iradg                     !
      integer :: icergh                          ! ice surface roughness
                                                 ! (1 = none, 2 = medium, 3 = high)

      real(wp) :: extliq_lw(nlayers,nbndlw)      ! LW liquid extinction coefficient
      real(wp) :: ssaliq_lw(nlayers,nbndlw)      ! LW liquid single scattering albedo
      real(wp) :: asyliq_lw(nlayers,nbndlw)      ! LW liquid asymmetry parameter

      real(wp) :: extice_lw(nlayers,nbndlw)      ! LW ice extinction coefficients
      real(wp) :: ssaice_lw(nlayers,nbndlw)      ! LW ice single scattering albedo
      real(wp) :: asyice_lw(nlayers,nbndlw)      ! LW ice asymmetry parameter

      real(wp) :: tauliq_lw                      ! LW liquid cloud extinction optical depth - no delta scaling
      real(wp) :: tauice_lw                      ! LW ice cloud extinction optical depth - no delta scaling

      real(wp) :: scatice                        ! Ice scattering term
      real(wp) :: scatliq                        ! Liquid scattering term

      real(wp) :: asycld                         ! LW asymmetry parameter - local

! ------- Definitions -------

! ------- Error checking -------
      error_msg = ''
      icergh = cloud_spec%icergh
      if (icergh < 1 .or. icergh > 3) then
         error_msg = 'cloud optics: cloud ice surface roughness flag is out of bounds'
         return
      endif

! Initialize
      extliq_lw(:,:) = 0.0_wp
      ssaliq_lw(:,:) = 1.0_wp
      asyliq_lw(:,:) = 0.0_wp
      extice_lw(:,:) = 0.0_wp
      ssaice_lw(:,:) = 1.0_wp
      asyice_lw(:,:) = 0.0_wp

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
   do icol = 1, ncol

! Main layer loop
      do ilyr = 1, nlayers
         cwp = ciwp(icol,ilyr) + clwp(icol,ilyr)
         if (cldfrac(icol,ilyr) .ge. cldmin .and. cwp .ge. cldmin) then

! Liquid OP
            radliq = rel(icol,ilyr)
            if (radliq .gt. 0.0_wp .and. clwp(icol,ilyr) .gt. 0.0_wp) then 
! For liquid OP, particle size is limited to 2.5 to 21.0 microns
               if (radliq .lt. 2.5_wp .or. radliq .gt. 21.0_wp) then 
                  error_msg = 'cloud optics: liquid effective radius is out of bounds'
                  return
               endif
               factor = radliq - 1.5_wp
               index = int(factor)
               if (index .eq. 0) index = 1
               fint = factor - real(index)
! Liquid
               do ibnd = 1, nbndlw
                  extliq_lw(ilyr,ibnd) = cloud_spec%lut_extliq(index,ibnd) + &
                                 fint * (cloud_spec%lut_extliq(index+1,ibnd) - &
                                         cloud_spec%lut_extliq(index,ibnd))
                  ssaliq_lw(ilyr,ibnd) = cloud_spec%lut_ssaliq(index,ibnd) + &
                                 fint * (cloud_spec%lut_ssaliq(index+1,ibnd) - &
                                         cloud_spec%lut_ssaliq(index,ibnd))
                  asyliq_lw(ilyr,ibnd) = cloud_spec%lut_asyliq(index,ibnd) + &
                                 fint * (cloud_spec%lut_asyliq(index+1,ibnd) - &
                                         cloud_spec%lut_asyliq(index,ibnd))
               enddo

            endif

! Ice OP for requested ice roughness (icergh)
            radice = rei(icol,ilyr)
            if (radice .gt. 0.0_wp .and. ciwp(icol,ilyr) .gt. 0.0_wp) then 
! For Yang (2013) ice OP, particle size is limited to 10.0 to 180.0 microns
               if (radice .lt. 10.0_wp .or. radice .gt. 180.0_wp) then
                  error_msg = 'cloud optics: ice effective radius is out of bounds'
                  return
               endif
               factor = radice / 10._wp
               index = int(factor)
               fint = factor - real(index)
! Ice
               do ibnd = 1, nbndlw
                  extice_lw(ilyr,ibnd) = cloud_spec%lut_extice(index,ibnd,icergh) + &
                                 fint * (cloud_spec%lut_extice(index+1,ibnd,icergh) - &
                                         cloud_spec%lut_extice(index,ibnd,icergh))
                  ssaice_lw(ilyr,ibnd) = cloud_spec%lut_ssaice(index,ibnd,icergh) + &
                                 fint * (cloud_spec%lut_ssaice(index+1,ibnd,icergh) - &
                                         cloud_spec%lut_ssaice(index,ibnd,icergh))
                  asyice_lw(ilyr,ibnd) = cloud_spec%lut_asyice(index,ibnd,icergh) + &
                                 fint * (cloud_spec%lut_asyice(index+1,ibnd,icergh) - &
                                         cloud_spec%lut_asyice(index,ibnd,icergh))
               enddo

            endif
         endif
      enddo

! End column loop
   enddo

! Combine liquid and ice contributions into total cloud optical properties
   do icol = 1, ncol
      do ilyr = 1, nlayers
         cwp = ciwp(icol,ilyr) + clwp(icol,ilyr)
         if (cldfrac(icol,ilyr) .gt. cldmin .and. cwp .gt. cldmin) then 

         do ibnd = 1, nbndlw
            tauice_lw = ciwp(icol,ilyr) * extice_lw(ilyr,ibnd)
            tauliq_lw = clwp(icol,ilyr) * extliq_lw(ilyr,ibnd)

            if (tauice_lw == 0.0_wp .and. tauliq_lw == 0.0_wp) then
               error_msg = 'cloud optics: longwave cloud optical depth is zero'
               return
            endif

            scatice = ssaice_lw(ilyr,ibnd) * tauice_lw
            scatliq = ssaliq_lw(ilyr,ibnd) * tauliq_lw

! No delta-scaling
              select type(cloud_optical_props)
                type is (ty_cloud_optical_props_1scl)
                  cloud_optical_props%taucld(icol,ilyr,ibnd) = tauice_lw + tauliq_lw
                type is (ty_cloud_optical_props_2str)
                  cloud_optical_props%taucld(icol,ilyr,ibnd) = tauice_lw + tauliq_lw
                  cloud_optical_props%ssacld(icol,ilyr,ibnd) = (scatice + scatliq) / &
                                        cloud_optical_props%taucld(icol,ilyr,ibnd)
                  cloud_optical_props%asycld(icol,ilyr,ibnd) = &
                       (scatice * asyice_lw(ilyr,ibnd) + scatliq * asyliq_lw(ilyr,ibnd)) / &
                       (scatice + scatliq)
                type is (ty_cloud_optical_props_nstr)
                  if (clwp(icol,ilyr) .gt. 0.0_wp) then 
                     error_msg = 'cloud optics: n-stream option not available in longwave for liquid clouds'
                     return
                  else
                    cloud_optical_props%taucld(icol,ilyr,ibnd) = tauice_lw + tauliq_lw
                    cloud_optical_props%ssacld(icol,ilyr,ibnd) = (scatice + scatliq) / &
                                          cloud_optical_props%taucld(icol,ilyr,ibnd)
                    asycld = &
                         (scatice * asyice_lw(ilyr,ibnd) + scatliq * asyliq_lw(ilyr,ibnd)) / &
                         (scatice + scatliq)
                    cloud_optical_props%pcld(1,icol,ilyr,ibnd) = 1.0_wp
                    cloud_optical_props%pcld(2,icol,ilyr,ibnd) = asycld
                    cloud_optical_props%pcld(3,icol,ilyr,ibnd) = asycld**2
                  endif
              end select

         enddo

         endif
      enddo
   enddo

  end function cloud_optics_lw

  ! ------------------------------------------------------------------------------
  !
  ! Derive SW cloud optical properties from provided cloud physical properties
  !
  !--------------------------------------------------------------------------------------------------------------------
  function cloud_optics_sw( &
  ! Input
                   cloud_spec, &
                   ncol, nlayers, nbndsw, &
                   cldfrac, clwp, ciwp, rel, rei, &
  ! Output
                   cloud_optical_props) &
                   result(error_msg)
  ! ------------------------------------------------------------------------------

  ! Purpose:  Compute the cloud optical properties for each cloudy layer.

  ! ------- Input -------

      class(ty_cloud_optics_lut), intent(in) :: cloud_spec
                                                 ! cloud specification data
      integer, intent(in) :: ncol                ! total number of columns
      integer, intent(in) :: nlayers             ! total number of layers
      integer, intent(in) :: nbndsw              ! number of SW bands

      real(wp), intent(in) :: cldfrac(:,:)       ! cloud fraction
                                                 !    Dimensions: (ncol,nlayers)
      real(wp), intent(in) :: ciwp(:,:)          ! cloud ice water path
                                                 !    Dimensions: (ncol,nlayers)
      real(wp), intent(in) :: clwp(:,:)          ! cloud liquid water path
                                                 !    Dimensions: (ncol,nlayers)
      real(wp), intent(in) :: rei(:,:)           ! cloud ice particle effective size (microns)
                                                 !    Dimensions: (ncol,nlayers)
      real(wp), intent(in) :: rel(:,:)           ! cloud liquid particle effective radius (microns)
                                                 !    Dimensions: (ncol,nlayers)

! ------- Output -------
      class(ty_cloud_optical_props_arry), intent(inout) :: cloud_optical_props
                                                 ! Dimensions: (ncol,nlayers,nbndlw)

! ------- Local -------

      character(len=128)    :: error_msg

      real(wp) :: cwp                            ! cloud water path (liquid + ice)
      real(wp), parameter :: cldmin = 1.e-20     ! minimum value for cloud quantities
      real(wp) :: radliq                         ! cloud liquid droplet radius (microns)
      real(wp) :: radice                         ! cloud ice effective size (microns)
      real(wp) :: factor                         ! 
      real(wp) :: fint                           ! 

      integer :: index                           !
      integer :: icol, ilyr, ibnd                !
      integer :: irad, iradg                     !
      integer :: icergh                          ! ice surface roughness
                                                 ! (1 = none, 2 = medium, 3 = high)

      real(wp) :: extliq_sw(nlayers,nbndsw)      ! SW liquid extinction coefficients
      real(wp) :: ssaliq_sw(nlayers,nbndsw)      ! SW liquid single scattering albedo
      real(wp) :: asyliq_sw(nlayers,nbndsw)      ! SW liquid asymmetry parameter

      real(wp) :: extice_sw(nlayers,nbndsw)      ! SW ice extinction coefficients
      real(wp) :: ssaice_sw(nlayers,nbndsw)      ! SW ice single scattering albedo
      real(wp) :: asyice_sw(nlayers,nbndsw)      ! SW ice asymmetry parameter

      real(wp) :: tauliq_sw                      ! SW liquid cloud extinction optical depth - no delta scaling
      real(wp) :: tauice_sw                      ! SW ice cloud extinction optical depth - no delta scaling

      real(wp) :: scatice                        ! Ice scattering term
      real(wp) :: scatliq                        ! Liquid scattering term

      real(wp) :: asycld                         ! LW asymmetry parameter - local

      real(wp) :: forwliq, forwice               ! Forward scattering terms
      real(wp) :: tauliq_del, tauice_del         ! Extinction optical depth - with delta scaling
      real(wp) :: ssaliq_del, ssaice_del         ! Single scattering albedo - with delta scaling

! ------- Definitions -------

! ------- Error checking -------
      error_msg = ''
      icergh = cloud_spec%icergh
      if (icergh < 1 .or. icergh > 3) then
         error_msg = 'cloud optics: cloud ice surface roughness flag is out of bounds'
         return
      endif

! Initialize
      extliq_sw(:,:) = 0.0_wp
      ssaliq_sw(:,:) = 1.0_wp
      asyliq_sw(:,:) = 0.0_wp
      extice_sw(:,:) = 0.0_wp
      ssaice_sw(:,:) = 1.0_wp
      asyice_sw(:,:) = 0.0_wp

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
   do icol = 1, ncol

! Main layer loop
      do ilyr = 1, nlayers
         cwp = ciwp(icol,ilyr) + clwp(icol,ilyr)
         if (cldfrac(icol,ilyr) .ge. cldmin .and. cwp .ge. cldmin) then

! Liquid OP
            radliq = rel(icol,ilyr)
            if (radliq .gt. 0.0_wp .and. clwp(icol,ilyr) .gt. 0.0_wp) then 
! For liquid OP, particle size is limited to 2.5 to 21.0 microns
               if (radliq .lt. 2.5_wp .or. radliq .gt. 21.0_wp) then 
                  error_msg = 'cloud optics: liquid effective radius is out of bounds'
                  return
               endif
               factor = radliq - 1.5_wp
               index = int(factor)
               if (index .eq. 0) index = 1
               fint = factor - real(index)
! Liquid
               do ibnd = 1, nbndsw
                  extliq_sw(ilyr,ibnd) = cloud_spec%lut_extliq(index,ibnd) + &
                                 fint * (cloud_spec%lut_extliq(index+1,ibnd) - &
                                         cloud_spec%lut_extliq(index,ibnd))
                  ssaliq_sw(ilyr,ibnd) = cloud_spec%lut_ssaliq(index,ibnd) + &
                                 fint * (cloud_spec%lut_ssaliq(index+1,ibnd) - &
                                         cloud_spec%lut_ssaliq(index,ibnd))
                  asyliq_sw(ilyr,ibnd) = cloud_spec%lut_asyliq(index,ibnd) + &
                                 fint * (cloud_spec%lut_asyliq(index+1,ibnd) - &
                                         cloud_spec%lut_asyliq(index,ibnd))
               enddo

            endif

! Ice OP for requested ice roughness (icergh)
            radice = rei(icol,ilyr)
            if (radice .gt. 0.0_wp .and. ciwp(icol,ilyr) .gt. 0.0_wp) then 
! For Yang (2013) ice OP, particle size is limited to 10.0 to 180.0 microns
               if (radice .lt. 10.0_wp .or. radice .gt. 180.0_wp) then
                  error_msg = 'cloud optics: ice effective radius is out of bounds'
                  return
               endif
               factor = radice / 10._wp
               index = int(factor)
               fint = factor - real(index)
! Ice
               do ibnd = 1, nbndsw
                  extice_sw(ilyr,ibnd) = cloud_spec%lut_extice(index,ibnd,icergh) + &
                                 fint * (cloud_spec%lut_extice(index+1,ibnd,icergh) - &
                                         cloud_spec%lut_extice(index,ibnd,icergh))
                  ssaice_sw(ilyr,ibnd) = cloud_spec%lut_ssaice(index,ibnd,icergh) + &
                                 fint * (cloud_spec%lut_ssaice(index+1,ibnd,icergh) - &
                                         cloud_spec%lut_ssaice(index,ibnd,icergh))
                  asyice_sw(ilyr,ibnd) = cloud_spec%lut_asyice(index,ibnd,icergh) + &
                                 fint * (cloud_spec%lut_asyice(index+1,ibnd,icergh) - &
                                         cloud_spec%lut_asyice(index,ibnd,icergh))
               enddo

            endif
         endif
      enddo

! End column loop
   enddo

! Combine liquid and ice contributions into total cloud optical properties
   do icol = 1, ncol
      do ilyr = 1, nlayers
         cwp = ciwp(icol,ilyr) + clwp(icol,ilyr)
         if (cldfrac(icol,ilyr) .gt. cldmin .and. cwp .gt. cldmin) then 

         do ibnd = 1, nbndsw
            tauice_sw = ciwp(icol,ilyr) * extice_sw(ilyr,ibnd)
            tauliq_sw = clwp(icol,ilyr) * extliq_sw(ilyr,ibnd)

            if (tauice_sw == 0.0_wp .and. tauliq_sw == 0.0_wp) then
               error_msg = 'cloud optics: shortwave cloud optical depth is zero'
               return
            endif

! Define final output with or without delta-scaling
            select case (cloud_spec%idelscl)

! No delta-scaling
            case(0)

              select type(cloud_optical_props)
                type is (ty_cloud_optical_props_1scl)
                  cloud_optical_props%taucld(icol,ilyr,ibnd) = tauice_sw + tauliq_sw
                type is (ty_cloud_optical_props_2str)
                  cloud_optical_props%taucld(icol,ilyr,ibnd) = tauice_sw + tauliq_sw
                  scatice = ssaice_sw(ilyr,ibnd) * tauice_sw
                  scatliq = ssaliq_sw(ilyr,ibnd) * tauliq_sw
                  cloud_optical_props%ssacld(icol,ilyr,ibnd) = (scatice + scatliq) / &
                                        cloud_optical_props%taucld(icol,ilyr,ibnd)
                  cloud_optical_props%asycld(icol,ilyr,ibnd) = &
                       (scatice * asyice_sw(ilyr,ibnd) + scatliq * asyliq_sw(ilyr,ibnd)) / &
                       (scatice + scatliq)
                type is (ty_cloud_optical_props_nstr)
                  cloud_optical_props%taucld(icol,ilyr,ibnd) = tauice_sw + tauliq_sw
                  scatice = ssaice_sw(ilyr,ibnd) * tauice_sw
                  scatliq = ssaliq_sw(ilyr,ibnd) * tauliq_sw
                  cloud_optical_props%ssacld(icol,ilyr,ibnd) = (scatice + scatliq) / &
                                        cloud_optical_props%taucld(icol,ilyr,ibnd)
                  asycld = &
                       (scatice * asyice_sw(ilyr,ibnd) + scatliq * asyliq_sw(ilyr,ibnd)) / &
                       (scatice + scatliq)
                  cloud_optical_props%pcld(1,icol,ilyr,ibnd) = 1.0_wp
                  cloud_optical_props%pcld(2,icol,ilyr,ibnd) = asycld
                  cloud_optical_props%pcld(3,icol,ilyr,ibnd) = asycld**2
              end select

! With delta-scaling
            case(1)

              forwliq = asyliq_sw(ilyr,ibnd)**2
              forwice = asyice_sw(ilyr,ibnd)**2
              tauliq_del = (1.0_wp - forwliq * ssaliq_sw(ilyr,ibnd)) * tauliq_sw
              tauice_del = (1.0_wp - forwice * ssaice_sw(ilyr,ibnd)) * tauice_sw

              select type(cloud_optical_props)
                type is (ty_cloud_optical_props_1scl)
                  cloud_optical_props%taucld(icol,ilyr,ibnd) = tauliq_del + tauice_del
                type is (ty_cloud_optical_props_2str)
                  cloud_optical_props%taucld(icol,ilyr,ibnd) = tauliq_del + tauice_del
                    ssaliq_del = ssaliq_sw(ilyr,ibnd) * (1._wp - forwliq) / &
                                 (1._wp - forwliq * ssaliq_sw(ilyr,ibnd))
                    ssaice_del = ssaice_sw(ilyr,ibnd) * (1._wp - forwice) / &
                                 (1._wp - forwice * ssaice_sw(ilyr,ibnd))
                    scatliq = ssaliq_del * tauliq_del
                    scatice = ssaice_del * tauice_del
                  cloud_optical_props%ssacld(icol,ilyr,ibnd) = (scatliq + scatice) / &
                                 cloud_optical_props%taucld(icol,ilyr,ibnd)
                  cloud_optical_props%asycld(icol,ilyr,ibnd) = &
                       (scatice * (asyice_sw(ilyr,ibnd) - forwice) / (1._wp - forwice) + &
                        scatliq * (asyliq_sw(ilyr,ibnd) - forwliq) / (1._wp - forwliq)) / &
                       (scatice + scatliq)
                type is (ty_cloud_optical_props_nstr)
                  cloud_optical_props%taucld(icol,ilyr,ibnd) = tauliq_del + tauice_del
                    ssaliq_del = ssaliq_sw(ilyr,ibnd) * (1._wp - forwliq) / &
                                 (1._wp - forwliq * ssaliq_sw(ilyr,ibnd))
                    ssaice_del = ssaice_sw(ilyr,ibnd) * (1._wp - forwice) / &
                                 (1._wp - forwice * ssaice_sw(ilyr,ibnd))
                    scatliq = ssaliq_del * tauliq_del
                    scatice = ssaice_del * tauice_del
                  cloud_optical_props%ssacld(icol,ilyr,ibnd) = (scatliq + scatice) / &
                                 cloud_optical_props%taucld(icol,ilyr,ibnd)
                    asycld = &
                       (scatliq * (asyliq_sw(ilyr,ibnd) - forwliq) / (1._wp - forwliq) + &
                        scatice * (asyice_sw(ilyr,ibnd) - forwice) / (1._wp - forwice)) / &
                       (scatliq + scatice)
                  cloud_optical_props%pcld(1,icol,ilyr,ibnd) = 1.0_wp
                  cloud_optical_props%pcld(2,icol,ilyr,ibnd) = asycld
                  cloud_optical_props%pcld(3,icol,ilyr,ibnd) = 0.0_wp
              end select

            end select

         enddo

         endif
      enddo
   enddo

  end function cloud_optics_sw

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Ancillary functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Does this file contain variables needed to do LW calculations ? 
  !
  function is_lw(fileName) 
    character(len=*), intent(in   ) :: fileName
    logical                         :: is_lw
    
    integer :: ncid, dimid, status

    is_lw = .false.
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) & 
      call stop_on_err("is_lw: can't find file " // trim(fileName))

    if(nf90_inq_dimid(ncid, 'nband_lw', dimid) == NF90_NOERR) is_lw = .true.
          
    ncid = nf90_close(ncid) 
  end function is_lw
!  ! ----------------------
!  function is_sw(fileName) 
!    character(len=*), intent(in   ) :: fileName
!    logical                         :: is_sw
!    
!    is_sw = .not. is_lw(fileName) 
!  end function is_sw                     

  !--------------------------------------------------------------------------------------------------------------------
  subroutine stop_on_err(msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg

    if(msg /= "") then
      write(error_unit, *) msg
      stop
    end if
  end subroutine

  !--------------------------------------------------------------------------------------------------------------------
  function get_dim_length(ncid, dimname)
    !
    ! Get the length of a dimension from an open netCDF file
    !  This is unfortunately a two-step process
    !
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer :: get_dim_length

    integer :: dimid

    if(nf90_inq_dimid(ncid, trim(dimname), dimid) == NF90_NOERR) then
      if(nf90_inquire_dimension(ncid, dimid, len=get_dim_length) /= NF90_NOERR) get_dim_length = 0
    else
      get_dim_length = 0
    end if

  end function get_dim_length
  !--------------------------------------------------------------------------------------------------------------------
  function get_data_size(ncid, varName, n)
    !
    ! Returns the extents of a netcdf variable on disk
    !
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: n
    integer                      :: get_data_size(n)

    integer :: i
    integer :: varid, ndims, dimids(n)

    get_data_size(n) = -1
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("get_data_size: can't find variable " // trim(varName))
    if(nf90_inquire_variable(ncid, varid, ndims = ndims) /= NF90_NOERR) &
      call stop_on_err("get_data_size: can't get information for variable " // trim(varName))
    if(ndims /= n) &
      call stop_on_err("get_data_size:  variable " // trim(varName) // " has the wrong number of dimensions" )
    if(nf90_inquire_variable(ncid, varid, dimids = dimids) /= NF90_NOERR) &
      call stop_on_err("get_data_size: can't read dimension ids for variable " // trim(varName))
    do i = 1, n
      if(nf90_inquire_dimension(ncid, dimids(i), len = get_data_size(i)) /= NF90_NOERR) &
        call stop_on_err("get_data_size: can't get dimension lengths for variable " // trim(varName))
    end do

  end function get_data_size
  !--------------------------------------------------------------------------------------------------------------------

  function read_scalar(ncid, varName)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    real(wp)                     :: read_scalar

    integer :: varid

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_scalar)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_scalar
  !--------------------------------------------------------------------------------------------------------------------
  function read_1d_field(ncid, varName, nx)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx
    real(wp), dimension(nx)      :: read_1d_field

    integer :: varid

    if(any(get_data_size(ncid, varName, 1) /= [nx])) &
      call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent.")
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_1d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_1d_field
  !--------------------------------------------------------------------------------------------------------------------
  function read_2d_field(ncid, varName, nx, ny)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx, ny
    real(wp), dimension(nx, ny)  :: read_2d_field

    integer :: varid
    if(any(get_data_size(ncid, varName, 2) /= [nx, ny])) &
      call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent.")
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_2d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_2d_field
  !--------------------------------------------------------------------------------------------------------------------
  function read_3d_field(ncid, varName, nx, ny, nz)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx, ny, nz
    real(wp), dimension(nx, ny, nz)  :: read_3d_field

    integer :: varid

    if(any(get_data_size(ncid, varName, 3) /= [nx, ny, nz])) &
      call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent.")
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_3d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_3d_field
  !--------------------------------------------------------------------------------------------------------------------
  function read_4d_field(ncid, varName, nw, nx, ny, nz)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nw, nx, ny, nz
    real(wp), dimension(nw, nx, ny, nz)  :: read_4d_field

    integer :: varid

    if(any(get_data_size(ncid, varName, 4) /= [nw, nx, ny, nz])) &
      call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent." )
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_4d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_4d_field

  !--------------------------------------------------------------------------------------------------------------------
  function var_exists(ncid, varName) 
    !
    ! Does this variable exist (have a valid var_id) in the open netCDF file? 
    ! 
    integer,          intent(in) :: ncid 
    character(len=*), intent(in) :: varName
    logical :: var_exists
    
    integer :: varId 
    var_exists = nf90_inq_varid(ncid, trim(varName), varid) == NF90_NOERR
  end function var_exists 
  !--------------------------------------------------------------------------------------------------------------------
  function dim_exists(ncid, dimName) 
    !
    ! Does this dimension exist (have a valid dim_id) in the open netCDF file? 
    ! 
    integer,          intent(in) :: ncid 
    character(len=*), intent(in) :: dimName
    logical :: dim_exists
    
    integer :: dimid
    dim_exists = nf90_inq_dimid(ncid, trim(dimName), dimid) == NF90_NOERR
  end function dim_exists 
  !--------------------------------------------------------------------------------------------------------------------
  subroutine create_dim(ncid, dimName, dimLength) 
    !
    ! Check to see if a dimiable with this name exists in the file
    !   If so, check against current size
    !   If not, create with specified dimensions  
    ! 
    integer,          intent(in) :: ncid 
    character(len=*), intent(in) :: dimName
    integer,          intent(in) :: dimLength
    
    integer                 :: i, dimid 
    
    if(dim_exists(ncid, dimName)) then 
      if (dimLength /= get_dim_length(ncid, trim(dimName))) & 
          call stop_on_err("dim " // trim(dimName) // " is present but incorrectly sized.") 
    else 
      if(nf90_redef(ncid) /= NF90_NOERR) & 
        call stop_on_err("create_dim: can't put file into redefine mode") 
      if(nf90_def_dim(ncid, dimName, dimLength, dimid) /= NF90_NOERR) & 
        call stop_on_err("create_dim: can't define dimension " // trim(dimName))
      if(nf90_enddef(ncid) /= NF90_NOERR) & 
        call stop_on_err("create_dim: can't end redefinition??")
    end if 
  end subroutine create_dim 
  !--------------------------------------------------------------------------------------------------------------------
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
    
    integer :: i, varid, xtype
    integer :: dimIds(size(dimNames))
    
    if(var_exists(ncid, varName)) then 
      do i = 1, size(dimNames) 
        if (dimLengths(i) /= get_dim_length(ncid, trim(dimNames(i)))) & 
          call stop_on_err("Variable " // trim(varName) // " is present but incorrectly sized.") 
      end do 
    else 
      do i = 1, size(dimNames) 
        if(nf90_inq_dimid(ncid, trim(dimNames(i)), dimIds(i)) /= NF90_NOERR) & 
          call stop_on_err("create_var: Can't get id for dimension " // trim(dimnames(i)))
      end do 
      if(nf90_redef(ncid) /= NF90_NOERR) & 
        call stop_on_err("create_var: can't put file into redefine mode") 
      xtype = NF90_DOUBLE
      if(present(dataType)) xtype = dataType
      if(nf90_def_var(ncid, varName, xtype, dimIds, varid) /= NF90_NOERR) & 
        call stop_on_err("create_var: can't define variable " // trim(varName))
      if(nf90_enddef(ncid) /= NF90_NOERR) & 
        call stop_on_err("create_dim: can't end redefinition??")
    end if 
  end subroutine create_var 
  !--------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------
  ! Writing functions 
  !--------------------------------------------------------------------------------------------------------------------
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
  !--------------------------------------------------------------------------------------------------------------------
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
  !--------------------------------------------------------------------------------------------------------------------
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
  !--------------------------------------------------------------------------------------------------------------------
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
  !--------------------------------------------------------------------------------------------------------------------
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
  !--------------------------------------------------------------------------------------------------------------------
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
  !--------------------------------------------------------------------------------------------------------------------

end module mo_cloud_optics_lut
