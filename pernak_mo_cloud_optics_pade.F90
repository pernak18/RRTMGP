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
! and return cloud optical properties by band using Pade formulations. 
!

module mo_cloud_optics_pade
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

  !----------------------------------------------------------------------------------------
!
! RP comments: 
!   It's true that ty_gas_optics extends ty_spectral_disc, but this is unique because the gas optics 
!   actually sets the spectral discretization. Instead we should extend ty_optical_props_1scl, _2str, and _nstr 
!   from mo_optical_props. 
!   I don't see a way to avoid extending each with the same data components.
!   We will need a load function to add the cloud optics data (coefficients or lookup tables); 
!   only one function is needed and all the types can point to it. 
!   We will need a compute function that fills in the optical properties arrays (tau/omega0/g or whatever)
!   validation() will need to be extended to ensure that the cloud optics data are present,  that 
!   cloud physical and optical properties arrays are allocated, the same size as the tau/omega0/g arrays, ... 
!   This code can be used generally -- it's just optical properties by band from data -- so maybe call it 
!   "ty_optics_1scl" etc.? 
!  


  type, extends(ty_spectral_disc), public :: ty_cloud_optics_pade
!    private

    ! User input
    ! Ice surface roughness category - needed for Yang (2013) ice optics parameterization
    integer :: icergh                                   ! (1 = none, 2 = medium, 3 = high)

! RP comments: 
!   Delta scaling should be omitted here. It's available from ty_optical_props_2scl etc 
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

! RP comments: 
!  Cloud optics types shouldn't have hardcoded LW/SW distinctions. Better to have 
!  explicit wavelength or wavenumber limits on bands provided at initialization. 
!  compute() should include a variable of ty_spectral_disc as an input argument. If 
!  the cloud optics variable doesn't know how to provide data on that spectral grid 
!  an error string should be returned. 
!  It would also be possible to have a single type that can provide both LW and SW values. 

    ! Model input: LUT
    ! Pade coefficient dimensions
    integer :: nband_lw,        &
               nband_sw,        &
               nrghice
    integer :: nsizereg,        & 
               ncoeff_ext,      &
               ncoeff_ssa_g
               
! RP comments: 
!   Magic numbers are undesirable
!   Size regimes should be parameterized here: number of regimes, size limits. 
!   
               
    ! Particle size regimes for Pade formulations
    integer, dimension(2,4) :: pade_sizreg_liqlw, pade_sizreg_icelw
    integer, dimension(2,4) :: pade_sizreg_liqsw, pade_sizreg_icesw
    ! Pade cloud coefficients
    real(wp), dimension(:,:,:  ), allocatable :: pade_extliq
    real(wp), dimension(:,:,:  ), allocatable :: pade_ssaliq
    real(wp), dimension(:,:,:  ), allocatable :: pade_asyliq
    real(wp), dimension(:,:,:,:), allocatable :: pade_extice
    real(wp), dimension(:,:,:,:), allocatable :: pade_ssaice
    real(wp), dimension(:,:,:,:), allocatable :: pade_asyice

! ------------------------------------------------------------------------------------------
  contains
! RP comments: 
!  Cloud optics types shouldn't have hardcoded LW/SW distinctions. 
    generic, public :: init_cldopt  => init_cldopt_lw, init_cldopt_sw
    generic, public :: cloud_optics => cloud_optics_lw, cloud_optics_sw
    ! Internal procedures
    procedure, private :: init_cldopt_lw
    procedure, private :: init_cldopt_sw
    procedure, private :: cloud_optics_lw
    procedure, private :: cloud_optics_sw
    procedure, private :: get_irad
    procedure, private :: pade_ext
    procedure, private :: pade_ssa
    procedure, private :: pade_asy

  end type ty_cloud_optics_pade

  contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Cloud optics initialization 
  !
  !--------------------------------------------------------------------------------------------------------------------
  ! Initialize object based on data read from netCDF file however the user desires. 
  !   Pade coefficient data used for longwave cloud property conversion
  ! 
  function init_cldopt_lw(cloud_spec, cld_coeff_file, is_lw) result(error_msg)

    ! Pade
    class(ty_cloud_optics_pade), intent(inout) :: cloud_spec

    ! Cloud coefficient optical property input file - Pade
    character(len=*), intent(in) :: cld_coeff_file

    logical, intent(in) :: is_lw

    character(len = 128) error_msg

    ! -----------------
    ! Local variables
    integer :: ncid
    integer :: nband_lw
    integer :: nrghice
    integer :: nsizereg
    integer :: ncoeff_ext
    integer :: ncoeff_ssa_g
! RP comments: 
!  The limits of the size regimes are data and should be provided the same way the Pade coefficients are. 
!  A axis of dimension 2 is confusing; better to have different variables (pade_sizreg_liq_g, pade_sizreg_liq_tau)
!  A possible refinement: having the same size regimes for 
 
    ! Particle size regimes for Pade formulations
    cloud_spec%pade_sizreg_liqlw(1,:) = (/2,10,35,60/)         ! for ext, ssa
    cloud_spec%pade_sizreg_liqlw(2,:) = (/2,8,20,60/)          ! for asy
    cloud_spec%pade_sizreg_icelw(1,:) = (/10,20,30,180/)       ! for all, except ssa, asy when icergh=0.00
    cloud_spec%pade_sizreg_icelw(2,:) = (/10,20,50,180/)       ! for ssa, asy, icergh=0.00

    ! -----------------
! RP comments: 
!  Types should be initialzed with data directly. No I/O in code we will distribute. 
!  See the initialization of gas_optics

    error_msg = ""

    ! Open cloud optical property coefficient file
    if(nf90_open(trim(cld_coeff_file), NF90_WRITE, ncid) /= NF90_NOERR) & 
       call stop_on_err("init_cldopt_lw(): can't open file " // trim(cld_coeff_file))

    cloud_spec%nband_lw     = get_dim_length(ncid,'nband_lw') 
    cloud_spec%nrghice      = get_dim_length(ncid,'nrghice') 
    cloud_spec%nsizereg     = get_dim_length(ncid,'nsizereg') 
    cloud_spec%ncoeff_ext   = get_dim_length(ncid,'ncoeff_ext') 
    cloud_spec%ncoeff_ssa_g = get_dim_length(ncid,'ncoeff_ssa_g') 

    nband_lw     = cloud_spec%nband_lw
    nrghice      = cloud_spec%nrghice
    nsizereg     = cloud_spec%nsizereg
    ncoeff_ext   = cloud_spec%ncoeff_ext
    ncoeff_ssa_g = cloud_spec%ncoeff_ssa_g

    ! Allocate cloud property Pade coefficient input arrays
    allocate(cloud_spec%pade_extliq(nband_lw, nsizereg, ncoeff_ext))
    allocate(cloud_spec%pade_ssaliq(nband_lw, nsizereg, ncoeff_ssa_g))
    allocate(cloud_spec%pade_asyliq(nband_lw, nsizereg, ncoeff_ssa_g))
    allocate(cloud_spec%pade_extice(nband_lw, nsizereg, ncoeff_ext, nrghice))
    allocate(cloud_spec%pade_ssaice(nband_lw, nsizereg, ncoeff_ssa_g, nrghice))
    allocate(cloud_spec%pade_asyice(nband_lw, nsizereg, ncoeff_ssa_g, nrghice))

    cloud_spec%pade_extliq  = read_field(ncid, 'pade_extliq_lw', nband_lw, nsizereg, ncoeff_ext) 
    cloud_spec%pade_ssaliq  = read_field(ncid, 'pade_ssaliq_lw', nband_lw, nsizereg, ncoeff_ssa_g) 
    cloud_spec%pade_asyliq  = read_field(ncid, 'pade_asyliq_lw', nband_lw, nsizereg, ncoeff_ssa_g) 
    cloud_spec%pade_extice  = read_field(ncid, 'pade_extice_lw', nband_lw, nsizereg, ncoeff_ext, nrghice) 
    cloud_spec%pade_ssaice  = read_field(ncid, 'pade_ssaice_lw', nband_lw, nsizereg, ncoeff_ssa_g, nrghice) 
    cloud_spec%pade_asyice  = read_field(ncid, 'pade_asyice_lw', nband_lw, nsizereg, ncoeff_ssa_g, nrghice) 

    ncid = nf90_close(ncid) 

  end function init_cldopt_lw

  !--------------------------------------------------------------------------------------------------------------------
  ! Initialize object based on data read from netCDF file however the user desires. 
  !   Pade coefficient data used for shortwave cloud property conversion
  ! 
  function init_cldopt_sw(cloud_spec, cld_coeff_file) result(error_msg)

    ! Pade
    class(ty_cloud_optics_pade), intent(inout) :: cloud_spec

    ! Cloud coefficient optical property input file - Pade
    character(len=*), intent(in) :: cld_coeff_file

    character(len = 128) error_msg

    ! -----------------
    ! Local variables
    integer :: ncid
    integer :: nband_sw
    integer :: nrghice
    integer :: nsizereg
    integer :: ncoeff_ext
    integer :: ncoeff_ssa_g

    ! Particle size regimes for Pade formulations
    cloud_spec%pade_sizreg_liqsw(1,:) = (/2,10,35,60/)         ! for ext, ssa
    cloud_spec%pade_sizreg_liqsw(2,:) = (/2,9,20,60/)          ! for asy
    cloud_spec%pade_sizreg_icesw(1,:) = (/10,20,30,180/)       ! for ext, ssa
    cloud_spec%pade_sizreg_icesw(2,:) = (/10,20,30,180/)       ! for asy

    ! -----------------
    error_msg = ""

    ! Open cloud optical property coefficient file
    if(nf90_open(trim(cld_coeff_file), NF90_WRITE, ncid) /= NF90_NOERR) & 
       call stop_on_err("init_cldopt_sw(): can't open file " // trim(cld_coeff_file))

    cloud_spec%nband_sw     = get_dim_length(ncid,'nband_sw') 
    cloud_spec%nrghice      = get_dim_length(ncid,'nrghice') 
    cloud_spec%nsizereg     = get_dim_length(ncid,'nsizereg') 
    cloud_spec%ncoeff_ext   = get_dim_length(ncid,'ncoeff_ext') 
    cloud_spec%ncoeff_ssa_g = get_dim_length(ncid,'ncoeff_ssa_g') 

    nband_sw     = cloud_spec%nband_sw
    nrghice      = cloud_spec%nrghice
    nsizereg     = cloud_spec%nsizereg
    ncoeff_ext   = cloud_spec%ncoeff_ext
    ncoeff_ssa_g = cloud_spec%ncoeff_ssa_g

    ! Allocate cloud property Pade coefficient input arrays
    allocate(cloud_spec%pade_extliq(nband_sw, nsizereg, ncoeff_ext))
    allocate(cloud_spec%pade_ssaliq(nband_sw, nsizereg, ncoeff_ssa_g))
    allocate(cloud_spec%pade_asyliq(nband_sw, nsizereg, ncoeff_ssa_g))
    allocate(cloud_spec%pade_extice(nband_sw, nsizereg, ncoeff_ext, nrghice))
    allocate(cloud_spec%pade_ssaice(nband_sw, nsizereg, ncoeff_ssa_g, nrghice))
    allocate(cloud_spec%pade_asyice(nband_sw, nsizereg, ncoeff_ssa_g, nrghice))

    cloud_spec%pade_extliq  = read_field(ncid, 'pade_extliq_sw', nband_sw, nsizereg, ncoeff_ext) 
    cloud_spec%pade_ssaliq  = read_field(ncid, 'pade_ssaliq_sw', nband_sw, nsizereg, ncoeff_ssa_g) 
    cloud_spec%pade_asyliq  = read_field(ncid, 'pade_asyliq_sw', nband_sw, nsizereg, ncoeff_ssa_g) 
    cloud_spec%pade_extice  = read_field(ncid, 'pade_extice_sw', nband_sw, nsizereg, ncoeff_ext, nrghice) 
    cloud_spec%pade_ssaice  = read_field(ncid, 'pade_ssaice_sw', nband_sw, nsizereg, ncoeff_ssa_g, nrghice) 
    cloud_spec%pade_asyice  = read_field(ncid, 'pade_asyice_sw', nband_sw, nsizereg, ncoeff_ssa_g, nrghice) 

    ncid = nf90_close(ncid) 

  end function init_cldopt_sw

  ! ------------------------------------------------------------------------------
  !
  ! Derive LW cloud optical properties from provided cloud physical pproperties
  !
  ! ------------------------------------------------------------------------------
  function cloud_optics_lw( &
  ! Input
! RP comments: the cloud physical properties and their dimensions should either be 
!   part of the type data (they're included above) or they should be arguments to this 
!   routine, but not both. 
!   Use "this" convection as elsewhere in RRTMGP to refer to the object being called. 
!   If ty_optics_1scl descends from ty_optical_props_1scl there won't be an output argument - 
!     the optical properties will go in this%tau etc. 

! Outstanding question: compute, increment vs. directly calling increment() ? 
!   Maybe that's an argument for passing arguments to routines. 
!   

                   cloud_spec, &
                   ncol, nlayers, nbndlw, is_lw, &
                   cldfrac, clwp, ciwp, rel, rei, &
  ! Output
                   cloud_optical_props) &
                   result(error_msg)
  ! ------------------------------------------------------------------------------

  ! Purpose:  Compute the cloud optical properties for each cloudy layer.

  ! ------- Input -------

      class(ty_cloud_optics_pade), intent(inout) :: cloud_spec
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

      real(wp) :: p_e(6)                         ! Local pade coefficients for extinction
      real(wp) :: p(5)                           ! Local pade coefficients for ssa and g

      real(wp) :: cwp                            ! cloud water path (liquid + ice)
      real(wp), parameter :: cldmin = 1.e-20     ! minimum value for cloud quantities
      real(wp) :: radliq                         ! cloud liquid droplet radius (microns)
      real(wp) :: radice                         ! cloud ice effective size (microns)
      real(wp) :: factor                         ! 
      real(wp) :: fint                           ! 

      integer :: index                           !
      integer :: icol, ilyr, ibnd                !
      integer :: irad, irade, irads, iradg       !
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

! RP comments: 
!   As per all other RRTMGP examples, the validity of input data should be checked 
!   outside the computational loop. That means checking that sizes are in range, 
!   values are positive, etc. 
!   One could also compute a logical "cloud present" mask in this loop and use it below.  
!   
!   No magic numbers. When checking validity of e.g. size refer to the size regime arrays. 

! RP comments
!   Some (many?) of these values will be replaced. Better to use merge() in loops below. 
!
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

! Cloud optical properties from Pade coefficients

      do ilyr = 1, nlayers
         cwp = ciwp(icol,ilyr) + clwp(icol,ilyr)
         if (cldfrac(icol,ilyr) .gt. cldmin .and. cwp .gt. cldmin) then 
! Derive optical properties from Pade functions 
! Original Pade formulations for EXT, SSA, G

! Liquid - LW
            radliq = rel(icol,ilyr)
            if (radliq .gt. 0.0_wp .and. clwp(icol,ilyr) .gt. 0.0_wp) then 
! For liquid OP, particle size is limited to 2.5 to 20.0 microns
               if (radliq .lt. 2.5_wp .or. radliq .gt. 20.0_wp) then 
                  error_msg = 'cloud optics: liquid effective radius is out of bounds'
                  return
               endif
! Define coefficient particle size regime for current size: extinction, ssa
               irade = cloud_spec%get_irad(radliq, 'liq', 'lw', 'ext')
               irads = cloud_spec%get_irad(radliq, 'liq', 'lw', 'ssa')
               iradg = cloud_spec%get_irad(radliq, 'liq', 'lw', 'asy')

!
! RP comments: copying coefficents is inefficient. Why not use coefficients directly? 
!   extliq_lw(ilyr,ibnd) = cloud_spec%pade_ext(cloud_spec%pade_extliq(ibnd,irade,:), radliq)

               do ibnd = 1, nbndlw
                  p_e(:) = cloud_spec%pade_extliq(ibnd,irade,:)
                  extliq_lw(ilyr,ibnd) = cloud_spec%pade_ext(p_e, radliq)
                  p(:) = cloud_spec%pade_ssaliq(ibnd,irads,:)
                  ssaliq_lw(ilyr,ibnd) = cloud_spec%pade_ssa(p, radliq)
                  p(:) = cloud_spec%pade_asyliq(ibnd,iradg,:)
                  asyliq_lw(ilyr,ibnd) = cloud_spec%pade_asy(p, radliq)
               enddo
            endif

! Ice - LW
            radice = rei(icol,ilyr)
            if (radice .gt. 0.0_wp .and. ciwp(icol,ilyr) .gt. 0.0_wp) then 
! For Yang (2013) ice OP, particle size is limited to 10.0 to 180.0 microns
               if (radice .lt. 10.0_wp .or. radice .gt. 180.0_wp) then
                  error_msg = 'cloud optics: ice effective radius is out of bounds'
                  return
               endif
! Define coefficient particle size regime for current size: extinction, ssa
               irade = cloud_spec%get_irad(radice, 'ice', 'lw', 'ext')
               irads = cloud_spec%get_irad(radice, 'ice', 'lw', 'ssa')
               iradg = cloud_spec%get_irad(radice, 'ice', 'lw', 'asy')


               do ibnd = 1, nbndlw
! Derive optical properties for selected ice roughness

! RP comments: 
!  As far as I can tell the case() statement here doesn't have an effect: 
!   the coefficients are determined by the value of icergh
                  select case (icergh)

! No ice roughness 
                  case(1)
                     p_e(:) = cloud_spec%pade_extice(ibnd,irade,:,icergh)
                     extice_lw(ilyr,ibnd) = cloud_spec%pade_ext(p_e, radice)
                     p(:) = cloud_spec%pade_ssaice(ibnd,irads,:,icergh)
                     ssaice_lw(ilyr,ibnd) = cloud_spec%pade_ssa(p, radice)
                     p(:) = cloud_spec%pade_asyice(ibnd,iradg,:,icergh)
                     asyice_lw(ilyr,ibnd) = cloud_spec%pade_asy(p, radice)

! Medium ice roughness 
                  case(2)
                     p_e(:) = cloud_spec%pade_extice(ibnd,irade,:,icergh)
                     extice_lw(ilyr,ibnd) = cloud_spec%pade_ext(p_e, radice)
                     p(:) = cloud_spec%pade_ssaice(ibnd,irads,:,icergh)
                     ssaice_lw(ilyr,ibnd) = cloud_spec%pade_ssa(p, radice)
                     p(:) = cloud_spec%pade_asyice(ibnd,iradg,:,icergh)
                     asyice_lw(ilyr,ibnd) = cloud_spec%pade_asy(p, radice)

! High ice roughness 
                  case(3)
                     p_e(:) = cloud_spec%pade_extice(ibnd,irade,:,icergh)
                     extice_lw(ilyr,ibnd) = cloud_spec%pade_ext(p_e, radice)
                     p(:) = cloud_spec%pade_ssaice(ibnd,irads,:,icergh)
                     ssaice_lw(ilyr,ibnd) = cloud_spec%pade_ssa(p, radice)
                     p(:) = cloud_spec%pade_asyice(ibnd,iradg,:,icergh)
                     asyice_lw(ilyr,ibnd) = cloud_spec%pade_asy(p, radice)

                  end select
               enddo
            endif
         endif
! End layer loop
      enddo

! End column loop
   enddo

! Combine liquid and ice contributions into total cloud optical properties
   do icol = 1, ncol
      do ilyr = 1, nlayers
         cwp = ciwp(icol,ilyr) + clwp(icol,ilyr)
         if (cldfrac(icol,ilyr) .gt. cldmin .and. cwp .gt. cldmin) then 
! Longwave
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
  ! Derive SW cloud optical properties from provided cloud physical pproperties
  !
  ! ------------------------------------------------------------------------------
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

      class(ty_cloud_optics_pade), intent(inout) :: cloud_spec
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

      real(wp) :: p_e(6)                         ! Local pade coefficients for extinction
      real(wp) :: p(5)                           ! Local pade coefficients for ssa and g

      real(wp) :: cwp                            ! cloud water path (liquid + ice)
      real(wp), parameter :: cldmin = 1.e-20     ! minimum value for cloud quantities
      real(wp) :: radliq                         ! cloud liquid droplet radius (microns)
      real(wp) :: radice                         ! cloud ice effective size (microns)
      real(wp) :: factor                         ! 
      real(wp) :: fint                           ! 

      integer :: index                           !
      integer :: icol, ilyr, ibnd                !
      integer :: irad, irade, irads, iradg       !
      integer :: icergh                          ! ice surface roughness
                                                 ! (1 = none, 2 = medium, 3 = high)
      integer :: idelscl                         ! delta-scaling of output
                                                 ! (0 = no delta-scaling, 1 = with delta-scaling)

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

! Cloud optical properties from Pade coefficients

      do ilyr = 1, nlayers
         cwp = ciwp(icol,ilyr) + clwp(icol,ilyr)
         if (cldfrac(icol,ilyr) .gt. cldmin .and. cwp .gt. cldmin) then 
! Derive optical properties from Pade functions 
! Original Pade formulations for EXT, SSA, G

! Liquid OP
            radliq = rel(icol,ilyr)
            if (radliq .gt. 0.0_wp .and. clwp(icol,ilyr) .gt. 0.0_wp) then 
! Define coefficient particle size regime for current size: extinction, ssa
               irade = cloud_spec%get_irad(radliq, 'liq', 'sw', 'ext')
               irads = cloud_spec%get_irad(radliq, 'liq', 'sw', 'ssa')
               iradg = cloud_spec%get_irad(radliq, 'liq', 'sw', 'asy')
               do ibnd = 1, nbndsw
                  p_e(:) = cloud_spec%pade_extliq(ibnd,irade,:)
                  extliq_sw(ilyr,ibnd) = cloud_spec%pade_ext(p_e, radliq)
                  p(:) = cloud_spec%pade_ssaliq(ibnd,irads,:)
                  ssaliq_sw(ilyr,ibnd) = cloud_spec%pade_ssa(p, radliq)
                  p(:) = cloud_spec%pade_asyliq(ibnd,iradg,:)
                  asyliq_sw(ilyr,ibnd) = cloud_spec%pade_asy(p, radliq)
               enddo
            endif

! Ice OP for requested ice roughness (icergh)
            radice = rei(icol,ilyr)
            if (radice .gt. 0.0_wp .and. ciwp(icol,ilyr) .gt. 0.0_wp) then 
! Define coefficient particle size regime for current size: extinction, ssa
               irade = cloud_spec%get_irad(radice, 'ice', 'sw', 'ext')
               irads = cloud_spec%get_irad(radice, 'ice', 'sw', 'ssa')
               iradg = cloud_spec%get_irad(radice, 'ice', 'sw', 'asy')
               do ibnd = 1, nbndsw
! Derive optical properties for selected ice roughness
                  select case (icergh)

! No ice roughness 
                  case(1)
                     p_e(:) = cloud_spec%pade_extice(ibnd,irade,:,icergh)
                     extice_sw(ilyr,ibnd) = cloud_spec%pade_ext(p_e, radice)
                     p(:) = cloud_spec%pade_ssaice(ibnd,irads,:,icergh)
                     ssaice_sw(ilyr,ibnd) = cloud_spec%pade_ssa(p, radice)
                     p(:) = cloud_spec%pade_asyice(ibnd,iradg,:,icergh)
                     asyice_sw(ilyr,ibnd) = cloud_spec%pade_asy(p, radice)

! Medium ice roughness 
                  case(2)
                     p_e(:) = cloud_spec%pade_extice(ibnd,irade,:,icergh)
                     extice_sw(ilyr,ibnd) = cloud_spec%pade_ext(p_e, radice)
                     p(:) = cloud_spec%pade_ssaice(ibnd,irads,:,icergh)
                     ssaice_sw(ilyr,ibnd) = cloud_spec%pade_ssa(p, radice)
                     p(:) = cloud_spec%pade_asyice(ibnd,iradg,:,icergh)
                     asyice_sw(ilyr,ibnd) = cloud_spec%pade_asy(p, radice)

! High ice roughness 
                  case(3)
                     p_e(:) = cloud_spec%pade_extice(ibnd,irade,:,icergh)
                     extice_sw(ilyr,ibnd) = cloud_spec%pade_ext(p_e, radice)
                     p(:) = cloud_spec%pade_ssaice(ibnd,irads,:,icergh)
                     ssaice_sw(ilyr,ibnd) = cloud_spec%pade_ssa(p, radice)
                     p(:) = cloud_spec%pade_asyice(ibnd,iradg,:,icergh)
                     asyice_sw(ilyr,ibnd) = cloud_spec%pade_asy(p, radice)

                  end select
               enddo

            endif
         endif
! End layer loop
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
                  cloud_optical_props%ssacld(icol,ilyr,ibnd) = (scatice + scatliq) / &
                                        cloud_optical_props%taucld(icol,ilyr,ibnd)
                  asycld = &
                       (scatice * asyice_sw(ilyr,ibnd) + scatliq * asyliq_sw(ilyr,ibnd)) / &
                       (scatice + scatliq)
                  cloud_optical_props%pcld(1,icol,ilyr,ibnd) = 1.0_wp
                  cloud_optical_props%pcld(2,icol,ilyr,ibnd) = asycld
                  cloud_optical_props%pcld(3,icol,ilyr,ibnd) = asycld**2
              end select

! RP comments: 
!   Delta scaling should be omitted here. It's available from ty_optical_props_2scl etc 

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
  !---------------------------------------------------------------------------
  function get_irad(cloud_spec,rad,phase,regime,param)
     class(ty_cloud_optics_pade), intent(inout) :: cloud_spec
     real(wp), intent(in) :: rad             ! particle radius
     character(len=3), intent(in) :: phase   ! liq/ice
     character(len=2), intent(in) :: regime  ! lw/sw
     character(len=3), intent(in) :: param   ! ext/ssa/asy
     integer :: get_irad                     ! irad index
     ! Local variables
     integer :: irad
     real(wp), dimension(4) :: sizreg

! Liq/LW
     if (phase .eq. 'liq' .and. regime .eq. 'lw') then
        if (param .eq. 'ext' .or. param .eq. 'ssa') sizreg = cloud_spec%pade_sizreg_liqlw(1,:)
        if (param .eq. 'asy') sizreg = cloud_spec%pade_sizreg_liqlw(2,:)
        do irad = 1, 2 
           if (rad .gt. sizreg(irad) .and. rad .le. sizreg(irad+1)) get_irad = irad
        enddo
     endif

! Ice/LW
     if (phase .eq. 'ice' .and. regime .eq. 'lw') then
        sizreg = cloud_spec%pade_sizreg_icelw(1,:)
        if (cloud_spec%icergh .eq. 1 .and. param .eq. 'ssa' .or. param .eq. 'asy') &
           sizreg = cloud_spec%pade_sizreg_icelw(2,:)
           do irad = 1, 3 
              if (rad .gt. sizreg(irad) .and. rad .le. sizreg(irad+1)) get_irad = irad
           enddo
        endif

! Liq/SW - ext, ssa
     if (phase .eq. 'liq' .and. regime .eq. 'sw') then
        if (param .eq. 'ext' .or. param .eq. 'ssa') sizreg = cloud_spec%pade_sizreg_liqsw(1,:)
        if (param .eq. 'asy') sizreg = cloud_spec%pade_sizreg_liqsw(2,:)
        do irad = 1, 2 
           if (rad .gt. sizreg(irad) .and. rad .le. sizreg(irad+1)) get_irad = irad
        enddo
     endif

! Ice/SW
     if (phase .eq. 'ice' .and. regime .eq. 'sw') then
        sizreg = cloud_spec%pade_sizreg_icesw(1,:)
        do irad = 1, 3 
           if (rad .gt. sizreg(irad) .and. rad .le. sizreg(irad+1)) get_irad = irad
        enddo
     endif

  end function get_irad

! RP comments: 
!  Much better to replace repeated code with a single routine to evaluate Pade approximants
!  since the formula is entirely general https://en.wikipedia.org/wiki/Padé_approximant

!   function pade(n_num, num, n_den, denom)
!     integer, intent(in) :: n_num, n_den
!     real(wp), intent(in) :: num(n_num), den(n_den) 
!
!  Direct exponentiation is very expensive. Please evaluate polynomals with Horner 
!  (p(1) + p(2)*reff + p(3)*reff**2  = p(1) + reff*(p(2) + reff * p3))
!  or better as a loop over the number of terms 

! If ssa is fit as coalbedo (1 - ssa) that can be done where this routine is called. 

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

end module mo_cloud_optics_pade
