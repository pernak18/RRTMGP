! ----------------------------------------------------------------------------------
! Test cloud optics with cloud properties derived from Pade formulations.
! ----------------------------------------------------------------------------------
program test_cloud_optics
  use mo_rte_kind,                     only: wp
  use mo_cloud_optics_pade,            only: ty_cloud_optics_pade, is_lw
  use mo_cloud_optical_props,          only: ty_cloud_optical_props_arry, &
                                             ty_cloud_optical_props_1scl, &
                                             ty_cloud_optical_props_2str, &
                                             ty_cloud_optical_props_nstr
  use mo_load_cloud_coefficients_pade, only: load_and_init_cldop, read_cldpp, write_cldop

  implicit none 
  ! ----------------------------------------------------------------------------------
  integer :: ncol, nlay, nbnd
  integer :: nbndlw, nbndsw
  integer :: irad

  integer :: b, nBlocks, colS, colE
  integer, parameter :: blockSize = 8
  integer, parameter :: nmom = 3

  integer :: i, ibnd
  character(len=128) :: filename

  type(ty_cloud_optics_pade) :: cloud_spec
  class(ty_cloud_optical_props_arry), allocatable :: &
           cloud_optical_props, cloud_optical_props_subset

! Cloud physical property inputs and optical property outputs
  character(len=128) :: cld_io_file = 'rrtmgp-cloud-optics-inputs-outputs.nc'

  ! ----------------------------------------------------------------------------------
  ! start cloud optics test
!  print *, 'Starting cloud optics test...'
!  print *, '  Using blocks of size ', blockSize

  !
  ! Input cloud coefficient data from netCDF file
  !
  ! ----------------------------------------------------------------------------------
  call load_and_init_cldop(cloud_spec, 'cld_coefficients.nc')

  !
  ! Input cloud physical properties from NetCDF file
  !
  call read_cldpp(cloud_spec, trim(cld_io_file))

  ! Define output array sizes
  ncol = size(cloud_spec%cldfrac, 1)
  nlay = size(cloud_spec%cldfrac, 2)
  ! Pade
  nbndlw = cloud_spec%nband_lw
  nbndsw = cloud_spec%nband_sw

  ! ----------------------------------------------------------------------------------
  filename = 'cld_coefficients.nc'

  if (is_lw(trim(filename))) then 
     nbnd = nbndlw
     print*, 'Calculating longwave cloud optical depths'
     allocate(ty_cloud_optical_props_1scl::cloud_optical_props)
     allocate(ty_cloud_optical_props_1scl::cloud_optical_props_subset)
  else
     nbnd = nbndsw
     print*, 'Calculating shortwave cloud optical depths'
     allocate(ty_cloud_optical_props_2str::cloud_optical_props)
     allocate(ty_cloud_optical_props_2str::cloud_optical_props_subset)
  endif

  select type (cloud_optical_props)
    type is (ty_cloud_optical_props_1scl)       ! two-stream calculation: tau only
      call stop_on_err(cloud_optical_props%init_1scl(ncol, nlay, nbnd))
    type is (ty_cloud_optical_props_2str)       ! two-stream calculation: tau, ssa, g
      call stop_on_err(cloud_optical_props%init_2str(ncol, nlay, nbnd))
    type is (ty_cloud_optical_props_nstr)       ! n-stream calculation: tau, ssa, p
      call stop_on_err(cloud_optical_props%init_nstr(nmom, ncol, nlay, nbnd))
  end select

  !
  ! Loop over subsets of the problem 
  !
  nBlocks = ncol/blockSize ! Integer division 
  print *, "Doing ", nBlocks, "blocks of size ", blockSize
  do b = 1, nBlocks  
    colS = (b-1) * blockSize + 1 
    colE =  b    * blockSize
    call stop_on_err(cloud_optical_props%get_subset(colS, colE-colS+1, cloud_optical_props_subset))
      
    if (is_lw(trim(filename))) then 
    ! Longwave
      call stop_on_err(cloud_spec%cloud_optics(blockSize, nlay, nbndlw, is_lw(trim(filename)), &
                                 cloud_spec%cldfrac(colS:colE,:), &
                                 cloud_spec%clwp(colS:colE,:), &
                                 cloud_spec%ciwp(colS:colE,:), &
                                 cloud_spec%rel(colS:colE,:), &
                                 cloud_spec%rei(colS:colE,:), &
                                 cloud_optical_props_subset ) )
    else 
    ! Shortwave
      call stop_on_err(cloud_spec%cloud_optics(blockSize, nlay, nbndsw, &
                                 cloud_spec%cldfrac(colS:colE,:), &
                                 cloud_spec%clwp(colS:colE,:), &
                                 cloud_spec%ciwp(colS:colE,:), &
                                 cloud_spec%rel(colS:colE,:), &
                                 cloud_spec%rei(colS:colE,:), &
                                 cloud_optical_props_subset ) )
    end if
    call stop_on_err(assign_subset(cloud_optical_props_subset, colS, colE, cloud_optical_props, nbnd))
  end do 
 
  if(mod(ncol, blockSize) /= 0) then 
    colS = ncol/blockSize * blockSize + 1  ! Integer arithmetic 
    colE = ncol
    print *, "Doing ", colE-colS+1, "extra columns" 
    call stop_on_err(cloud_optical_props%get_subset(colS, colE-colS+1, cloud_optical_props_subset))

    if (is_lw(trim(filename))) then 
    ! Longwave
      call stop_on_err(cloud_spec%cloud_optics(colE-colS+1, nlay, nbndlw, is_lw(trim(filename)), &
                                 cloud_spec%cldfrac(colS:colE,:), &
                                 cloud_spec%clwp(colS:colE,:), &
                                 cloud_spec%ciwp(colS:colE,:), &
                                 cloud_spec%rel(colS:colE,:), &
                                 cloud_spec%rei(colS:colE,:), &
                                 cloud_optical_props_subset ) )
    else
    ! Shortwave
      call stop_on_err(cloud_spec%cloud_optics(colE-colS+1, nlay, nbndsw, &
                                 cloud_spec%cldfrac(colS:colE,:), &
                                 cloud_spec%clwp(colS:colE,:), &
                                 cloud_spec%ciwp(colS:colE,:), &
                                 cloud_spec%rel(colS:colE,:), &
                                 cloud_spec%rei(colS:colE,:), &
                                 cloud_optical_props_subset ) )
    end if 
    call stop_on_err(assign_subset(cloud_optical_props_subset, colS, colE, cloud_optical_props, nbnd))
  end if 

  call write_cldop(cld_io_file, is_lw(trim(filename)), cloud_optical_props)

  ! all done
  print *, 'cloud optics pade test end.'

contains
! -----------------------------------------------------------------------------------
  subroutine stop_on_err(msg)
    !
    ! Print error message and stop  
    ! 
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg
    if(len_trim(msg) > 0) then 
      write (error_unit,*) trim(msg)
      write (error_unit,*) "test_cloud_optics stopping"
      stop
    end if 
  end subroutine

! -----------------------------------------------------------------------------------
!
! Assign cloud optical properties from a set of columns to a specified position in a larger set
!   Could be bound to ty_cloud_optical_props but this would require more careful error checking and
!   it doesn't seem likely users will need this
!
  function assign_subset(subset, colS, colE, full, nbnd) result(error_msg)
    class(ty_cloud_optical_props_arry), intent(in   ) :: subset
    integer,                            intent(in   ) :: colS, colE
    class(ty_cloud_optical_props_arry), intent(inout) :: full
    integer,                            intent(in   ) :: nbnd
    character(len=128)                                :: error_msg

    real(wp), dimension(colE-colS+1, size(subset%taucld,2), size(subset%taucld,3)) :: taucld, ssacld, asycld
    real(wp), dimension(:,:,:,:), allocatable :: pcld
    integer :: nmom
    ! ------------------------------------------
    error_msg = ""
    if(colS > size(full%taucld,1) .or.  colE > size(full%taucld,1) .or. &
       colS < 1 .or.  colE <1) then
      error_msg = "  Subset, colS, colE not consistent with full cloud optical properties arrays???"
      return
    end if

    nmom = 0
    select type (subset)
      type is (ty_cloud_optical_props_nstr) ! n-stream calculation
        nmom = size(subset%pcld,1)
        allocate(pcld(nmom, colE-colS+1, nlay, nbnd))
    end select

    full%taucld(colS:colE,:,:) = subset%taucld
    ! For whatever reason the Intel compiler, at least, can't tell that full and subset
    !   have the same type, so we copy to intermediate storage.
    select type (subset)
      type is (ty_cloud_optical_props_2str) ! two-stream calculation
        ssacld = subset%ssacld
        asycld = subset%asycld
      type is (ty_cloud_optical_props_nstr) ! n-stream calculation
        ssacld = subset%ssacld
        pcld = subset%pcld
    end select
    select type (full)
      type is (ty_cloud_optical_props_2str) ! two-stream calculation
        full%ssacld(colS:colE,:,:) = ssacld
        full%asycld(colS:colE,:,:) = asycld
      type is (ty_cloud_optical_props_nstr) ! n-stream calculation
        full%ssacld(colS:colE,:,:) = ssacld
        full%pcld(1:nmom,colS:colE,:,:) = pcld
    end select

  end function assign_subset
! -----------------------------------------------------------------------------------

end program test_cloud_optics
