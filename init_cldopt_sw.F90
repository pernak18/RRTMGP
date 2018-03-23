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

