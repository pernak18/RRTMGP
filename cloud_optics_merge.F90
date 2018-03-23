function cloud_optics_lw(cloud_spec, cloud_optical_props) &
  result(error_msg)
  ! ------------------------------------------------------------------------------
  ! Purpose:  Compute the cloud optical properties for each cloudy layer.
  ! ------- Input -------

  ! cloud specification data
  class(ty_optics_1scl), intent(inout) :: cloud_spec

  ! ------- Output -------
  ! Dimensions: (ncol x nlayers x nbndlw)
  class(ty_cloud_optical_props_arry), intent(inout) :: &
    cloud_optical_props

  ! ------- Local -------

  character(len=128)    :: error_msg

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
  real(wp) :: tauliq_lw, tauice_lw

  ! ice and liquid scattering terms, (liquid + ice?) asymmetry
  real(wp) :: scatice, scatliq, asycld

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
    enddo ! End layer loop
  enddo ! End column loop

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

