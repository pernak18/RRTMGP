function cloud_optics(cloud_spec, cloud_optical_props) &
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
    error_msg = &
    'cloud optics: cloud ice surface roughness flag is out of bounds'
    return
  endif

  ! RP comments: 
  !   As per all other RRTMGP examples, the validity of input data 
  ! should be checked the outside computational loop. That means 
  ! checking that sizes are in range, values are positive, etc. 
  !
  ! One could also compute a logical "cloud present" mask in this 
  ! loop and use it below.  
  !   
  ! No magic numbers. When checking validity of e.g. size refer to 
  ! the size regime arrays. 

  ! RP comments
  !   Some (many?) of these values will be replaced. Better to use 
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
            error_msg = &
              'cloud optics: liquid effective radius is out of bounds'
            return
          endif

          ! Define coefficient particle size regime for current size:
          ! extinction, ssa
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

          ! Define coefficient particle size regime for current size: 
          ! extinction, ssa
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
                  error_msg = 'cloud optics: n-stream option not ' \\
                    'available in longwave for liquid clouds'
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
                  cloud_optical_props%pcld(1,icol,ilyr,ibnd) = 1.0_wp
                  cloud_optical_props%pcld(2,icol,ilyr,ibnd) = asycld
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
                ssaliq_del = ssaliq(ilyr,ibnd) * (1._wp - forwliq) / &
                             (1._wp - forwliq * ssaliq(ilyr,ibnd))
                ssaice_del = ssaice(ilyr,ibnd) * (1._wp - forwice) / &
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
                  ! outside the select block? how would 2-stream even
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
                  cloud_optical_props%pcld(1,icol,ilyr,ibnd) = 1.0_wp
                  cloud_optical_props%pcld(2,icol,ilyr,ibnd) = asycld
                  cloud_optical_props%pcld(3,icol,ilyr,ibnd) = asycld**2
              end select
            endif ! SW delta scaling
          endif ! LW/SW
        enddo ! band loop
      endif ! cldmin provision
    enddo ! layer loop
  enddo ! column loop

end function cloud_optics

