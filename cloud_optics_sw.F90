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

