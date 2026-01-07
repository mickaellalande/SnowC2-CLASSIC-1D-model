!> \file
!! Calculates the N litter generated from leaves, stem, and root.
!> @author A. Asaadi
!!

subroutine nlitter(il1, il2, sort, fcancmx, lfstatus, iday, & !In
                   tltrleaf, tltrstem, tltrroot, gl2bl_grass_cflux, & !In
                   bl2ltr_grass_cflux, gleafmas, bleafmas, stemmass, & !In
                   rootmass, & !In
                   ngleafmas, ngleafmas_ns, ngleafmas_s, nbleafmas, & !In/Out
                   nstemmass, nstemmass_ns, nstemmass_s, nrootmass, & !In/Out
                   nrootmass_ns, nrootmass_s, nlitrmass, & ! In/Out
                   gl2bl_grass_nflux, nlitrveg_l, nlitrveg_s, nlitrveg_r) !Out

  use classicParams, only : ilg, ignd, icc, iccp1, ctempfts, nresorp_coeff, &
                            ntolerance, c2n_lmin, deltat

  implicit none

  integer,intent(in) :: il1      !< il1=1
  integer,intent(in) :: il2      !< il2=ilg (no. of grid cells in latitude circle)
  integer,intent(in) :: iday                !< day of year
  integer,intent(in),dimension(ilg,icc) :: lfstatus  !< integer indicating leaf status or mode (max. growth=1, normal growth=2, fall=3, no leaves=4)
  integer,intent(in),dimension(icc) :: sort          !< index for correspondence between PFTs and the 12 values in parameters vectors
  real,intent(in),dimension(ilg,icc)  :: fcancmx            !< max. fractional coverage of PFTs
  real,intent(in),dimension(ilg,icc)  :: tltrleaf           !< total leaf carbon litter fall rate, \f$umol CO2 m^{-2} sec^{-1}\f$ 
  real,intent(in),dimension(ilg,icc)  :: tltrstem           !< total stem carbon litter fall rate, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in),dimension(ilg,icc)  :: tltrroot           !< total root carbon litter fall rate, \f$umol CO2 m^{-2} sec^{-1}\f$
  real,intent(in),dimension(ilg,icc)  :: gleafmas           !< green leaf carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in),dimension(ilg,icc)  :: bleafmas           !< brown leaf carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in),dimension(ilg,icc)  :: stemmass           !< stem carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in),dimension(ilg,icc)  :: rootmass           !< root carbon for each PFT, \f$kg C m^{-2}\f$
  real,intent(in),dimension(ilg,icc)  :: gl2bl_grass_cflux  !< carbon flux from green to brown leaf carbon, \f$kg C m^{-2} day^{-1}\f$ 
  real,intent(in),dimension(ilg,icc)  :: bl2ltr_grass_cflux !< carbon flux from brown leaf carbon to litter carbon, \f$kg C m^{-2} day^{-1}\f$ 

  real,intent(out),dimension(ilg,icc) :: nlitrveg_l         !< total leaf nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nlitrveg_s         !< total stem nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: nlitrveg_r         !< total root nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,icc) :: gl2bl_grass_nflux  !< nitrogen flux from green to brown leaf nitrogen, \f$g N m^{-2} day^{-1}\f$

  real,intent(inout),dimension(ilg,icc) :: ngleafmas        !< green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: ngleafmas_ns     !< non-structural green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: ngleafmas_s      !< structural green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nbleafmas        !< brown leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nstemmass        !< stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nstemmass_ns     !< non-structural stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nstemmass_s      !< structural stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nrootmass        !< root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nrootmass_ns     !< non-structural root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,icc) :: nrootmass_s      !< structural root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,iccp1) :: nlitrmass      !< litter nitrogen for individual PFTs + bare soil, \f$g N m^{-2}\f$

  ! local variables
  integer :: i, j
  character(8) :: pftkind
  real :: temp_leaf    !< temporary tltrleaf in \f$kg C m^{-2}\f$
  real :: temp_stem    !< temporary tltrstem in \f$kg C m^{-2}\f$
  real :: temp_root    !< temporary tltrroot in \f$kg C m^{-2}\f$
  real :: gleafltr_fr  !< fraction of gleafmas becoming litter 
  real :: stemltr_fr   !< fraction of stemmass becoming litter 
  real :: rootltr_fr   !< fraction of rootmass becoming litter 
  real,dimension(ilg,icc) :: ngleafmas_ns_in  !< previous non-structural green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: ngleafmas_s_in   !< previous structural green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: ngleafmas_in     !< previous green leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nbleafmas_in     !< previous brown leaf nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nstemmass_ns_in  !< previous non-structural stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nstemmass_s_in   !< previous structural stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nstemmass_in     !< previous stem nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nrootmass_ns_in  !< previous non-structural root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nrootmass_s_in   !< previous structural root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nrootmass_in     !< previous root nitrogen for each PFT, \f$g N m^{-2}\f$
  real,dimension(ilg,iccp1) :: nlitrmass_in   !< previous litter nitrogen for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,dimension(ilg,icc) :: nlitrveg_l_ns      !< non-structural leaf nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: nlitrveg_l_s       !< structural leaf nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: gl2bl_grass_nflux_ns  !< non-structural nitrogen flux from green to brown leaf nitrogen, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: gl2bl_grass_nflux_s   !< structural nitrogen flux from green to brown leaf nitrogen, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: nlitrveg_s_ns      !< non-structural stem nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: nlitrveg_s_s       !< structural stem nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: nlitrveg_r_ns      !< non-structural root nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
  real,dimension(ilg,icc) :: nlitrveg_r_s       !< structural root nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$

  temp_leaf = 0.0
  temp_stem = 0.0
  temp_root = 0.0
  gleafltr_fr = 0.0
  stemltr_fr  = 0.0
  rootltr_fr  = 0.0
  !
  do i = il1,il2
     do j = 1,icc
        nlitrveg_l(i,j) = 0.0
        nlitrveg_s(i,j) = 0.0
        nlitrveg_r(i,j) = 0.0
        nlitrveg_l_ns(i,j) = 0.0
        nlitrveg_l_s(i,j) = 0.0
        nlitrveg_s_ns(i,j) = 0.0
        nlitrveg_s_s(i,j) = 0.0
        nlitrveg_r_ns(i,j) = 0.0
        nlitrveg_r_s(i,j) = 0.0
        ngleafmas_ns_in(i,j) = 0.0
        ngleafmas_s_in(i,j) = 0.0
        ngleafmas_in(i,j) = 0.0
        nbleafmas_in(i,j) = 0.0
        nstemmass_ns_in(i,j) = 0.0
        nstemmass_s_in(i,j) = 0.0
        nstemmass_in(i,j) = 0.0
        nrootmass_ns_in(i,j) = 0.0
        nrootmass_s_in(i,j) = 0.0
        nrootmass_in(i,j) = 0.0
        nlitrmass_in(i,j) = 0.0
        gl2bl_grass_nflux(i,j) = 0.0
        gl2bl_grass_nflux_ns(i,j) = 0.0
        gl2bl_grass_nflux_s(i,j) = 0.0
     end do ! loop 110
        nlitrmass_in(i,iccp1) = 0.0
  end do ! loop 100

  !> Litter fall nitrogen is calculated based on litter fall carbon multiplied by the C:N ratio:
  !! (see equation 11 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ LF_i = \frac{(1-r_i)LF_{i,C}}{C:N_i} \f$
  do j = 1,icc
     pftkind = ctempfts(j)
     do i = il1,il2

       if (fcancmx(i,j) > 0.0) then

          !! Leaf litter fall nitrogen:    
          select case (pftkind)
          case ('GrassC3 ','GrassC4 ','Sedge   ')

             if (ngleafmas(i,j) > 0.0 .and. gleafmas(i,j) > 0.0) gl2bl_grass_nflux(i,j) = &
                       (gl2bl_grass_cflux(i,j) / (gl2bl_grass_cflux(i,j) + gleafmas(i,j))) &
                       * ngleafmas(i,j)
             if (bleafmas(i,j) > 0.0)  nlitrveg_l(i,j) = &
                                  (bl2ltr_grass_cflux(i,j) / bleafmas(i,j)) * nbleafmas(i,j)           

          case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr','BdlDDrTr',&
                'CropC3  ','CropC4  ','BdlEvgSh','BdlDCoSh') ! trees and crops   
   
             if (gleafmas(i,j) > 0.0) then
                temp_leaf = (deltat/963.62) * tltrleaf(i,j)
                gleafltr_fr = min(temp_leaf / (temp_leaf + gleafmas(i,j)), 1.0)
             end if

             !! resorption of N for deciduous PFTs during leaf fall
             if (lfstatus(i,j) == 3  .and. (pftkind == 'BdlDCoTr' .or. pftkind == 'BdlDDrTr' &
                                          .or. pftkind == 'NdlDcTr ' .or. pftkind == 'BdlDCoSh')) then
                nlitrveg_l(i,j) = (1.0 - nresorp_coeff(sort(j))) * gleafltr_fr &
                                  * ngleafmas(i,j)

             else
                nlitrveg_l(i,j) = gleafltr_fr * ngleafmas(i,j)
             end if
             temp_leaf = 0.0
             gleafltr_fr = 0.0

          case default
             print * ,'Unknown CTEM PFT in nlitter ',pftkind
             call errorHandler('nlitter', - 1)
          end select
   
          !! Stem litter fall nitrogen: 
          if (stemmass(i,j) > 0.0) then
             temp_stem = (deltat/963.62) * tltrstem(i,j)
             stemltr_fr = min(temp_stem / (temp_stem + stemmass(i,j)), 1.0)
          end if
          nlitrveg_s(i,j) = stemltr_fr * nstemmass(i,j)
          stemltr_fr = 0.0
          temp_stem = 0.0

          !! Root litter fall nitrogen: 
          if (rootmass(i,j).ne.0.0) then
             temp_root = (deltat/963.62) * tltrroot(i,j)
             rootltr_fr = min(temp_root / (temp_root + rootmass(i,j)), 1.0)
          end if
          nlitrveg_r(i,j) = rootltr_fr * nrootmass(i,j)
          temp_root = 0.0
          rootltr_fr = 0.0

      end if !for fcancmx
     end do ! loop 130
  end do ! loop 120

  !> Updating pools and checking conservation 
  do j= 1,icc
     pftkind = ctempfts(j)
     do i = il1,il2
      
        ngleafmas_ns_in(i,j) = ngleafmas_ns(i,j)
        ngleafmas_s_in(i,j) = ngleafmas_s(i,j)
        ngleafmas_in(i,j) = ngleafmas(i,j)
        nbleafmas_in(i,j) = nbleafmas(i,j)
        nstemmass_ns_in(i,j) = nstemmass_ns(i,j)
        nstemmass_s_in(i,j) = nstemmass_s(i,j)
        nstemmass_in(i,j) = nstemmass(i,j)
        nrootmass_ns_in(i,j) = nrootmass_ns(i,j)
        nrootmass_s_in(i,j) = nrootmass_s(i,j)
        nrootmass_in(i,j) = nrootmass(i,j)
        nlitrmass_in(i,j) = nlitrmass(i,j)

        !! 1) ngleafmas

        if (ngleafmas_in(i,j) > 0.0) then

           select case (pftkind)

           case ('GrassC3 ','GrassC4 ','Sedge   ')
              gl2bl_grass_nflux_ns(i,j) = gl2bl_grass_nflux(i,j) * (ngleafmas_ns_in(i,j) / ngleafmas_in(i,j))
              gl2bl_grass_nflux_s(i,j) = gl2bl_grass_nflux(i,j) * (ngleafmas_s_in(i,j) / ngleafmas_in(i,j))
              ngleafmas_ns(i,j) = ngleafmas_ns_in(i,j) - gl2bl_grass_nflux_ns(i,j) 
              ngleafmas_s(i,j) = ngleafmas_s_in(i,j) - gl2bl_grass_nflux_s(i,j) 

           case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr','BdlDDrTr',&
                 'CropC3  ','CropC4  ','BdlEvgSh','BdlDCoSh')
              nlitrveg_l_ns(i,j) = nlitrveg_l(i,j) * (ngleafmas_ns_in(i,j) / ngleafmas_in(i,j))
              nlitrveg_l_s(i,j) = nlitrveg_l(i,j) * (ngleafmas_s_in(i,j) / ngleafmas_in(i,j))
              ngleafmas_ns(i,j) = ngleafmas_ns_in(i,j) - nlitrveg_l_ns(i,j)
              ngleafmas_s(i,j) = ngleafmas_s_in(i,j) - nlitrveg_l_s(i,j) 
            
           case default
              print * ,'Unknown CTEM PFT in nlitter ',pftkind
              call errorHandler('nlitter', - 2)
           end select

        end if !for ngleafmas

        ngleafmas(i,j) = ngleafmas_ns(i,j) + ngleafmas_s(i,j)

        !! Adjusting to prevent negative ngleafmas_ns(i,j)
        if (ngleafmas_ns(i,j) < 0.0) then
           write(*,*)'ngleafmas_ns lt zero at i=',i,' for pft=',j
           ! write(*,*)'ngleafmas_ns(i,j)     = ',ngleafmas_ns(i,j)
           ! write(*,*)'nlitrveg_l(i,j)       = ',nlitrveg_l(i,j)
           ! write(*,*)'gl2bl_grass_nflux(i,j)= ',gl2bl_grass_nflux(i,j)
           write(*,*)'adjusting to prevent negative ngleafmas_ns(i,j)'

           select case (pftkind)

           case ('GrassC3 ','GrassC4 ','Sedge   ')
              gl2bl_grass_nflux(i,j) = gl2bl_grass_nflux(i,j) - gl2bl_grass_nflux(i,j) &
                                       * abs(ngleafmas_ns(i,j) / gl2bl_grass_nflux_ns(i,j))
              gl2bl_grass_nflux_ns(i,j) = gl2bl_grass_nflux(i,j) * (ngleafmas_ns_in(i,j) / ngleafmas_in(i,j))
              gl2bl_grass_nflux_s(i,j) = gl2bl_grass_nflux(i,j) * (ngleafmas_s_in(i,j) / ngleafmas_in(i,j))
              ngleafmas_s(i,j) = ngleafmas_s_in(i,j) - gl2bl_grass_nflux_s(i,j) 
           case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr','BdlDDrTr',&
                 'CropC3  ','CropC4  ','BdlEvgSh','BdlDCoSh')
              nlitrveg_l(i,j) = nlitrveg_l(i,j) - nlitrveg_l(i,j) &
                                * abs(ngleafmas_ns(i,j) / nlitrveg_l_ns(i,j))
              nlitrveg_l_ns(i,j) = nlitrveg_l(i,j) * (ngleafmas_ns_in(i,j) / ngleafmas_in(i,j))
              nlitrveg_l_s(i,j) = nlitrveg_l(i,j) * (ngleafmas_s_in(i,j) / ngleafmas_in(i,j))
              ngleafmas_s(i,j) = ngleafmas_s_in(i,j) - nlitrveg_l_s(i,j)
           case default
              print * ,'Unknown CTEM PFT in nlitter ',pftkind
              call errorHandler('nlitter', - 3)
           end select
         
           ngleafmas_ns(i,j) = 0.0
           ngleafmas(i,j) = ngleafmas_ns(i,j) + ngleafmas_s(i,j)
         
        end if ! for ngleafmas_ns

        !! Adjusting to prevent negative ngleafmas_s(i,j)
        if (ngleafmas_s(i,j) < 0.0) then
           write(*,*)'ngleafmas_s lt zero at i=',i,' for pft=',j
           ! write(*,*)'ngleafmas_s(i,j)      = ',ngleafmas_s(i,j)
           ! write(*,*)'nlitrveg_l(i,j)       = ',nlitrveg_l(i,j)
           ! write(*,*)'gl2bl_grass_nflux(i,j)= ',gl2bl_grass_nflux(i,j)
           write(*,*)'adjusting to prevent negative ngleafmas_s(i,j)'
         
           select case (pftkind)

           case ('GrassC3 ','GrassC4 ','Sedge   ')
              gl2bl_grass_nflux(i,j) = gl2bl_grass_nflux(i,j) - gl2bl_grass_nflux(i,j) &
                                       * abs(ngleafmas_s(i,j) / gl2bl_grass_nflux_s(i,j))
              gl2bl_grass_nflux_ns(i,j) = gl2bl_grass_nflux(i,j) * (ngleafmas_ns_in(i,j) / ngleafmas_in(i,j))
              gl2bl_grass_nflux_s(i,j) = gl2bl_grass_nflux(i,j) * (ngleafmas_s_in(i,j) / ngleafmas_in(i,j))
              ngleafmas_ns(i,j) = ngleafmas_ns_in(i,j) - gl2bl_grass_nflux_ns(i,j) 
           case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr','BdlDDrTr',&
                 'CropC3  ','CropC4  ','BdlEvgSh','BdlDCoSh')
              nlitrveg_l(i,j) = nlitrveg_l(i,j) - nlitrveg_l(i,j) &
                                * abs(ngleafmas_s(i,j) / nlitrveg_l_s(i,j))
              nlitrveg_l_ns(i,j) = nlitrveg_l(i,j) * (ngleafmas_ns_in(i,j) / ngleafmas_in(i,j))
              nlitrveg_l_s(i,j) = nlitrveg_l(i,j) * (ngleafmas_s_in(i,j) / ngleafmas_in(i,j))
              ngleafmas_ns(i,j) = ngleafmas_ns_in(i,j) - nlitrveg_l_ns(i,j)
           case default
              print * ,'Unknown CTEM PFT in nlitter ',pftkind
              call errorHandler('nlitter', - 4)
           end select
         
           ngleafmas_s(i,j) = 0.0
           ngleafmas(i,j) = ngleafmas_ns(i,j) + ngleafmas_s(i,j)
         
        end if ! for ngleafmas_s

        !! Checking conservation for ngleafmas(i,j) 
        select case (pftkind)

        case ('GrassC3 ','GrassC4 ','Sedge   ')
           if ((ngleafmas_in(i,j) - ngleafmas(i,j) - gl2bl_grass_nflux(i,j)) > &
               ntolerance) then
               write(*,*)'imbalance in ngleafmas at (i)=',i,'pft=',j
               call errorHandler('nlitter', - 5)
           end if
           if ((ngleafmas_ns_in(i,j) - ngleafmas_ns(i,j) - gl2bl_grass_nflux_ns(i,j)) > &
               ntolerance) then
               write(*,*)'imbalance in ngleafmas_ns at (i)=',i,'pft=',j
               call errorHandler('nlitter', - 5)
           end if
           if ((ngleafmas_s_in(i,j) - ngleafmas_s(i,j) - gl2bl_grass_nflux_s(i,j)) > &
               ntolerance) then
               write(*,*)'imbalance in ngleafmas_s at (i)=',i,'pft=',j
               call errorHandler('nlitter', - 5)
           end if

        case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr','BdlDDrTr',&
              'CropC3  ','CropC4  ','BdlEvgSh','BdlDCoSh')
           if ((ngleafmas_in(i,j) - ngleafmas(i,j) - nlitrveg_l(i,j)) > ntolerance &
                                                         ) then
              write(*,*)'imbalance in ngleafmas at (i)=',i,'pft=',j
              call errorHandler('nlitter', - 6)
           end if
           if ((ngleafmas_ns_in(i,j) - ngleafmas_ns(i,j) - nlitrveg_l_ns(i,j)) > ntolerance &
                                                         ) then
              write(*,*)'imbalance in ngleafmas_ns at (i)=',i,'pft=',j
              call errorHandler('nlitter', - 6)
           end if
           if ((ngleafmas_s_in(i,j) - ngleafmas_s(i,j) - nlitrveg_l_s(i,j)) > ntolerance &
                                                         ) then
              write(*,*)'imbalance in ngleafmas_s at (i)=',i,'pft=',j
              call errorHandler('nlitter', - 6)
           end if

        case default
           print * ,'Unknown CTEM PFT in nlitter ',pftkind
           call errorHandler('nlitter', - 7)
        end select

        !! 2) nbleafmas

        select case (pftkind)

        case ('GrassC3 ','GrassC4 ','Sedge   ')
           nbleafmas(i,j) = nbleafmas(i,j) + gl2bl_grass_nflux(i,j) - nlitrveg_l(i,j)
         
        case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr','BdlDDrTr',&
              'CropC3  ','CropC4  ','BdlEvgSh','BdlDCoSh')
           nbleafmas(i,j) = 0.0
         
        case default
           print * ,'Unknown CTEM PFT in nlitter ',pftkind
           call errorHandler('nlitter', - 8)
        end select

        !! Adjusting to prevent negative nbleafmas(i,j)
        if (nbleafmas(i,j) < 0.0) then
           write(*,*)'nbleafmas lt zero at i=',i,' for pft=',j
           ! write(*,*)'nbleafmas(i,j)        = ',nbleafmas(i,j)
           ! write(*,*)'nlitrveg_l(i,j)       = ',nlitrveg_l(i,j)
           ! write(*,*)'gl2bl_grass_nflux(i,j)= ',gl2bl_grass_nflux(i,j)
           write(*,*)'adjusting to prevent negative nbleafmas(i,j)'

           select case (pftkind)

           case ('GrassC3 ','GrassC4 ','Sedge   ')
              nlitrveg_l(i,j) = nlitrveg_l(i,j) - nlitrveg_l(i,j) &
                                * abs(nbleafmas(i,j) / nlitrveg_l(i,j))
            
           case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr','BdlDDrTr',&
                 'CropC3  ','CropC4  ','BdlEvgSh','BdlDCoSh')

              ! do nothing.
            
           case default
              print * ,'Unknown CTEM PFT in nlitter ',pftkind
              call errorHandler('nlitter', - 9)
           end select
         
           nbleafmas(i,j) = 0.0
        end if

        !! Checking conservation for nbleafmas(i,j)
        if ((nbleafmas_in(i,j) - nbleafmas(i,j) + gl2bl_grass_nflux(i,j) &
             - nlitrveg_l(i,j)) > ntolerance) then
           write(*,*)'imbalance in nbleafmas at (i)=',i,'pft=',j
           call errorHandler('nlitter', - 10)
        end if

        !! 3) nstemmass

        if (nstemmass(i,j) > 0.0) then
           nlitrveg_s_ns(i,j) = nlitrveg_s(i,j) * (nstemmass_ns_in(i,j) / nstemmass_in(i,j))
           nlitrveg_s_s(i,j) = nlitrveg_s(i,j) * (nstemmass_s_in(i,j) / nstemmass_in(i,j)) 
           nstemmass_ns(i,j) = nstemmass_ns_in(i,j) - nlitrveg_s_ns(i,j)
           nstemmass_s(i,j) = nstemmass_s_in(i,j) - nlitrveg_s_s(i,j)
        end if

        nstemmass(i,j) = nstemmass_ns(i,j) + nstemmass_s(i,j)

        !! Adjusting to prevent negative nstemmass_ns(i,j)
        if (nstemmass_ns(i,j) < 0.0) then
           write(*,*)'nstemmass_ns lt zero at i=',i,' for pft=',j
           ! write(*,*)'nstemmass_ns(i,j)= ',nstemmass_ns(i,j)
           ! write(*,*)'nlitrveg_s(i,j)  = ',nlitrveg_s(i,j)
           write(*,*)'adjusting to prevent negative nstemmass_ns(i,j)'

           nlitrveg_s(i,j) = nlitrveg_s(i,j) - nlitrveg_s(i,j) &
                             * abs(nstemmass_ns(i,j) / nlitrveg_s_ns(i,j))
           nlitrveg_s_ns(i,j) = nlitrveg_s(i,j) * (nstemmass_ns_in(i,j) / nstemmass_in(i,j))
           nlitrveg_s_s(i,j) = nlitrveg_s(i,j) * (nstemmass_s_in(i,j) / nstemmass_in(i,j))
           nstemmass_ns(i,j) = 0.0
           nstemmass_s(i,j) = nstemmass_s_in(i,j) - nlitrveg_s_s(i,j) 
           nstemmass(i,j) = nstemmass_ns(i,j) + nstemmass_s(i,j)
        end if

        !! Adjusting to prevent negative nstemmass_s(i,j)
        if (nstemmass_s(i,j) < 0.0) then
           write(*,*)'nstemmass_s lt zero at i=',i,' for pft=',j
           !write(*,*)'nstemmass_s(i,j) = ',nstemmass_s(i,j)
           !write(*,*)'nlitrveg_s(i,j)  = ',nlitrveg_s(i,j)
           write(*,*)'adjusting to prevent negative nstemmass_s(i,j)'

           nlitrveg_s(i,j) = nlitrveg_s(i,j) - nlitrveg_s(i,j) &
                             * abs(nstemmass_s(i,j) / nlitrveg_s_s(i,j))
           nlitrveg_s_ns(i,j) = nlitrveg_s(i,j) * (nstemmass_ns_in(i,j) / nstemmass_in(i,j))
           nlitrveg_s_s(i,j) = nlitrveg_s(i,j) * (nstemmass_s_in(i,j) / nstemmass_in(i,j))
           nstemmass_s(i,j) = 0.0
           nstemmass_ns(i,j) = nstemmass_ns_in(i,j) - nlitrveg_s_ns(i,j)
           nstemmass(i,j) = nstemmass_ns(i,j) + nstemmass_s(i,j)
        end if

        !! Checking conservation for nstemmass(i,j)
        if ((nstemmass_in(i,j) - nstemmass(i,j) - nlitrveg_s(i,j)) > ntolerance &
           ) then
           write(*,*)'imbalance in nstemmass at (i)=',i,'pft=',j
           call errorHandler('nlitter', - 11)
        end if

        !! 4) nrootmass

        if (nrootmass(i,j) > 0.0) then
           nlitrveg_r_ns(i,j) = nlitrveg_r(i,j) * (nrootmass_ns_in(i,j) / nrootmass_in(i,j))
           nlitrveg_r_s(i,j) = nlitrveg_r(i,j) * (nrootmass_s_in(i,j) / nrootmass_in(i,j)) 
           nrootmass_ns(i,j) = nrootmass_ns(i,j) - nlitrveg_r_ns(i,j)
           nrootmass_s(i,j) = nrootmass_s(i,j) - nlitrveg_r_s(i,j) 
        end if

        nrootmass(i,j) = nrootmass_ns(i,j) + nrootmass_s(i,j)

        !! Adjusting to prevent negative nrootmass_ns(i,j)
        if (nrootmass_ns(i,j) < 0.0) then
           write(*,*)'nrootmass_ns lt zero at i=',i,' for pft=',j
           ! write(*,*)'nrootmass_ns(i,j)= ',nrootmass_ns(i,j)
           ! write(*,*)'nlitrveg_r(i,j)  = ',nlitrveg_r(i,j)
           write(*,*)'adjusting to prevent negative nrootmass_ns(i,j)'

           nlitrveg_r(i,j) = nlitrveg_r(i,j) - nlitrveg_r(i,j) &
                             * abs(nrootmass_ns(i,j)/nlitrveg_r_ns(i,j))
           nlitrveg_r_ns(i,j) = nlitrveg_r(i,j) * (nrootmass_ns_in(i,j) / nrootmass_in(i,j))
           nlitrveg_r_s(i,j) = nlitrveg_r(i,j) * (nrootmass_s_in(i,j) / nrootmass_in(i,j))
           nrootmass_ns(i,j) = 0.0
           nrootmass_s(i,j) = nrootmass_s(i,j) - nlitrveg_r_s(i,j) 
           nrootmass(i,j) = nrootmass_ns(i,j) + nrootmass_s(i,j)
        end if

        !! Adjusting to prevent negative nrootmass_s(i,j)
        if (nrootmass_s(i,j) < 0.0) then
           write(*,*)'nrootmass_s lt zero at i=',i,' for pft=',j
           ! write(*,*)'nrootmass_s(i,j) = ',nrootmass_s(i,j)
           ! write(*,*)'nlitrveg_r(i,j)  = ',nlitrveg_r(i,j)
           write(*,*)'adjusting to prevent negative nrootmass_s(i,j)'

           nlitrveg_r(i,j) = nlitrveg_r(i,j) - nlitrveg_r(i,j) &
                             * abs(nrootmass_s(i,j)/nlitrveg_r_s(i,j))
           nlitrveg_r_ns(i,j) = nlitrveg_r(i,j) * (nrootmass_ns_in(i,j) / nrootmass_in(i,j))
           nlitrveg_r_s(i,j) = nlitrveg_r(i,j) * (nrootmass_s_in(i,j) / nrootmass_in(i,j))
           nrootmass_s(i,j) = 0.0
           nrootmass_ns(i,j) = nrootmass_ns(i,j) - nlitrveg_r_ns(i,j) 
           nrootmass(i,j) = nrootmass_ns(i,j) + nrootmass_s(i,j)
        end if

        !! Checking conservation for nrootmass(i,j)
        if ((nrootmass_in(i,j) - nrootmass(i,j) - nlitrveg_r(i,j)) > ntolerance &
           ) then
           write(*,*)'imbalance in nrootmass at (i)=',i,'pft=',j
           call errorHandler('nlitter', - 12)
        end if

        !! 5) nlitrmass

        nlitrmass(i,j) = nlitrmass(i,j) + nlitrveg_l(i,j) + nlitrveg_s(i,j) + nlitrveg_r(i,j)

        if (nlitrmass(i,j) < 0.0) then
           write(*,*)'nlitrmass lt zero at i=',i,' for pft=',j
           write(*,*)'nlitrmass(i,j)  = ',nlitrmass(i,j)
           write(*,*)'nlitrveg_l(i,j) = ',nlitrveg_l(i,j)
           write(*,*)'nlitrveg_s(i,j) = ',nlitrveg_s(i,j)
           write(*,*)'nlitrveg_r(i,j) = ',nlitrveg_r(i,j)
           call errorHandler('nlitter', - 13)
        end if

        !! Checking conservation for nlitrmass(i,j)
        if ((nlitrmass_in(i,j) - nlitrmass(i,j) + nlitrveg_l(i,j)    &
                                                + nlitrveg_s(i,j)    &
                                                + nlitrveg_r(i,j) )  &
                                                > ntolerance         &
                                               ) then
           write(*,*)'imbalance in nlitrmass at (i)=',i,'pft=',j
           write(*,*)'nlitrmass_in(i,j) = ',nlitrmass_in(i,j)
           write(*,*)'nlitrmass(i,j)    = ',nlitrmass(i,j)
           write(*,*)'nlitrveg_l(i,j)   = ',nlitrveg_l(i,j)
           write(*,*)'nlitrveg_s(i,j)   = ',nlitrveg_s(i,j)
           write(*,*)'nlitrveg_r(i,j)   = ',nlitrveg_r(i,j)
           call errorHandler('nlitter', - 14)
        end if

     end do ! loop 150
  end do ! loop 140

  return
end subroutine nlitter
