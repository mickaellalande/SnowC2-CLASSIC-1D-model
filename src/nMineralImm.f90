!> \file
!! Calculates N mineralization and immobilization
!> @author A. Asaadi
!!

subroutine nmineralimm(il1, il2, spinfast, soilcmas_bulk, c2nveg_litr, c2nveg_humus, & ! In
                       ltrestep, screstep, litrmass_bulk, co2conc, & ! In
                       nlitrmass, soilnmas, nh4_mass, no3_mass, & ! In/Out
                       nmineralveg_litr, nmineralveg_humus, nimmobilveg_nh4, & ! Out
                       nimmobilveg_no3) ! Out

  use classicParams, only : ilg, ignd, icc, iccp1, iccp2, ctempfts, ntolerance, &
                            const_c2n_humus, initial_c2n

  implicit none

  integer,intent(in) :: il1       !< il1=1
  integer,intent(in) :: il2       !< il2=ilg (no. of grid cells in latitude circle)
  integer,intent(in) :: spinfast  !< spinup factor for soil carbon & nitrogen pools
  real,intent(in),dimension(ilg,iccp1) :: ltrestep       !< litter carbon respiration rate for individual PFTs + bare soil, \f$kg C m^{-2} day^{-1}\f$ 
  real,intent(in),dimension(ilg,iccp1) :: screstep       !< soil carbon respiration rate for individual PFTs + bare soil, \f$kg C m^{-2} day^{-1}\f$
  real,intent(in),dimension(ilg,iccp1) :: c2nveg_litr    !< simulated C:N ratio for litter, \f$g C g N^{-1}\f$
  real,intent(in),dimension(ilg,iccp1) :: c2nveg_humus   !< simulated C:N ratio for soil, \f$g C g N^{-1}\f$
  real,intent(in),dimension(ilg,iccp2) :: soilcmas_bulk  !< soil carbon for each PFT + bare soil + LUC product, \f$kg C m^{-2}\f$
  real,intent(in),dimension(ilg,iccp2) :: litrmass_bulk  !< litter carbon for each PFT + bare soil + LUC product, \f$kg C m^{-2}\f$
  real,intent(in),dimension(ilg) :: co2conc              !< atmospheric CO2 concentration, ppm

  real,intent(out),dimension(ilg,iccp1) :: nmineralveg_litr  !< nitrogen mineralization from litter, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,iccp1) :: nmineralveg_humus !< nitrogen mineralization from soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,iccp1) :: nimmobilveg_nh4   !< nitrogen immobilization from NH4+ to soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,iccp1) :: nimmobilveg_no3   !< nitrogen immobilization from NO3- to soil, \f$g N m^{-2} day^{-1}\f$

  real,intent(inout),dimension(ilg,iccp2) :: nlitrmass       !< litter nitrogen for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,iccp2) :: soilnmas        !< soil nitrogen for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,iccp1) :: nh4_mass        !< NH4+ for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,iccp1) :: no3_mass        !< NO3- for individual PFTs + bare soil, \f$g N m^{-2}\f$

  ! local variables
  integer :: i, j
  character(8) :: pftkind
  real :: frac_of_nh4, frac_of_no3               !< fractions of NH4+ & NO3- in total inorganic nitrogen 
  real :: temp,scale_factor
  real,dimension(ilg,iccp1)   :: nlitrmass_in    !< previous litter nitrogen for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,dimension(ilg,iccp1)   :: soilnmas_in     !< previous soil nitrogen for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,dimension(ilg,iccp1)   :: nh4_mass_in     !< previous NH4+ for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,dimension(ilg,iccp1)   :: no3_mass_in     !< previous NO3- for individual PFTs + bare soil, \f$g N m^{-2}\f$

  !
  temp = 0.0
  scale_factor = 0.0
  frac_of_nh4 = 0.0
  frac_of_no3 = 0.0

  do i = il1,il2
     do j = 1,iccp1
        nmineralveg_litr(i,j) = 0.0
        nmineralveg_humus(i,j) = 0.0
        nimmobilveg_nh4(i,j) = 0.0
        nimmobilveg_no3(i,j) = 0.0
        nlitrmass_in(i,j) = 0.0
        soilnmas_in(i,j) = 0.0
        nh4_mass_in(i,j) = 0.0
        no3_mass_in(i,j) = 0.0
    end do ! loop 110
  end do ! loop 100

  !> N mineralization from litter and soil N pools to the NH4 pool:
  !! (see equation 16 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ M_{D,NH_4} = \frac{R_{h,D}}{C:N_D} \f$
  !! \f$ M_{H,NH_4} = \frac{R_{h,H}}{C:N_H} \f$
  do i = il1,il2
     do j = 1,iccp1

        if (nlitrmass(i,j) > 0.0 .and. litrmass_bulk(i,j) > 0.0) nmineralveg_litr(i,j) = ltrestep(i,j) / &
                                   (litrmass_bulk(i,j) / nlitrmass(i,j))

        if (soilnmas(i,j) > 0.0 .and. soilcmas_bulk(i,j) > 0.0) nmineralveg_humus(i,j) = screstep(i,j) / &
                                  (soilcmas_bulk(i,j) / soilnmas(i,j))

        nmineralveg_litr(i,j) = min(nmineralveg_litr(i,j), nlitrmass(i,j))
        nmineralveg_humus(i,j) = min(nmineralveg_humus(i,j), soilnmas(i,j)/real(spinfast))

        if (nmineralveg_litr(i,j) < 0.0)then
            write(*,*)'nmineralveg_litr lt zero at i=',i,' for pft=',j
            write(*,*)'nmineralveg_litr(i,j) = ',nmineralveg_litr(i,j)
            call errorHandler('nmineralization', - 1)
        end if
        if (nmineralveg_humus(i,j) < 0.0)then
            write(*,*)'nmineralveg_humus lt zero at i=',i,' for pft=',j
            write(*,*)'nmineralveg_humus(i,j) = ',nmineralveg_humus(i,j)
            call errorHandler('nmineralization', - 2)
        end if

     end do ! loop 130
  end do ! loop 120

  !> N immobilization from NH4 and NO3 pools to soil organic N pool:
  !! (see equation 17 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ O_{NH_4} = max(0,(\frac{C_H}{C:N_H}-N_H)\frac{N_{NH_4}}{N_{NH_4}+N_{NO_3}}) k_O \f$
  !! \f$ O_{NO_3} = max(0,(\frac{C_H}{C:N_H}-N_H)\frac{N_{NO_3}}{N_{NH_4}+N_{NO_3}}) k_O \f$
  do i = il1,il2
     do j = 1,iccp1

        if (nh4_mass(i,j)+no3_mass(i,j) > 0.0) then
            frac_of_nh4 = nh4_mass(i,j) / (nh4_mass(i,j) + no3_mass(i,j))
            frac_of_no3 = no3_mass(i,j) / (nh4_mass(i,j) + no3_mass(i,j))

            ! We let the C:N ratio of the humus pool to increase slightly above the constant value of 13.
            ! The increase in humus pool's C:N ratio under elevated co2 cocentration has been mentioned
            ! in other studies as well (e.g., Luo et al., 2006).
            ! As we have learnt from our co2-only runs, this adjustment will prevent too much reduction
            ! seen in the NH4+ pool which leads to reduction in N uptake by plant despite the increase
            ! seen in the plant N demand under elevated co2.

            if (nlitrmass(i,j) > 0.0) then
               scale_factor= 1.0 !(c2nveg_litr(i,j)/52.)**0.05
            else
               scale_factor= 1.0
            end if

            nimmobilveg_nh4(i,j) = ((1000. * soilcmas_bulk(i,j) / (const_c2n_humus * scale_factor)) &
                                   - soilnmas(i,j)) * frac_of_nh4
            nimmobilveg_no3(i,j) = ((1000. * soilcmas_bulk(i,j) / (const_c2n_humus * scale_factor)) &
                                   - soilnmas(i,j)) * frac_of_no3

            nimmobilveg_nh4(i,j) = max (nimmobilveg_nh4(i,j), 0.0)
            nimmobilveg_no3(i,j) = max (nimmobilveg_no3(i,j), 0.0)

            ! Prevent negative NH4+ and NO3-: 
            if (nimmobilveg_nh4(i,j).gt.(nmineralveg_litr(i,j)+nmineralveg_humus(i,j)) ) then
              nimmobilveg_nh4(i,j) = min(nimmobilveg_nh4(i,j),nh4_mass(i,j)/real(spinfast))
            end if
            nimmobilveg_no3(i,j) = min(nimmobilveg_no3(i,j),no3_mass(i,j)/real(spinfast))

        end if
     end do ! loop 150
  end do ! loop 140


  !> Updating NH4+, NO3-, litter nitrogen, and soil nitrogen and checking conservation 
  do j = 1,iccp1
     do i = il1,il2

       nh4_mass_in(i,j) = nh4_mass(i,j) 
       no3_mass_in(i,j) = no3_mass(i,j)
       nlitrmass_in(i,j)= nlitrmass(i,j)
       soilnmas_in(i,j) = soilnmas(i,j)

       !! 1) nh4_mass
       nh4_mass(i,j)= nh4_mass(i,j) + nmineralveg_litr(i,j) &
                               + nmineralveg_humus(i,j) - nimmobilveg_nh4(i,j)

       !! Adjusting to prevent negative nh4_mass(i,j)
       if (nh4_mass(i,j) < 0.0) then
          write(*,*)'nh4_mass lt zero in nMineralImm.f90 at i=',i,' for pft=',j
          ! write(*,*)'nh4_mass(i,j)         = ',nh4_mass(i,j)
          ! write(*,*)'nmineralveg_litr(i,j) = ',nmineralveg_litr(i,j)
          ! write(*,*)'nmineralveg_humus(i,j)= ',nmineralveg_humus(i,j)
          ! write(*,*)'nimmobilveg_nh4(i,j)  = ',nimmobilveg_nh4(i,j)
          write(*,*)'adjusting to prevent negative nh4_mass(i,j)'
          temp = nimmobilveg_nh4(i,j)
          nimmobilveg_nh4(i,j) = nimmobilveg_nh4(i,j) - nimmobilveg_nh4(i,j) &
                                 * abs(nh4_mass(i,j) / temp)
          nh4_mass(i,j) = 0.0
          temp = 0.0
       end if

      !! Checking conservation for nh4_mass(i,j)
       if ((nh4_mass_in(i,j) - nh4_mass(i,j) + nmineralveg_litr(i,j)   &
          + nmineralveg_humus(i,j) - nimmobilveg_nh4(i,j)) &
          > ntolerance) then
          write(*,*)'imbalance in nh4_mass pool at (i)=',i,'pft=',j
          call errorHandler('nmineralization', - 3)
       end if

       !! 2) no3_mass
       no3_mass(i,j) = no3_mass(i,j) - nimmobilveg_no3(i,j)

       !! Adjusting to prevent negative no3_mass(i,j)
       if (no3_mass(i,j) < 0.0) then
          write(*,*)'no3_mass lt zero in nMineralImm.f90 at i=',i,' for pft=',j
          ! write(*,*)'no3_mass        = ',no3_mass(i,j)
          ! write(*,*)'nimmobilveg_no3 = ',nimmobilveg_no3(i,j)
          write(*,*)'adjusting to prevent negative no3_mass(i,j)'
          temp = nimmobilveg_no3(i,j)
          nimmobilveg_no3(i,j) = nimmobilveg_no3(i,j) - nimmobilveg_no3(i,j) &
                                 * abs(no3_mass(i,j) / temp)
          no3_mass(i,j) = 0.0
          temp = 0.0
       end if

       !! Checking conservation for no3_mass(i,j)
       if ((no3_mass_in(i,j) - no3_mass(i,j) - nimmobilveg_no3(i,j)) &
                                > ntolerance) then
          write(*,*)'imbalance in no3_mass pool at (i)=',i,'pft=',j
          call errorHandler('nmineralization', - 4)
       end if


       !! 3) nlitrmass
       nlitrmass(i,j) = nlitrmass(i,j) - nmineralveg_litr(i,j)

      !! Adjusting to prevent negative nlitrmass(i,j)
      if (nlitrmass(i,j) .lt. 0.0) then
         write(*,*)'nlitrmass lt zero in nMineralImm.f90 at i=',i,' for pft=',j
         ! write(*,*)'nlitrmass(i,j)        = ',nlitrmass(i,j)
         ! write(*,*)'nmineralveg_litr(i,j) = ',nmineralveg_litr(i,j)
         write(*,*)'adjusting to prevent negative nlitrmass(i,j)'
         temp = nmineralveg_litr(i,j)
         nmineralveg_litr(i,j) = nmineralveg_litr(i,j) - nmineralveg_litr(i,j) &
                                 * abs(nlitrmass(i,j) / temp)
         nlitrmass(i,j) = 0.0
         temp = 0.0
      end if

      !! Checking conservation for nlitrmass(i,j)
      if ((nlitrmass_in(i,j) - nlitrmass(i,j) - nmineralveg_litr(i,j)) &
                                 > ntolerance) then
         write(*,*)'imbalance in nlitrmass pool at (i)=',i,'pft=',j
         call errorHandler('nmineralization', - 5)
      end if

      !! 4) soilnmas
      soilnmas(i,j) = soilnmas(i,j) + real(spinfast) * (- nmineralveg_humus(i,j) &
                                    + nimmobilveg_nh4(i,j) + nimmobilveg_no3(i,j))

      if (spinfast > 1 .and. soilnmas(i,j) < 0.0 .and. &
         abs(soilnmas(i,j)) < 1.E-15) soilnmas(i,j) = 0.0

      !! Adjusting to prevent negative soilnmas(i,j)
      if (soilnmas(i,j) < 0.0) then
         write(*,*)'soilnmas lt zero in nMineralImm.f90 at i=',i,' for pft=',j
         ! write(*,*)'soilnmas(i,j)         = ',soilnmas(i,j)
         ! write(*,*)'nmineralveg_humus(i,j)= ',nmineralveg_humus(i,j)
         ! write(*,*)'nimmobilveg_nh4(i,j)  = ',nimmobilveg_nh4(i,j)
         ! write(*,*)'nimmobilveg_no3(i,j)  = ',nimmobilveg_no3(i,j)
         write(*,*)'adjusting to prevent negative soilnmas(i,j)'
         temp = real(spinfast) * nmineralveg_humus(i,j)
         nmineralveg_humus(i,j) = nmineralveg_humus(i,j) - nmineralveg_humus(i,j) &
                                  * abs(soilnmas(i,j) / temp)
         soilnmas(i,j) = 0.0
         temp = 0.0
      end if

      !! Checking conservation for soilnmas(i,j)
      if ((soilnmas_in(i,j) - soilnmas(i,j) - (nmineralveg_humus(i,j) &
          + nimmobilveg_nh4(i,j) + nimmobilveg_no3(i,j))) &
          > ntolerance .and. spinfast == 1) then
         write(*,*)'imbalance in soilnmas pool at (i)=',i,'pft=',j
         call errorHandler('nmineralization', - 6)
      end if

     end do ! loop 170
  end do ! loop 160

  return
end subroutine nmineralimm
