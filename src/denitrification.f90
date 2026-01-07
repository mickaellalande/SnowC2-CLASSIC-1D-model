!> \file
!! Calculate denitrification.
!> @author A. Asaadi


subroutine denitrific(il1, il2, zbotw, isand, tbar, thliq, & ! In
                      THFC, THLW, solctempavg, & ! In
                      no3_mass, & ! In/Out
                      no_denitveg, n2o_denitveg, n2_denitveg, denitveg) ! Out


  use classicParams, only : ilg, ignd, icc, iccp1, mu_denit, &
                            topt_denit, ntolerance, tanhq10, &
                            denitswthrsh, nitdenitdpth, &
                            mu_no_denit, mu_n2o_denit, ctempfts

  implicit none

  integer,intent(in) :: il1        !< il1=1
  integer,intent(in) :: il2        !< il2=ilg (no. of grid cells in latitude circle)
  integer,intent(in),dimension(ilg,ignd) :: isand        !< flag for bedrock or ice in a soil layer
  real,intent(in),dimension(ilg,ignd) :: tbar            !< soil temperature, K 
  real,intent(in),dimension(ilg,ignd) :: thliq           !< volumetric soil moisture, \f$m^3 m^{-3}\f$ 
  real,intent(in),dimension(ilg,ignd) :: zbotw           !< bottom of each soil layer, \f$m\f$ 
  real,intent(in),dimension(ilg,ignd) :: THLW            !< volumetric soil moisture at wilting point, \f$m^3 m^{-3}\f$
  real,intent(in),dimension(ilg,ignd) :: THFC            !< volumetric soil moisture at field capacity, \f$m^3 m^{-3}\f$
  real,intent(in),dimension(ilg) :: solctempavg          !< average soil temperature over soil layers in which denitrification occurs, K
  real,intent(inout),dimension(ilg,iccp1) :: no3_mass    !< NO3- for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,intent(out),dimension(ilg,iccp1) :: n2_denitveg   !< denitrification N2 loss for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$ 
  real,intent(out),dimension(ilg,iccp1) :: denitveg      !< denitrification flux from NO3- for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,iccp1) :: no_denitveg   !< denitrification NO loss for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,iccp1) :: n2o_denitveg  !< denitrification N2O loss for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$

  ! local variables
  integer :: i, j, k
  integer,dimension(ilg) :: ntsl            !< number of soil layers in which denitrification occurs 
  real,dimension(ilg,iccp1) :: soilcq10     !< Q10 for denitrification
  real,dimension(ilg,iccp1) :: tempq10s     !< average soil temperature over soil layers in which denitrification occurs, C 
  real,dimension(ilg,iccp1) :: fdenit_t     !< soil temperature dependence over soil layers in which denitrification occurs 
  real,dimension(ilg) :: fdenit_sw          !< volumetric soil moisture dependence over soil layers in which denitrification occurs
  real,dimension(ilg,ignd) :: soilw         !< soil wetness 
  real,dimension(ilg) :: soilwavg           !< average soil wetness over soil layers in which denitrification occurs 
  real,dimension(ilg,iccp1) :: no3_mass_in  !< previous NO3- for individual PFTs + bare soil, \f$g N m^{-2}\f$

  real :: temp

  temp = 0.0

  do i = il1,il2
       fdenit_sw(i) = 0.0
       ntsl(i) = 0
       soilwavg(i) = 0.0
    do k = 1,ignd
       soilw(i,k) = 0.0
    end do ! loop 110
    do j = 1,iccp1
       denitveg(i,j) = 0.0
       no_denitveg(i,j) = 0.0
       n2o_denitveg(i,j) = 0.0
       n2_denitveg(i,j) = 0.0
       fdenit_t(i,j) = 0.0
       no3_mass_in(i,j) = 0.0
       tempq10s(i,j) = 0.0
       soilcq10(i,j) = 0.0
    end do ! loop 120
  end do ! loop 100


  !> **A.** Find number of soil layers in which denitrification occurs: 
  do k = 1,ignd-1
   do i = il1,il2
      if (nitdenitdpth > zbotw(i,k) .and. nitdenitdpth <= zbotw(i,k+1)) ntsl(i) = k
   end do ! loop 140
  end do ! loop 130

  !> **B.** Using average soil temperature found in nitrification.f90,
  !! calculate soil temperature dependence over soil layers in which denitrification occurs:
  !! \f$ f(T) = Q_{10}^{(T_{0.5m}-20.)/10.} \f$ <br>
  !! where \f$ Q_{10}\f$ itself is found as a function of temperature following CTEM's logic used in the
  !! heterotrophic respiration module. <br>
  !! \f$ Q_{10}=1.44 +  0.56(tanh(0.075(46-T_{0.5m})))\f$ .
  do j = 1,iccp1
   do i = il1,il2
      tempq10s(i,j) = solctempavg(i) - 273.16
      soilcq10(i,j) = tanhq10(1) + tanhq10(2) * &
                      (tanh(tanhq10(3) * (tanhq10(4) - tempq10s(i,j))))
      fdenit_t(i,j) = soilcq10(i,j)**((solctempavg(i) - 273.16 - topt_denit) / 10.)
   end do ! loop 160
  end do ! loop 150

  !> **C.** Calculate volumetric soil moisture dependence over soil layers in which denitrification occurs:
  !! (see equation 25 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! Soil wetness as a function of soil moisture, <br>
  !! \f$ w(\theta) = \max(0, \min(1, \frac{(\theta - \theta_{w})}{(\theta_{f}-\theta_{w})} )) \f$ <br>
  !! Dependence on soil wetness, <br>
  !! \f$ f(\theta) = 1 - tanh( 2.5 ( \frac{1-w(\theta)}{1 - w_{d}} )^{2} ) \f$ <br>
  !! where \f$ w_{d} \f$ is the soil wetness threshold below which denitification doesn't occur.
  do i = il1,il2
    do k = 1,ignd
      if (isand(i,k) /= -4 .and. isand(i,k) /= -3) then !exclude ice & bedrock
         soilw(i,k) = (thliq(i,k) - THLW(i,k)) / (THFC(i,k) - THLW(i,k))
         soilw(i,k) = max(0.,min(1., soilw(i,k)))
      endif
    end do ! loop 180
       soilwavg(i) = sum(soilw(i,1:ntsl(i))) / ntsl(i)
       fdenit_sw(i) = 1.-tanh(2.5* ((1 - soilwavg(i)) / (1. - denitswthrsh))**2)
  end do ! loop 170

  !> **D.** With the soil moisture and soil temperature scalars calculated the 
  !!  \f$ NO, N_{2}O \f$ and \f$ N_{2} \f$ gaseous losses can now be found. <br><br>
  !!  
  !! \f$ E_{NO} = \mu_{NO} f(T)  f(\theta) N_{NO3} \f$ <br>
  !! \f$ E_{N2O} = \mu_{N2O} f(T)  f(\theta) N_{NO3} \f$       <br>
  !! \f$ E_{N2} = \mu_{N2} f(T)  f(\theta) N_{NO3} \f$ <br>
  !!
  !! <br> where \f$\mu_{NO}, \mu_{N2O}, \f$ and \f$ \mu_{N2} \f$ are coefficients. In the model
  !! code, however, things are done slighly differently. First a total denitification flux is calculated which
  !! is then divided into gaseous losses associated with \f$ NO, N_{2}O \f$ and \f$ N_{2}\f$.
  !! (see equation 24 and 25 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  do i = il1,il2
    do j = 1,icc
      
      ! PFTs
      denitveg(i,j) = mu_denit * fdenit_t(i,j) * fdenit_sw(i) * no3_mass(i,j) ! gross denitrifcation

      if (denitveg(i,j) < 0.0) then
         write(*,*)'denitveg lt zero at i=',i,' for pft=',j
         write(*,*)'denitveg = ',denitveg(i,j)
         call errorHandler('denitrification', - 1)
      end if
      no_denitveg(i,j) = mu_no_denit * denitveg(i,j)   ! NO emission during denitrification 
      n2o_denitveg(i,j) = mu_n2o_denit * denitveg(i,j) ! N2O emission during denitrification
      n2_denitveg(i,j) = denitveg(i,j) - no_denitveg(i,j) - n2o_denitveg(i,j)  ! N2 emission during denitrification

    end do ! loop 200
    
    ! bare soil 
    denitveg(i,iccp1) = mu_denit * fdenit_t(i,iccp1) * fdenit_sw(i) * no3_mass(i,iccp1)

    if (denitveg(i,iccp1) < 0.0) then
       write(*,*)'denitveg lt zero at i=',i,' for pft=',iccp1
       write(*,*)'denitveg = ',denitveg(i,iccp1)
       call errorHandler('denitrification', - 2)
    end if

    no_denitveg(i,iccp1) = mu_no_denit * denitveg(i,iccp1)
    n2o_denitveg(i,iccp1) = mu_n2o_denit * denitveg(i,iccp1)
    n2_denitveg(i,iccp1) = denitveg(i,iccp1) - no_denitveg(i,iccp1) - n2o_denitveg(i,iccp1)

  end do ! loop 190

  !> **E.** Update NO3- and check conservation:
  do j = 1,iccp1
     do i = il1,il2

        no3_mass_in(i,j) = no3_mass(i,j)
        no3_mass(i,j) = no3_mass_in(i,j) - denitveg(i,j)

        !! Adjusting to prevent negative no3_mass(i,j)
        if (no3_mass(i,j) < 0.0) then
           write(*,*)'no3_mass lt zero in denitrification.f90 at i=',i,' for pft=',j
           write(*,*)'no3_mass(i,j) = ',no3_mass(i,j)
           write(*,*)'denitveg(i,j) = ',denitveg(i,j)
           write(*,*)'adjusting to prevent negative no3_mass(i,j)'
           temp = denitveg(i,j)
           denitveg(i,j) = denitveg(i,j) - denitveg(i,j) * abs(no3_mass(i,j) / temp)
           no3_mass(i,j) = 0.0
           temp = 0.0
        end if

        !! Checking conservation for no3_mass(i,j)
        if ( (no3_mass_in(i,j) - no3_mass(i,j) - denitveg(i,j)) &
                                          >  ntolerance) then
           write(*,*)'imbalance in NO3- in denitrification  &
                      subroutine at (i)=',i,'pft=',j  !j=10 is bare soil
           call errorHandler('denitrification', - 3)
        end if

     end do ! loop 220
  end do ! loop 210

  return
end subroutine denitrific
