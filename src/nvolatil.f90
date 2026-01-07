!> \file
!! Calculates NH3 (ammonia) volatilization from NH4+.
!> @author A. Asaadi

subroutine volatilnh3(il1, il2, tbar, CFLUX_GA, USTARBS_GA, & ! In
                      soilpH, thliq, delzw, & ! In
                      nh4_mass, & ! In/Out
                      nvolveg) ! Out

  use classicParams, only : ilg, ignd, icc, ctempfts, iccp1, &
                            nh3atm, nvol_coeff, ntolerance

  implicit none

  integer,intent(in) :: il1       !< il1=1
  integer,intent(in) :: il2       !< il2=ilg (no. of grid cells in latitude circle)
  real,dimension(ilg,ignd),intent(in) :: tbar          !< soil temperature, K 
  real,dimension(ilg),intent(in) :: soilpH             !< soil PH
  real,dimension(ilg),intent(in) :: CFLUX_GA           !< aerodynamic conductance, inverse of aerodynamic resistance, \f$m s^{-1}\f$ 
  real,dimension(ilg),intent(in) :: USTARBS_GA         !< friction velocity, \f$m s^{-1}\f$ 
  real,dimension(ilg,ignd),intent(in) :: thliq         !< volumetric soil moisture, \f$m^3 m^{-3}\f$
  real,dimension(ilg,ignd),intent(in) :: delzw         !< thickness of each soil layer,\f$m\f$
  real,dimension(ilg,iccp1),intent(inout) :: nh4_mass  !< NH4+ for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,dimension(ilg,iccp1),intent(out) :: nvolveg     !< NH3 volatilization loss from NH4+, \f$g N m^{-2} day^{-1}\f$ 

  ! local variables
  integer :: i, j
  real :: temp
  real,dimension(ilg) :: rb                 !< boundary layer resistance, \f$s m^{-1}\f$
  real,dimension(ilg) :: ra                 !< aerodynamic resistance, \f$s m^{-1}\f$
  real,dimension(ilg) :: KH                 !< Henry's law constant for NH3, following Riddick et al. (2016) 
  real,dimension(ilg) :: KNH4               !< dissociation equilibrium constant for NH4, following Riddick et al. (2016), \f$mol L^{-1}\f$
  real,dimension(ilg) :: D                  !< denominator variable 
  real,dimension(ilg,iccp1) :: nh3g         !< NH3 concentration at the soil-atmosphere interface, \f$g N m^{-3}\f$
  real,dimension(ilg,iccp1) :: nh4_mass_in  !< previous NH4+, \f$g N m^{-2}\f$
  !
  temp = 0.0


  !> **A.** Calculate K_H, K_NH4 and the denominator 1+ K_H + [H^+]*(K_H/K_NH4)
  !! (see equation 29, 30 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ \chi = 0.26 \frac{N_{NH_4}}{1+K_H+K_H[H^+] / K_{NH_4}} \f$
  !! \f$ K_H = 4.59 T_{0.1} exp(4092 (\frac{1}{T_{0.1}} - \frac{1}{T_{ref,v}})) \f$
  !! \f$ K_{NH_4} = 5.67 \times 10^{-10} exp(-6286(\frac{1}{T_{0.1}} - \frac{1}{T_{ref,v}})) \f$
  do i = il1,il2
     KH(i) = 4.59 * tbar(i,1)*exp(4092.* ( (1./tbar(i,1)) - (1./298.15) ) )
     KNH4(i) = 5.67*1.E-10*exp(-6286.*((1. / tbar(i,1)) - (1. / 298.15)))
     D(i) = 1.0 + KH(i) + 10.**(- soilpH(i)) * (KH(i) / KNH4(i))
     do j = 1,iccp1
        nvolveg(i,j) = 0.0
        nh4_mass_in(i,j) = 0.0
     end do ! loop 90
  end do ! loop 80

  do j = 1,iccp1
   do i = il1,il2

      !> **B.** Calculate r_b and r_a  
      !! (see equation 28 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
      !! \f$ r_b = 6.2u_*^{-0.67}\f$
      rb(i) = 6.2 * USTARBS_GA(i)**(-0.67)
      ra(i) = CFLUX_GA(i)**(-1.)

      !> **C.** Calculate chi
      !! (see equation 29 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
      !! \f$ \chi = 0.26 \frac{N_{NH_4}}{1+K_H+K_H[H^+] / K_{NH_4}} \f$
      nh3g(i,j) = 0.259*nh4_mass(i,j)/D(i)

      !> **D.** Calculate NH3 volatilization
      !! (see equation 27 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
      !! \f$ V_{NH_4} = \vartheta 86400 \frac{1}{r_a + r_b}(\chi - [NH_{3,a}]) \f$
      if (ra(i) + rb(i) > 0.0) then
        nvolveg(i,j) = (nh3g(i,j) - nh3atm) * (1. / (ra(i) + rb(i))) * 86400.
        nvolveg(i,j) = nvol_coeff * (1.0 - 0.7) * nvolveg(i,j) !fcapture=0.7 following Table A1 of Riddick et al. (2016)
        nvolveg(i,j) = min(max(nvolveg(i,j), 0.), nh4_mass(i,j))
      else
        nvolveg(i,j) = 0.0
      end if

      ! Updating NH4+
      nh4_mass_in(i,j) = nh4_mass(i,j)
      nh4_mass(i,j) = nh4_mass(i,j) - nvolveg(i,j)

      ! Adjusting to prevent negative nh4_mass(i,j)
      if (nh4_mass(i,j) < 0.0) then
          write(*,*)'nh4_mass lt zero in nvolatil.f90 at i=',i,' for pft=',j,''
          write(*,*)'nh4_mass(i,j) = ',nh4_mass(i,j)
          write(*,*)'nvolveg(i,j)  = ',nvolveg(i,j)
          write(*,*)'adjusting nvolveg(i,j) to prevent negative nh4_mass(i,j)'
          temp = nvolveg(i,j)
          nvolveg(i,j) = nvolveg(i,j) - nvolveg(i,j) * abs(nh4_mass(i,j)) / temp
          nh4_mass(i,j) = 0.0
          temp = 0.0
      end if

      ! Checking conservation of N
      if ( (nh4_mass_in(i,j) - nh4_mass(i,j) - nvolveg(i,j)) &
                                         > ntolerance) then
           write(*,*)'imbalance in NH4+ at (i)=',i,'pft=',j  !j=10 is bare soil
           call errorHandler('nvolatil', - 1)
      end if

  end do ! loop 110
 end do ! loop 100

  return
end subroutine volatilnh3
