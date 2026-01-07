!> \file
!! Calculates nitrification.
!> @author A. Asaadi
!!


subroutine nitrific(il1, il2, zbotw, tbar, thliq, thice, thpor, & ! In
                    THFC, sort, isand, ipeatland, & ! In
                    psisat, bi,   & ! In
                    nh4_mass, no3_mass, & ! In/Out
                    solctempavg, nitrifveg, no_nitveg, n2o_nitveg) ! Out


  use classicParams, only : ilg, ignd, icc, iccp1, mu_nit, topt_nit, &
                            nitdenitdpth, mu_no_nit, mu_n2o_nit, &
                            ctempfts, tanhq10, ntolerance

  implicit none


  integer,intent(in) :: il1         !< il1=1
  integer,intent(in) :: il2         !< il2=ilg (no. of grid cells in latitude circle)
  integer,intent(in),dimension(icc) :: sort        !< index for correspondence between PFTs and the 12 values in parameters vectors
  integer,intent(in),dimension(ilg) :: ipeatland   !< Peatland flag: 0 = not a peatland, 1= bog, 2 = fen
  integer,intent(in),dimension(ilg,ignd) :: isand  !< flag for bedrock or ice in a soil layer
  real,intent(in),dimension(ilg,ignd) :: tbar      !< soil temperature, K 
  real,intent(in),dimension(ilg,ignd) :: thliq     !< volumetric soil moisture, \f$m^3 m^{-3}\f$ 
  real,intent(in),dimension(ilg,ignd) :: thice     !< frozen volumetric soil moisture, \f$m^3 m^{-3}\f$ 
  real,intent(in),dimension(ilg,ignd) :: thpor     !< soil total porosity, \f$m^3 m^{-3}\f$
  real,intent(in),dimension(ilg,ignd) :: THFC      !< volumetric soil moisture at field capacity, \f$m^3 m^{-3}\f$
  real,intent(in),dimension(ilg,ignd) :: zbotw     !< bottom of each soil layer, \f$m\f$
  real,intent(in),dimension(ilg,ignd) :: psisat    !< soil matric potential at saturation, Mpa
  real,intent(in),dimension(ilg,ignd) :: bi        !< parameter B of Clapp and Hornberger 1978 

  real,intent(inout),dimension(ilg,iccp1) :: nh4_mass  !< NH4+ for individual PFTs + bare soil, \f$g N m^{-2}\f$ 
  real,intent(inout),dimension(ilg,iccp1) :: no3_mass  !< NO3- for individual PFTs + bare soil, \f$g N m^{-2}\f$

  real,intent(out),dimension(ilg,iccp1) :: nitrifveg   !< nitrification flux from NH4+ to NO3- for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg,iccp1) :: no_nitveg   !< nitrification NO loss for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$ 
  real,intent(out),dimension(ilg,iccp1) :: n2o_nitveg  !< nitrification N2O loss for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  real,intent(out),dimension(ilg) :: solctempavg       !< average soil temperature over soil layers in which nitrification occurs, K 

  ! local variables
  integer :: i, j, k
  integer,dimension(ilg) :: ntsl             !< number of soil layers in which nitrification occurs
  real,dimension(ilg,iccp1) :: soilcq10      !< Q10 for nitrification 
  real,dimension(ilg,ignd) :: scmotrm        !< soil matric potential dependence for each soil layer in which nitrification occurs 
  real,dimension(ilg,iccp1) :: tempq10s      !< average soil temperature over soil layers in which nitrification occurs, C
  real,dimension(ilg,iccp1) :: fnit_t        !< soil temperature dependence over soil layers in which nitrification occurs 
  real,dimension(ilg) :: fnit_sw             !< average soil matric potential dependence over soil layers in which nitrification occurs
  real,dimension(ilg,ignd) :: psi            !< soil matric potential, MPa 
  real,dimension(ilg,iccp1) :: no3_mass_in   !< previous NO3- for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,dimension(ilg,iccp1) :: nh4_mass_in   !< previous NH4+ for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real :: temp

  temp = 0.0
  do i = il1,il2
    ntsl(i) = 0
    do k = 1,ignd
       scmotrm(i,k) = 0.0
    end do ! loop 105
    do j = 1,iccp1
       nitrifveg(i,j) = 0.0
       no_nitveg(i,j) = 0.0
       n2o_nitveg(i,j) = 0.0
       solctempavg(i) = 0.0
       tempq10s(i,j) = 0.0
       soilcq10(i,j) = 0.0
       no3_mass_in(i,j) = 0.0
       nh4_mass_in(i,j) = 0.0
    end do ! loop 115
  end do ! loop 110

  !> **A.** Find number of soil layers in which nitrification occurs:
  do k = 1,ignd-1
    do i = il1,il2
      if (nitdenitdpth > zbotw(i,k) .and. nitdenitdpth <= zbotw(i,k+1)) ntsl(i) = k
    end do ! loop 140
  end do ! loop 130

  !> **B.** Calculate average temperature over soil layers in which nitrification occurs
  do i = il1,il2
       solctempavg(i) = sum(tbar(i,1:ntsl(i))) / ntsl(i)
  end do ! loop 180

  !> **C.** Calculate soil temperature dependence over soil layers in which nitrification occurs: 
  !! (see equation 20 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ f_I(T_{0.5}) = Q_{10,I}^{\frac{T_{0.5}-20}{10}} \f$
  !! \f$ Q_{10,I} = 1.44 + 0.56(tanh(0.075(46-T_{0.5}))) \f$
  do j = 1,iccp1
     do i = il1,il2
        tempq10s(i,j) = solctempavg(i) - 273.16
        soilcq10(i,j) = tanhq10(1) + tanhq10(2) * (tanh(tanhq10(3) * (tanhq10(4) - tempq10s(i,j))))
        fnit_t(i,j) = soilcq10(i,j)**((tempq10s(i,j) - topt_nit) / 10.)
     end do ! loop 200
  end do ! loop 190

  !> **D.** Calculate soil matric potential dependence for each soil layer in which nitrification occurs:
  !! (see equation 21 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ f_I(\psi) = 0.5 \f$ if \f$ \psi \leq \psi_{sat} \f$
  !! \f$ f_I(\psi) = 1-0.5\frac{log(0.4)-log(\psi)}{log(0.4)-log(\psi_{sat})} \f$ if \f$ 0.4 > \psi \geq \psi_{sat} \f$
  !! \f$ f_I(\psi) = 1-0.8\frac{log(\psi)-log(0.6)}{log(100)-log(0.6)} \f$ if \f$ 100 > \psi > 0.6 \f$
  !! \f$ f_I(\psi) = 0.2 \f$ if \f$ \psi > 100 \f$
  do k = 1, ignd
     do i = il1, il2

        if (isand(i,k) == -3 .or. isand(i,k) == -4) then
           scmotrm (i,k) = 0.2
           psi (i,k) = 10000.0 ! arbitrary large number
        else
           if (ipeatland(i) > 0) then
             if (thliq(i,k) + thice(i,k) + 0.01 < thpor(i,k) .and. tbar(i,k) < 273.16) then
                psi(i,k) = 0.001
             else if ( thice(i,k) > thpor(i,k) ) then
                psi(i,k) = 0.001
             else
                psi(i,k) = psisat(i,k) * (thliq(i,k) / (thpor(i,k) + 0.005 - thice(i,k)))**(-bi(i,k))
             end if
           else
            if (thice(i,k) < thpor(i,k)) then
               psi(i,k) = psisat(i,k) * (thliq(i,k) / (thpor(i,k) + 0.005 - thice(i,k)))**(-bi(i,k))
            else
               psi(i,k) = 10000.0
            end if
           end if

           if (psi(i,k) >= 10000.0) then
              scmotrm(i,k) = 0.2
           else if( psi(i,k) < 10000.0 .and.  psi(i,k) > 6.0 ) then
              scmotrm(i,k) = 1.0 - 0.8 * ( (log10(psi(i,k)) - log10(6.0)) / (log10(10000.0) - log10(6.0)) )
           else if( psi(i,k) <= 6.0 .and. psi(i,k) >= 4.0 ) then
              scmotrm(i,k) = 1.0
           else if( psi(i,k) < 4.0 .and. psi(i,k) > psisat(i,k) )then
              scmotrm(i,k) = 1.0 - 0.5 * ( (log10(4.0) - log10(psi(i,k))) /(log10(4.0) - log10(psisat(i,k))) )
           else if( psi(i,k) <= psisat(i,k) ) then
              scmotrm(i,k) = 0.5
           end if

           scmotrm(i,k)= max(0.2,min(1.0,scmotrm(i,k)))

        end if

     end do ! loop 220
  end do ! loop 210
  !> Calculate average soil matric potential dependence over soil layers in which nitrification occurs:
  do i = il1,il2
       fnit_sw(i) = sum(scmotrm(i,1:ntsl(i))) / ntsl(i)
  end do ! loop 230

  !> **E.** Find net nitrification and NO and N2O emissions during nitrification:
  !! (see equation 19 and 23 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ I_{NO_3} = \eta f_I(T_{0.5}) f_I(\psi) N_{NH_4} \f$
  !! \f$ I_{NO} = \eta_{NO} I_{NO_3} \f$
  !! \f$ I_{N_2O} = \eta_{N_2O} I_{NO_3} \f$
  do i = il1,il2
    do j = 1,icc

      nitrifveg(i,j) = mu_nit * fnit_t(i,j) * fnit_sw(i)*nh4_mass(i,j)  !gross nitrification

      if (nitrifveg(i,j) < 0.0) then
         write(*,*)'nitrifveg lt zero at i=',i,' for pft=',j
         write(*,*)'nitrifveg = ',nitrifveg(i,j)
         call errorHandler('nitrification', - 1)
      end if

      no_nitveg(i,j) = mu_no_nit * nitrifveg(i,j)    ! NO emission during nitrification
      n2o_nitveg(i,j) = mu_n2o_nit * nitrifveg(i,j)  ! N2O emission during nitrification

      nitrifveg(i,j) = nitrifveg(i,j) - no_nitveg(i,j) - n2o_nitveg(i,j)  !net nitrification

    end do ! loop 250

    ! bare soil
    nitrifveg(i,iccp1) = min(max(mu_nit * fnit_t(i,iccp1) * fnit_sw(i) &
                             * nh4_mass(i,iccp1), 0.0), nh4_mass(i,iccp1))
    no_nitveg(i,iccp1) = mu_no_nit * nitrifveg(i,iccp1)
    n2o_nitveg(i,iccp1) = mu_n2o_nit * nitrifveg(i,iccp1)
    nitrifveg(i,iccp1) = nitrifveg(i,iccp1) - no_nitveg(i,iccp1) - n2o_nitveg(i,iccp1)

  end do ! loop 240

  !> **F.** Update NH4+ and NO3- and check conservation: 
  do j = 1,iccp1
   do i = il1,il2

      no3_mass_in(i,j) = no3_mass(i,j)
      nh4_mass_in(i,j) = nh4_mass(i,j)
      nh4_mass(i,j) = nh4_mass(i,j) - &
                      (nitrifveg(i,j) + no_nitveg(i,j) + n2o_nitveg(i,j))
      no3_mass(i,j) = no3_mass(i,j) + nitrifveg(i,j)

      !! Adjusting to prevent negative nh4_mass(i,j)
      if (nh4_mass(i,j) < 0.0) then
          write(*,*)'nh4_mass lt zero in nitrification.f90 at i=',i,' for pft=',j
          write(*,*)'nh4_mass(i,j)  = ',nh4_mass(i,j)
          write(*,*)'nitrifveg(i,j) = ',nitrifveg(i,j)
          write(*,*)'no_nitveg(i,j) = ',no_nitveg(i,j)
          write(*,*)'n2o_nitveg(i,j)= ',n2o_nitveg(i,j)
          write(*,*)'adjusting to prevent negative nh4_mass(i,j)'
          temp = nitrifveg(i,j) + no_nitveg(i,j) + n2o_nitveg(i,j)
          nitrifveg(i,j) = nitrifveg(i,j) - nitrifveg(i,j) &
                           * abs(nh4_mass(i,j) / temp)
          no_nitveg(i,j) = no_nitveg(i,j) - no_nitveg(i,j) &
                           * abs(nh4_mass(i,j) / temp)
          n2o_nitveg(i,j) = n2o_nitveg(i,j) - n2o_nitveg(i,j) &
                            * abs(nh4_mass(i,j) / temp)
          nh4_mass(i,j) = 0.0
          temp = 0.0
      end if

      !! Checking conservation for nh4_mass(i,j)
      if ( (nh4_mass_in(i,j) - nh4_mass(i,j) - &
           (nitrifveg(i,j) + no_nitveg(i,j) + n2o_nitveg(i,j)))  &
                            > ntolerance) then
           write(*,*)'imbalance in NH4+ in nitrification &
                      subroutine at (i)=',i,'pft=',j  !j=10 is bareground
           call errorHandler('nitrification', - 2)
      end if
      
      !! Adjusting to prevent negative no3_mass(i,j)
      if (no3_mass(i,j) < 0.0) then
          write(*,*)'no3_mass lt zero in nitrification.f90 at i=',i,' for pft=',j
          write(*,*)'no3_mass(i,j)  = ',no3_mass(i,j)
          write(*,*)'nitrifveg(i,j) = ',nitrifveg(i,j)
          call errorHandler('nitrification', - 3)
      end if

      !! Checking conservation for no3_mass(i,j)
      if ( (no3_mass_in(i,j) - no3_mass(i,j) + nitrifveg(i,j)) &
                                           > ntolerance) then
           write(*,*)'imbalance in NO3- in nitrification &
                      subroutine at (i)=',i,'pft=',j  !j=10 is bareground
           call errorHandler('nitrification', - 4)
      end if

   end do ! loop 270
  end do ! loop 260

  return
end subroutine nitrific
