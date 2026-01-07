!> \file
!! Calculates N leaching from NO3-. 
!> @author A. Asaadi

subroutine leaching(il1, il2, ROFB, & ! In
                    no3_mass, & ! In/Out
                    nleachveg) ! Out


  use classicParams, only : ilg, ignd, icc, iccp1, &
                            ctempfts, nleach_coeff, &
                            ntolerance

  implicit none

  integer,intent(in) :: il1       !< il1=1
  integer,intent(in) :: il2       !< il2=ilg (no. of grid cells in latitude circle)
  real,intent(in),dimension(ilg) :: ROFB                !< base flow from bottom of soil column, \f$kg H2O m^{-2} s^{-1}\f$ 
  real,intent(inout),dimension(ilg,iccp1) :: no3_mass   !< NO3- for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,intent(out),dimension(ilg,iccp1)   :: nleachveg  !< leaching from NO3-, \f$g N m^{-2} day^{-1}\f$

  ! local variables
  integer :: i, j
  real :: temp
  real,dimension(ilg,iccp1) :: no3_mass_in !< previous NO3- for individual PFTs + bare soil, \f$g N m^{-2}\f$

  do j = 1,iccp1
    do i = il1,il2
      nleachveg(i,j) = 0.0
      no3_mass_in(i,j) = 0.0
    end do ! loop 120
  end do ! loop 110

  do j = 1,iccp1
     do i = il1,il2

        !! Calculate NO3- leaching
        !! (see equation 26 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
        !! \f$ L_{NO_3} = 86400 \varphi b_t N_{NO_3} \f$
        nleachveg(i,j) = nleach_coeff * ROFB(i) * 86400. * no3_mass(i,j)

        if (nleachveg(i,j) .lt. 0.0) then
           write(*,*)'nleachveg lt zero at i=',i,' for pft=',j
           write(*,*)'nleachveg = ',nleachveg(i,j)
           call errorHandler('nleach', - 1)
        end if

        !! Updating NO3- and check conservation
        no3_mass_in(i,j) = no3_mass(i,j)
        no3_mass(i,j) = no3_mass(i,j) - nleachveg(i,j)

        !! Adjusting to prevent negative no3_mass(i,j)
        if (no3_mass(i,j) .lt. 0.0) then
           write(*,*)'no3_mass lt zero in nleach.f90 at i=',i,' for pft=',j
           write(*,*)'no3_mass(i,j)  = ',no3_mass(i,j)
           write(*,*)'nleachveg(i,j) = ',nleachveg(i,j)
           write(*,*)'adjusting to prevent negative no3_mass(i,j)'
           temp = nleachveg(i,j)
           nleachveg(i,j) = nleachveg(i,j) - nleachveg(i,j) &
                            * abs(no3_mass(i,j) / temp)
           no3_mass(i,j)  = 0.0
           temp = 0.0
        end if

        !! Checking conservation for no3_mass(i,j)
        if ( (no3_mass_in(i,j) - no3_mass(i,j) - nleachveg(i,j)) >  ntolerance) then
           write(*,*)'imbalance in NO3- in leaching &
                      subroutine at (i)=',i,'pft=',j  !j=10 is bareground
           call errorHandler('nleach', - 2)
        end if

     end do ! loop 140
  end do ! loop 130

  return
end subroutine leaching
