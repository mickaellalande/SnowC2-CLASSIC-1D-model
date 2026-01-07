!> \file
!! Calculates natural and anthropogenic biological nitrogen fixation (BNF)
!> @author A. Asaadi

subroutine bnfix(il1, il2, spinfast, thliq, zbotw, tbar, isand, & ! In
                 THFC, THLW, fcancmx, soilcmas_bulk, & ! In
                 soilnmas, & ! In/Out
                 bnf_free) ! Out

  use classicParams, only : ilg, ignd, icc, iccp1, iccp2, bnfdpth, &
                            ntolerance, r_bnf_f

  implicit none

  integer,intent(in) :: il1               !< il1=1
  integer,intent(in) :: il2               !< il2=ilg (no. of grid cells in latitude circle)
  integer,intent(in) :: spinfast          !< spinup factor for soil carbon & nitrogen pools
  integer,intent(in),dimension(ilg,ignd) :: isand      !< flag for bedrock or ice in a soil layer
  real,intent(in),dimension(ilg,icc) :: fcancmx        !< max. fractional coverage of PFTs
  real,intent(in),dimension(ilg,ignd) :: tbar          !< soil temperature, K 
  real,intent(in),dimension(ilg,ignd) :: thliq         !< volumetric soil moisture, \f$m^3 m^{-3}\f$
  real,intent(in),dimension(ilg,ignd) :: zbotw         !< bottom of each soil layer, \f$m\f$ 
  real,intent(in),dimension(ilg,ignd) :: THFC          !< volumetric soil moisture at field capacity, \f$m^3 m^{-3}\f$
  real,intent(in),dimension(ilg,ignd) :: THLW          !< volumetric soil moisture at wilting point, \f$m^3 m^{-3}\f$
  real,intent(in),dimension(ilg,iccp2):: soilcmas_bulk  !< SOM carbon for each PFT + bare soil + LUC product, \f$kg C m^{-2}\f$
  real,intent(inout),dimension(ilg,iccp2) :: soilnmas  !< SOM nitrogen for PFTs + bare soil, \f$g N m^{-2}\f$
  real,intent(out),dimension(ilg,iccp1) :: bnf_free    !< free-living biological nitrogen fixation for PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
  
  ! local variables
  integer :: i, j, k
  integer,dimension(ilg)   :: ntsl           !< number of soil layers in which BNF occurs 
  real,dimension(ilg)      :: fbnf_swavg     !< average volumetric soil moisture dependence over soil layers in which BNF occurs 
  real,dimension(ilg)      :: fbnf_tavg      !< average temperature dependence over soil layers in which BNF occurs
  real,dimension(ilg)      :: barefrac       !< bare soil fractional coverage 
  real,dimension(ilg,ignd) :: fbnf_t         !< soil temperature dependence over soil layers in which BNF occurs 
  real,dimension(ilg,ignd) :: fbnf_sw        !< volumetric soil moisture dependence 
  real,dimension(ilg,iccp2):: soilnmas_in    !< initial SOM nitrogen for PFTs + bare soil, \f$g N m^{-2}\f$
  
  do i = il1,il2
    fbnf_swavg(i) = 0.0
    fbnf_tavg(i) = 0.0
    ntsl(i) = 0
    barefrac(i) = 1.0
    do j = 1,iccp1
       bnf_free(i,j) = 0.0
    end do ! loop 110
    do k = 1,ignd
       fbnf_sw(i,k) = 0.0
       fbnf_t(i,k) = 0.0
    end do ! loop 120
  end do ! loop 100

  !> **A.** Calculate bare soil fractional coverage of each grid cell: 
  !!
  do i = il1,il2
   do j = 1,icc
        barefrac(i) = barefrac(i) - fcancmx(i,j)
   end do ! loop 140
  end do ! loop 130

  !> **B.** Find number of soil layers in which free-living BNF occurs: 
  !!
  do k = 1,ignd-1
   do i = il1,il2
     if (bnfdpth > zbotw(i,k) .and. bnfdpth <= zbotw(i,k+1)) ntsl(i) = k
   end do ! loop 160
  end do ! loop 150

  !> **C.** Calculate volumetric soil moisture dependence over soil layers in which free-living BNF occurs: 
  !! 
  do i = il1,il2
   !do k = 1,ignd
   ! if (isand(i,k) /= -4 .and. isand(i,k) /= -3) then !exclude bedrock or ice
   !    fbnf_sw(i,k) = (thliq(i,k) - THLW(i,k)) /(THFC(i,k) - THLW(i,k))
   !    fbnf_sw(i,k) = max(0.,min(1., fbnf_sw(i,k)))
   ! endif
   !end do ! loop 180
    fbnf_swavg(i) = 1. !sum(fbnf_sw(i,1:ntsl(i))) / ntsl(i)
  end do ! loop 170

  !> **D.** Calculate temperature dependence over soil layers in which free-living BNF occurs:
  !!
  do i = il1,il2
    do k = 1,ignd
     if (isand(i,k) /= -4 .and. isand(i,k) /= -3) then !exclude bedrock or ice
        fbnf_t(i,k) = exp(-2.6 + 0.21 * (tbar(i,k) - 273.16) * (1 - 0.5 * (tbar(i,k) - 273.16)/24.4))
     endif
    end do
    fbnf_tavg(i) = sum(fbnf_t(i,1:ntsl(i))) / ntsl(i)
  end do
  
  !> **E.** Calculate free-living BNF 
  !!
  do j = 1,icc
    do i = il1,il2
      bnf_free(i,j) = r_bnf_f * fbnf_tavg(i) * fbnf_swavg(i) * soilcmas_bulk(i,j) 
    end do
  end do ! loop 250
  do i = il1,il2
     bnf_free(i,iccp1) = r_bnf_f * fbnf_tavg(i) * fbnf_swavg(i) * soilcmas_bulk(i,iccp1)
  end do

  !> **F.** Update SOM nitrogen for PFTs + bare soil and check conservation: 
  !!
  do j = 1,iccp1
   do i = il1,il2

     soilnmas_in(i,j) = soilnmas(i,j)
     soilnmas(i,j) = soilnmas(i,j) + real(spinfast) * bnf_free(i,j)

     if (soilnmas(i,j) < 0.0) then
        write(*,*)'soilnmas lt zero in bnf.f90 at i = ',i,' for pft=',j
        write(*,*)'soilnmas(i,j) = ',soilnmas(i,j)
        write(*,*)'bnf_free(i,j)   = ',bnf_free(i,j)
        call errorHandler('bnf', - 3)
     end if

     if ((soilnmas_in(i,j) - soilnmas(i,j) + bnf_free(i,j)) > ntolerance &
         .and. spinfast == 1) then
        write(*,*)'imbalance in soilnmas in bnf.f90 subroutine at (i)=',i,'pft=',j !j=10 is bare soil
        call errorHandler('bnf', - 4)
     end if

   end do ! loop 300
  end do ! loop 290

  return
end subroutine bnfix
