!> \file
!> Module for all snow albedo and transmissivity subroutines

module snowAlbedoTransmissMod

  implicit none

  public :: snowAlbedoTransmiss
  public :: snow_albval
  public :: snow_tranval

contains

  !> \ingroup snowAlbedoTransmissMod_snowAlbedoTransmiss
  !! @{
  !! Diagnoses snowpack visible and near-IR albedos given the
  !! all-wave albedo at the current time step. Calculates snowpack
  !! transmissivity for shortwave radiation.
  !> @author D. Verseghy, J. Cole, M. Lazare
  !!
  !! In subroutine snowAging, called at the end of waterBudgetDriver, the change of
  !! total snow albedo over the current time step is calculated using
  !! an empirical exponential decay function, which has different
  !! coefficients depending on whether the snow is dry or melting. In
  !! this subroutine, if the ISNOALB switch is set to 0, the visible and
  !! near-IR components of the snow albedo are diagnosed from the total albedo.
  !! According to the literature (Aguado, 1985 \cite Aguado1985-fv ; Robinson and Kukla, 1984; Dirmhirn and
  !! Eaton, 1975 \cite Dirmhirn1975-vx), the following represent typical snow albedos for
  !! fresh snow, old dry snow and melting snow:
  !!
  !! \f[
  !! \begin{array} { | l | c | c | c | }
  !! \hline
  !!              & \text{Total albedo}  & \text{Visible albedo} & \text{Near-IR albedo} \\ \hline
  !! \text{Fresh snow}   &      0.84     &      0.95      &      0.73      \\ \hline
  !! \text{Old dry snow} &      0.70     &      0.84      &      0.56      \\ \hline
  !! \text{Melting snow} &      0.50     &      0.62      &      0.38      \\ \hline
  !! \end{array}
  !! \f]
  !!
  !! The same decay function is assumed to apply to all three albedo
  !! ranges, so the relative location of the visible and near-IR
  !! albedos, \f$\alpha_{s, VIS}\f$ and \f$\alpha_{s, NIR}\f$, on the decay curve will be analogous
  !! to that of the total albedo, \f$\alpha_{s, T}\f$. Thus, for dry snow:
  !!
  !! \f$[\alpha_{s, VIS} - 0.84]/[0.95-0.84] = [\alpha_{s, T} - 0.70]/[0.84-0.70]\f$
  !! \f$[\alpha_{s, NIR} - 0.56]/[0.73-0.56] = [\alpha_{s, T} - 0.70]/[0.84-0.70]\f$
  !!
  !! or, simplifying:
  !!
  !! \f$\alpha_{s, VIS} = 0.7857 \alpha_{s, T} + 0.2900\f$
  !! \f$\alpha_{s, NIR} = 1.2142 \alpha_{s, T} - 0.2900\f$
  !!
  !! For melting snow:
  !!
  !! [\f$\alpha_{s, VIS} - 0.62]/[0.95-0.62] = [\alpha_{s, T} - 0.50]/[0.84-0.50]\f$
  !! [\f$\alpha_{s, NIR} - 0.38]/[0.73-0.38] = [\alpha_{s, T} - 0.50]/[0.84-0.50]\f$
  !!
  !! or, simplifying:
  !!
  !! \f$\alpha_{s, VIS} = 0.9706 \alpha_{s, T} + 0.1347\f$
  !! \f$\alpha_{s, NIR} = 1.0294 \alpha_{s, T} - 0.1347\f$
  !!
  !! The above calculations are performed if the flag IALS is set to
  !! zero. If IALS is set to one, indicating that assigned snow
  !! albedos are to be used instead of calculated values, \f$\alpha_{s, VIS}\f$ and
  !! \f$\alpha_{s, NIR}\f$ are set to the assigned values ASVDAT and ASIDAT
  !! respectively. The sub-canopy values of visible and near-IR albedo
  !! are currently set equal to the open snowpack values (this is
  !! expected to change if a canopy litterfall parametrization is
  !! developed).
  !!
  !! The transmissivity of snow under vegetation \f$\tau_{s, c}\f$ is
  !! then calculated from the snow depth
  !! ZSNOW using Beer’s law, with an empirical extinction coefficient
  !! of \f$25.0 m^{-1}\f$ derived from the literature (Grenfell and Maykut,
  !! 1977 \cite Grenfell1977-pi ; Thomas, 1963):
  !!
  !! \f$\tau_{s, c} = exp[-25.0 z_s]\f$
  !!
  !! If the ISNOALB switch is set to zero, the value of ALSNO in the first
  !! wavelength band is set to the previously calculated value of \f$\alpha_{s, VIS}\f$
  !! and the values for the remaining bands are set to the previously calculated value
  !! of \f$\alpha_{s, NIR}\f$; the value of \f$\tau_{s, g}\f$  is set to \f$\tau_{s, c}\f$.
  !! If the ISNOALB switch is set to 1, a new parameterization for the snow albedo and
  !! transmissivity in four shortwave radiation bands (one visible and three near-IR)
  !! is used, according to Cole et al. (2017).  This parameterization incorporates
  !! the effects of snow grain size and black carbon content, and makes use of lookup
  !! tables contained in the CCCma subroutines SNOW_ALBVAL and SNOW_TRANVAL.
  subroutine snowAlbedoTransmiss (ALVSSN, ALIRSN, ALVSSC, ALIRSC, ALBSNO, & ! Formerly SNOWALBA
                                  TRSNOWC, ALSNO, TRSNOWG, FSDB, FSFB, RHOSNO, &
                                  REFSN, BCSN, SNO, CSZ, ZSNOW, FSNOW, ASVDAT, ASIDAT, &
                                  ALVSG, ALIRG, &
                                  ILG, IG, IL1, IL2, JL, IALS, NBS, ISNOALB)
    !
    !     * JAN 27/16 - D.VERSEGHY. REFINE CALCULATIONS OF ALVSSN AND ALIRSN.
    !     * NOV 16/13 - J.COLE.     Final version for gcm17:
    !     *                         - Fixes to get the proper BC mixing ratio in
    !     *                           snow, which required passing in and using
    !     *                           the snow density RHON.
    !     * JUN 22/13 - J.COLE.     ADD CODE FOR "ISNOALB" OPTION,
    !     *                         WHICH IS BASED ON 4-BAND SOLAR.
    !     * FEB 05/07 - D.VERSEGHY. STREAMLINE CALCULATIONS OF
    !     *                         ALVSSN AND ALIRSN.
    !     * APR 13/06 - D.VERSEGHY. SEPARATE ALBEDOS FOR OPEN AND
    !     *                         CANOPY-COVERED SNOW.
    !     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
    !     * MAR 18/02 - D.VERSEGHY. UPDATES TO ALLOW ASSIGNMENT OF
    !     *                         USER-SPECIFIED VALUES TO SNOW
    !     *                         ALBEDO.
    !     * JUN 05/97 - D.VERSEGHY. CLASS - VERSION 2.7.
    !     *                         SPECIFY LOCATION OF ICE SHEETS
    !     *                         BY SOIL TEXTURE ARRAY RATHER
    !     *                         THAN BY SOIL COLOUR INDEX.
    !     * NOV 29/94 - M.LAZARE.   CLASS - VERSION 2.3.
    !     *                         CALL ABORT CHANGED TO CALL errorHandler TO
    !     *                         ENABLE RUNNING ON PC'S.
    !     * MAR 13/92 - M.LAZARE.   CODE FOR MODEL VERSION GCM7 -
    !     *                         DIVIDE PREVIOUS SUBROUTINE
    !     *                         "SNOALB" INTO "snowAlbedoTransmiss" AND
    !     *                         "snowAging" AND VECTORIZE.
    !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
    !     *                         CLASS VERSION 2.0 (WITH CANOPY).
    !     * APR 11/89 - D.VERSEGHY. DISAGGREGATE SNOW ALBEDO INTO
    !     *                         VISIBLE AND NEAR-IR PORTIONS
    !     *                         CALCULATE TRANSMISSIVITY TO
    !     *                         SHORTWAVE RADIATION.
    !
    implicit none
    !
    !     * INTEGER CONSTANTS.
    !
    integer, intent(in) :: ILG, IG, IL1, IL2, JL, IALS
    integer, intent(in) :: NBS       !< Number of modelled shortwave radiation wavelength bands
    integer, intent(in) :: ISNOALB   !< Switch to model snow albedo in two or more wavelength bands
    integer             :: IPTBAD, I, IB
    !
    !     * OUTPUT ARRAYS.
    !
    real, intent(out)   :: ALVSSC(ILG)  !< Visible albedo of snow on ground under vegetation canopy [ ]
    real, intent(out)   :: ALIRSC(ILG)  !< Near-IR albedo of snow on ground under vegetation canopy [ ]
    real, intent(out)   :: ALSNO (ILG,NBS) !< Albedo of snow in each modelled wavelength band  [  ]
    real, intent(out)   :: TRSNOWG(ILG,NBS)!< Transmissivity of snow in bare areas to shortwave radiation \f$[ ] (\tau_{s,g})\f$
    !
    real, intent(inout) :: TRSNOWC(ILG) !< Transmissivity of snow under vegetation to shortwave radiation
    !< \f$[ ] (\tau_{s, c})\f$
    real, intent(inout) :: ALVSSN(ILG)  !< Visible albedo of snow pack on bare ground \f$[ ]
    !< (alpha_{s, VIS})\f$
    real, intent(inout) :: ALIRSN(ILG)  !< Near-IR albedo of snow pack on bare ground \f$[ ]
    !< (alpha_{s, NIR})\f$
    !
    !     * INPUT ARRAYS.
    !
    real, intent(in)    :: ALVSG  (ILG)    !< Near-IR albedo of bare ground  [  ]
    real, intent(in)    :: ALIRG  (ILG)    !< Visible albedo of bare ground  [  ]
    real, intent(in)    :: FSDB(ILG,NBS)!< Direct solar radiation in each modelled wavelength band \f$[W m^{-2}]\f$
    real, intent(in)    :: FSFB(ILG,NBS)!< Diffuse solar radiation in each modelled wavelength band \f$[W m^{-2}]\f$
    real, intent(in)    :: ZSNOW (ILG)  !< Depth of snow \f$[m] (z_s)\f$
    real, intent(in)    :: FSNOW (ILG)  !< Fractional coverage of snow on grid cell [ ]
    real, intent(in)    :: ASVDAT(ILG)  !< Assigned value of visible albedo of snow pack –
    !< optional [ ]
    real, intent(in)    :: ASIDAT(ILG)  !< Assigned value of near-IR albedo of snow pack –
    !< optional [ ]
    real, intent(in)    :: REFSN (ILG)  !< Snow grain size  [m]
    real, intent(in)    :: BCSN  (ILG)  !< Black carbon mixing ratio \f$[kg m^{-3}]\f$
    real, intent(in)    :: CSZ   (ILG)  !< Cosine of solar zenith angle  [  ]
    real, intent(in)    :: SNO   (ILG)  !< Mass of snow pack \f$[kg m^{-2}]\f$
    real, intent(in)    :: RHOSNO(ILG)  !< Density of snow pack  \f$[kg m^{-3}]\f$
    real, intent(inout) :: ALBSNO(ILG)  !< All-wave albedo of snow pack \f$[ ] (\alpha_{s, T})\f$

    !
    !     * LOCAL ARRAYS
    !
    real    :: SALBG(ILG,NBS), ALDIR(ILG,NBS), ALDIF (ILG,NBS), &
              TRDIR(ILG,NBS), TRDIF (ILG,NBS), REFSNO(ILG), BCSNO(ILG)
    integer :: C_FLAG(ILG)
    !
    !     * CONSTANTS.
    !
    real    :: WDIRCT, WDIFF
    integer :: SUM_C_FLAG
    !------------------------------------------------------------------

    IPTBAD = 0
    do I = IL1,IL2 ! loop 100
      if (ALBSNO(I) < 0.50 .and. ALBSNO(I) > 0.499) ALBSNO(I) = 0.50
      if (FSNOW(I) > 0.0 .and. IALS == 0) then
        if (ALBSNO(I) > 0.70) then
          ALVSSN(I) = 0.7857 * ALBSNO(I) + 0.2900
          ALIRSN(I) = 1.2142 * ALBSNO(I) - 0.2900
        else
          ALVSSN(I) = 0.9706 * ALBSNO(I) + 0.1347
          ALIRSN(I) = 1.0294 * ALBSNO(I) - 0.1347
        end if
        if (ALVSSN(I) > 0.999 .or. ALVSSN(I) < 0.001) IPTBAD = I
        if (ALIRSN(I) > 0.999 .or. ALIRSN(I) < 0.001) IPTBAD = I
      else if (FSNOW(I) > 0.0 .and. IALS == 1) then
        ALVSSN(I) = ASVDAT(I)
        ALIRSN(I) = ASIDAT(I)
      end if
      ALVSSC(I) = ALVSSN(I)
      ALIRSC(I) = ALIRSN(I)
      TRSNOWC(I) = EXP( - 25.0 * ZSNOW(I))
    end do ! loop 100
    !
    if (IPTBAD /= 0) then
      write(6,6100) IPTBAD,JL,ALVSSN(IPTBAD),ALIRSN(IPTBAD)
  6100 format('0AT (I,J) = (',I3,',',I3,'),ALVSSN,ALIRSN = ',2F10.5)
      call errorHandler('snowAlbedoTransmiss', - 1)
    end if
    !
    if (ISNOALB == 0) then
      do I = IL1,IL2
        ALSNO(I,1) = ALVSSN(I)
        ALSNO(I,2) = ALIRSN(I)
        ALSNO(I,3) = ALIRSN(I)
        ALSNO(I,4) = ALIRSN(I)

        TRSNOWG(I,1:NBS) = TRSNOWC(I)
      end do ! I
    else if (ISNOALB == 1) then
      do IB = 1,NBS
        do I = IL1,IL2
          if (IB == 1) then
            SALBG(I,IB) = ALVSG(I)
            ALSNO(I,IB) = ALVSSN(I)
          else
            SALBG(I,IB) = ALIRG(I)
            ALSNO(I,IB) = ALIRSN(I)
          end if
        end do ! I
      end do ! IB
      SUM_C_FLAG = 0
      do I = IL1,IL2
        if (ZSNOW(I) > 0.0) then
          C_FLAG(I) = 1
        else
          C_FLAG(I) = 0
        end if
        SUM_C_FLAG = SUM_C_FLAG + C_FLAG(I)
      end do ! I

      if (IALS == 0) then
        if (SUM_C_FLAG > 0) then
          !>
          !! Convert the units of the snow grain size and BC mixing ratio
          !! Snow grain size from meters to microns and BC from \f$kg BC/m^3\f$ to ng BC/kg SNOW
          !!
          do I = IL1,IL2
            if (C_FLAG(I) == 1) then
              REFSNO(I) = REFSN(I) * 1.0E6
              BCSNO(I)  = (BCSN(I) / RHOSNO(I)) * 1.0E12
            end if
          end do ! I

          call SNOW_ALBVAL(ALDIF, & ! OUTPUT
                          ALDIR, &
                          CSZ, & ! INPUT
                          SALBG, &
                          BCSNO, &
                          REFSNO, &
                          SNO, &
                          C_FLAG, &
                          IL1, &
                          IL2, &
                          ILG, &
                          NBS)

          call SNOW_TRANVAL(TRDIF, & ! OUTPUT
                            TRDIR, &
                            CSZ, &  ! INPUT
                            SALBG, &
                            BCSNO, &
                            REFSNO, &
                            SNO, &
                            C_FLAG, &
                            IL1, &
                            IL2, &
                            ILG, &
                            NBS)

          do IB = 1,NBS
            do I = IL1,IL2
              if (C_FLAG(I) == 1) then
                WDIRCT = FSDB(I,IB) &
                        / (FSDB(I,IB) + FSFB(I,IB) + 1.E-10)
                WDIFF  = 1.0 - WDIRCT
                ALSNO(I,IB) = ALDIF (I,IB) * WDIFF &
                              + ALDIR(I,IB) * WDIRCT
                TRSNOWG(I,IB) = TRDIF (I,IB) * WDIFF &
                                + TRDIR(I,IB) * WDIRCT
              end if ! C_FLAG
            end do ! I
          end do ! IB
        else ! SUM_C_FLAG == 0
          do I = IL1,IL2
            ALSNO(I,1)     = ALVSSN(I)
            ALSNO(I,2:NBS) = ALIRSN(I)
            TRSNOWG(I,1:NBS) = TRSNOWC(I)
          end do ! I
        end if ! SUM_C_FLAG
      else if (IALS == 1) then
        do I = IL1,IL2
          ALSNO(I,1) = ASVDAT(I)
          ALSNO(I,2) = ASIDAT(I)
          ALSNO(I,3) = ASIDAT(I)
          ALSNO(I,4) = ASIDAT(I)

          TRSNOWG(I,1:NBS) = TRSNOWC(I)
        end do ! I
      end if ! IALS
    end if ! ISNOALB
    return
  end subroutine snowAlbedoTransmiss
  !! @}
  ! ------------------------------------------------------------------------------
  !> \ingroup snowAlbedoTransmissMod_snow_albval
  !! @{
  !> \brief Computes the direct and diffuse snow albedo using a
  !! lookup table and information about the current snow pack state.
  !!
  !! Albedos are computed for each solar radiation wavelength intervals
  !! so a total of 8 albedos will be returned. These albedos can then be
  !! used to compute the total snow albedo based on the by weighting
  !! the results by the direct beam fraction of the incident solar radiation.
  !!
  !! @author J. Cole
  !!
  subroutine snow_albval(albdif, & ! output
                        albdir, &
                        smu, & ! input
                        salb, &
                        bc_conc, &
                        snow_reff, &
                        swe, &
                        c_ind, &
                        il1, &
                        il2, &
                        ilg, &
                        nbnd)
    !
    !     * feb 10/2015 - j.cole. new version for gcm18:
    !                             - nbc increased from 12 to 20.
    !                               therefor lbc_conc data statement
    !                               changed accordingly.
    !     * jan 24/2013 - j.cole. previous version for gcm17:
    !                    - computes the direct and diffuse snow albedo
    !                      using lookup table and current snow conditions.
    !
    use classicParams, only : albdif_lut, albdir_lut
    implicit none
    !
    ! input
    !
    real, intent(in), dimension(ilg) :: smu !< cosine of the solar zenith angle [unitless]
    real, intent(in), dimension(ilg) :: bc_conc !< concentration of black carbon in the snow pack \f$[ng (bc)/kg (snow)]\f$
    real, intent(in), dimension(ilg) :: snow_reff !< effective radius of the snow grain [microns]
    real, intent(in), dimension(ilg) :: swe !< snow water equivalent (snowpack density*snow pack depth) \f$[kg/m^2]\f$

    real, intent(in), dimension(ilg,nbnd) :: salb !< albedo of the underlying surface [unitless]
    integer, intent(in), dimension(ilg) :: c_ind !< indicator that a calculation should be performed for this point
                                                 !! 1-yes, 0-no
    integer, intent(in) :: il1   !< Index of first atmospheric column for calculations \f$[unitless]\f$
    integer, intent(in) :: il2   !< Index of last atmospheric column for calculations \f$[unitless]\f$
    integer, intent(in) :: ilg   !< Total number of atmospheric columns \f$[unitless]\f$
    integer, intent(in) :: nbnd  !<  number of wavelength intervals for which to compute the albedos
    !
    ! output
    !
    real, intent(out), dimension(ilg,nbnd) :: albdif !< diffuse snow albedo (aka white sky albedo) [unitless]
    real, intent(out), dimension(ilg,nbnd) :: albdir !< direct beam snow albedo (aka black sky albedo) [unitless]
    !==================================================================
    ! physical (adjustable) parameters
    !
    ! define and document here any adjustable parameters.
    ! this should be variable described using the doxygen format above as
    ! well as a description of its minimum/default/maximum.
    !
    ! here is an example,
    !
    ! real :: beta !< This is the adjustable factor for computing liquid cloud effective radius \f$[]\f$
    !           !! It is compute differently when using bulk or PAM aerosols.
    !           !! For bulk aerosols its minimum/default/maximum is (1.0/1.3/1.5).
    !==================================================================

    !
    ! local
    !
    real, dimension(ilg,2) :: wsmu
    real, dimension(ilg,2) :: wbc
    real, dimension(ilg,2) :: wreff
    real, dimension(ilg,2) :: wswe

    real, dimension(2) :: wsalb

    real :: wtt
    real :: snow_reff_l

    integer, dimension(ilg) :: ismu
    integer, dimension(ilg) :: ibc
    integer, dimension(ilg) :: ireff
    integer, dimension(ilg) :: iswe

    integer :: ib
    integer :: i
    integer :: isalb

    integer :: iismu
    integer :: iisalb
    integer :: iibc
    integer :: iireff
    integer :: iiswe

    integer :: mvidx

    !
    ! constants
    !
    integer, parameter :: nsmu = 10
    integer, parameter :: nsalb    = 11
    integer, parameter :: nbc      = 20
    integer, parameter :: nreff    = 10
    integer, parameter :: nswe     = 11
    integer, parameter :: nbnd_lut = 4

    !ignoreLint(6)
    real, parameter :: lsalb (nsalb) = (/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
    real, parameter :: lsmu(nsmu) = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
    real, parameter :: lsnow_reff(nreff) = (/50.0,75.0,100.0,150.0,200.0,275.0,375.0,500.0,700.0,1000.0/)
    real, parameter :: lswe(nswe) = (/0.1,0.25,0.65,1.7,4.4,12.0,30.0,75.0,200.0,500.0,5000.0/)
    real, parameter :: lbc_conc(nbc) = (/0.0,1.0,5.0,10.0,50.0,100.0,500.0,1000.0,5000.0,10000.0,50000.0,100000.0, &
                                      250000.0,500000.0,750000.0,1000000.0,2500000.0,5000000.0,7500000.0,10000000.0/)

    !real, dimension(nbc,nswe,nreff,nsmu,nsalb,nbnd_lut) :: albdif_lut ! offline, kept in classicParams
    !real, dimension(nbc,nswe,nreff,nsmu,nsalb,nbnd_lut) :: albdir_lut ! offline, kept in classicParams

    !integer :: snow_alb_lut_init ! Not used offline

    !common /snowalblut/ albdif_lut,albdir_lut,snow_alb_lut_init ! offline, kept in classicParams

    ! abort if the lut has not been read in
    ! if (snow_alb_lut_init /= 1) then
    !   write(6, * ) 'SNOW ALBEDO LUT HAS NOT BEEN INITIALIZED', &
    !             snow_alb_lut_init
    !   call errorHandler('SNOW_ALBVAL', - 1)
    ! end if

    ! abort if the number of bands in the lut does not match size passed in
    if (nbnd_lut /= nbnd) then
      write(6, * ) 'MISMATCH IN NUMBER OF WAVELENGTH INTERVALS'
      call errorHandler('SNOW_ALBVAL', - 2)
    end if

    ! compute the albedos using linear interpolation

    ! compute the interpolation weights and points once and reuse for
    ! albedo interpolation for each band.
    ! have a check to set the weights depending if the input is
    ! outside or inside the lookup table range

    do i = il1, il2
      if (c_ind(i) == 1) then
        snow_reff_l = snow_reff(i)
        ismu(i)     = mvidx(lsmu,       nsmu,  smu(i))
        ibc(i)      = mvidx(lbc_conc,   nbc,   bc_conc(i))
        ireff(i)    = mvidx(lsnow_reff, nreff, snow_reff(i))
        iswe(i)     = mvidx(lswe,       nswe,  swe(i))

        if (smu(i) <= lsmu(1)) then
          wsmu(i,2) = 0.0
          wsmu(i,1) = 1.0 - wsmu(i,2)
        else if (smu(i) > lsmu(nsmu)) then
          wsmu(i,2) = 1.0
          wsmu(i,1) = 1.0 - wsmu(i,2)
        else
          wsmu(i,2) = (smu(i) - lsmu(ismu(i))) &
                      / (lsmu(ismu(i) + 1) - lsmu(ismu(i)))
          wsmu(i,1) = 1.0 - wsmu(i,2)
        end if

        if (bc_conc(i) <= lbc_conc(1)) then
          wbc(i,2) = 0.0
          wbc(i,1) = 1.0 - wbc(i,2)
        else if (bc_conc(i) > lbc_conc(nbc)) then
          wbc(i,2) = 1.0
          wbc(i,1) = 1.0 - wbc(i,2)
        else
          wbc(i,2) = (bc_conc(i) - lbc_conc(ibc(i))) &
                    / (lbc_conc(ibc(i) + 1) - lbc_conc(ibc(i)))
          wbc(i,1) = 1.0 - wbc(i,2)
        end if

        if (snow_reff_l <= lsnow_reff(1)) then
          wreff(i,2) = 0.0
          wreff(i,1) = 1.0 - wreff(i,2)
        else if (snow_reff_l > lsnow_reff(nreff)) then
          wreff(i,2) = 1.0
          wreff(i,1) = 1.0 - wreff(i,2)
        else
          wreff(i,2) = (snow_reff_l - lsnow_reff(ireff(i))) &
                      / (lsnow_reff(ireff(i) + 1) &
                      - lsnow_reff(ireff(i)))
          wreff(i,1) = 1.0 - wreff(i,2)
        end if

        if (swe(i) <= lswe(1)) then
          wswe(i,2) = 0.0
          wswe(i,1) = 1.0 - wswe(i,2)
        else if (swe(i) > lswe(nswe)) then
          wswe(i,2) = 1.0
          wswe(i,1) = 1.0 - wswe(i,2)
        else
          wswe(i,2) = (swe(i) - lswe(iswe(i))) &
                      / (lswe(iswe(i) + 1) - lswe(iswe(i)))
          wswe(i,1) = 1.0 - wswe(i,2)
        end if
      end if
    end do ! i

    do ib = 1, nbnd
      do i = il1, il2
        if (c_ind(i) == 1) then

          isalb = mvidx(lsalb,    nsalb, salb(i,ib))

          if (salb(i,ib) <= lsalb(1)) then
            wsalb(2) = 0.0
            wsalb(1) = 1.0 - wsalb(2)
          else if (salb(i,ib) > lsalb(nsalb)) then
            wsalb(2) = 1.0
            wsalb(1) = 1.0 - wsalb(2)
          else
            wsalb(2) = (salb(i,ib) - lsalb(isalb)) &
                      / (lsalb(isalb + 1) - lsalb(isalb))
            wsalb(1) = 1.0 - wsalb(2)
          end if

          albdir(i,ib) = 0.0
          albdif (i,ib) = 0.0

          do iisalb = isalb,isalb + 1
            do iismu = ismu(i),ismu(i) + 1
              do iireff = ireff(i),ireff(i) + 1
                do iiswe = iswe(i), iswe(i) + 1
                  do iibc = ibc(i), ibc(i) + 1

                    wtt = wsmu(i,iismu - ismu(i) + 1) &
                          * wreff(i,iireff - ireff(i) + 1) &
                          * wswe(i,iiswe - iswe(i) + 1) &
                          * wbc(i,iibc - ibc(i) + 1) &
                          * wsalb(iisalb - isalb + 1)

                    albdif (i,ib) = albdif (i,ib) + wtt &
                                    * albdif_lut(iibc,iiswe,iireff,iismu,iisalb,ib)
                    albdir(i,ib) = albdir(i,ib) + wtt &
                                  * albdir_lut(iibc,iiswe,iireff,iismu,iisalb,ib)

                  end do ! iibc
                end do  ! iiswe
              end do     ! iireff
            end do        ! iismu
          end do           ! iisalb

          if (albdif (i,ib) > 1.0 .or. albdif (i,ib) < 0.0) then
            write(6, * ) 'Bad albdif ',i,ib,smu(i),bc_conc(i), &
                        snow_reff(i),swe(i),salb(i,ib),albdif (i,ib)
            write(6, * ) i,ib,ismu(i),ibc(i),ireff(i),iswe(i),isalb
            call errorHandler('SNOW_ALBVAL', - 3)
          end if
          if (albdir(i,ib) > 1.0 .or. albdir(i,ib) < 0.0) then
            write(6, * ) 'Bad albdir ',i,ib,smu(i),bc_conc(i), &
                          snow_reff(i),swe(i),salb(i,ib),albdir(i,ib)
            write(6, * ) i,ib,ismu(i),ibc(i),ireff(i),iswe(i),isalb
            call errorHandler('SNOW_ALBVAL', - 3)
          end if
        else
          albdif (i,ib) = - 999.0
          albdir(i,ib) = - 999.0
        end if
      end do                 ! i
    end do                    ! ib

    return
  end subroutine snow_albval
  !! @}
  ! ------------------------------------------------------------------------------
  !> \ingroup snowAlbedoTransmissMod_snow_tranval
  !! @{
  !> \brief Computes the direct and diffuse snow transmission using a
  !! lookup table and information about the current snow pack state.
  !!
  !! Transmission are computed for each solar radiation wavelength intervals
  !! so a total of 8 albedos will be returned.  These transmissions can then be
  !! used to compute the total snow tramission by weighting
  !! the results by the direct beam fraction of the incident solar radiation.
  !! @author J. Cole
  !
  subroutine snow_tranval(trandif, & ! output
                          trandir, &
                          smu, & ! input
                          salb, &
                          bc_conc, &
                          snow_reff, &
                          swe, &
                          c_ind, &
                          il1, &
                          il2, &
                          ilg, &
                          nbnd)
    !
    !     * feb 10/2015 - j.cole. new version for gcm18:
    !                             - nbc increased from 12 to 20.
    !                               therefor lbc_conc data statement
    !                               changed accordingly.
    !     * jan 24/2013 - j.cole. previous version for gcm17:
    !                    - computes the direct and diffuse snow transmission
    !                      using lookup table and current snow conditions.
    !
    use classicParams, only : trandif_lut, trandir_lut

    implicit none
    !
    ! input
    !
    real, intent(in), dimension(ilg) :: smu !< cosine of the solar zenith angle [unitless]
    real, intent(in), dimension(ilg) :: bc_conc !< concentration of black carbon in the snow pack \f$[ng (bc)/kg (snow)]\f$
    real, intent(in), dimension(ilg) :: snow_reff !< effective radius of the snow grain [microns]
    real, intent(in), dimension(ilg) :: swe !< snow water equivalent (snowpack density*snow pack depth) \f$[kg/m^2]\f$

    real, intent(in), dimension(ilg,nbnd) :: salb !< albedo of the underlying surface [unitless]
    integer, intent(in), dimension(ilg) :: c_ind !< indicator that a calculation should be performed for this point
                                                 !! 1-yes, 0-no
    integer, intent(in) :: il1   !< Index of first atmospheric column for calculations \f$[unitless]\f$
    integer, intent(in) :: il2   !< Index of last atmospheric column for calculations \f$[unitless]\f$
    integer, intent(in) :: ilg   !< Total number of atmospheric columns \f$[unitless]\f$
    integer, intent(in) :: nbnd  !<  number of wavelength intervals for which to compute the albedos
    !
    ! output
    !
    real, intent(out), dimension(ilg,nbnd) :: trandif !< diffuse snow transmission (aka white sky transmission)
    real, intent(out), dimension(ilg,nbnd) :: trandir !< direct beam snow transmission (aka black sky transmission)
    !==================================================================
    ! physical (adjustable) parameters
    !
    ! define and document here any adjustable parameters.
    ! this should be variable described using the doxygen format above as
    ! well as a description of its minimum/default/maximum.
    !
    ! here is an example,
    !
    ! real :: beta !< This is the adjustable factor for computing liquid cloud effective radius \f$[]\f$
    !           !! It is compute differently when using bulk or PAM aerosols.
    !           !! For bulk aerosols its minimum/default/maximum is (1.0/1.3/1.5).
    !==================================================================

    !
    ! local
    !
    real, dimension(ilg,2) :: wsmu
    real, dimension(ilg,2) :: wbc
    real, dimension(ilg,2) :: wreff
    real, dimension(ilg,2) :: wswe

    real, dimension(2) :: wsalb

    real :: wtt

    integer, dimension(ilg) :: ismu
    integer, dimension(ilg) :: ibc
    integer, dimension(ilg) :: ireff
    integer, dimension(ilg) :: iswe

    integer :: ib
    integer :: i
    integer :: isalb

    integer :: iismu
    integer :: iisalb
    integer :: iibc
    integer :: iireff
    integer :: iiswe

    integer :: mvidx

    !
    ! constants
    !
    integer, parameter :: nsmu = 10
    integer, parameter :: nsalb    = 11
    integer, parameter :: nbc      = 20
    integer, parameter :: nreff    = 10
    integer, parameter :: nswe     = 11
    integer, parameter :: nbnd_lut = 4

    !ignoreLint(6)
    real, parameter :: lsalb (nsalb) = (/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
    real, parameter :: lsmu(nsmu) = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
    real, parameter :: lsnow_reff(nreff) = (/50.0,75.0,100.0,150.0,200.0,275.0,375.0,500.0,700.0,1000.0/)
    real, parameter :: lswe(nswe) = (/0.1,0.25,0.65,1.7,4.4,12.0,30.0,75.0,200.0,500.0,5000.0/)
    real, parameter :: lbc_conc(nbc) = (/0.0,1.0,5.0,10.0,50.0,100.0,500.0,1000.0,5000.0,10000.0,50000.0,100000.0, &
                                      250000.0,500000.0,750000.0,1000000.0,2500000.0,5000000.0,7500000.0,10000000.0/)

    !real, dimension(nbc,nswe,nreff,nsmu,nsalb,nbnd_lut) :: trandif_lut  ! offline, kept in classicParams
    !real, dimension(nbc,nswe,nreff,nsmu,nsalb,nbnd_lut) :: trandir_lut  ! offline, kept in classicParams

    !integer :: snow_tran_lut_init ! not used offline

    !common /snowtranlut/ trandif_lut,trandir_lut,snow_tran_lut_init

    ! abort if the lut has not been read in
    ! if (snow_tran_lut_init /= 1) then
    !   write(6, * ) 'SNOW TRANSMISSION LUT HAS NOT BEEN INITIALIZED', &
    !             snow_tran_lut_init
    !   call errorHandler('SNOW_TRANVAL', - 1)
    ! end if

    ! abort if the number of bands in the lut does not match size passed in
    if (nbnd_lut /= nbnd) then
      write(6, * ) 'MISMATCH IN NUMBER OF WAVELENGTH INTERVALS'
      call errorHandler('SNOW_TRANVAL', - 2)
    end if

    ! compute the transmissions using linear interpolation

    ! compute the interpolation weights and points once and reuse for
    ! transmission interpolation for each band.
    ! have a check to set the weights depending if the input is
    ! outside or inside the lookup table range

    do i = il1, il2
      if (c_ind(i) == 1) then
        ismu(i)  = mvidx(lsmu,       nsmu,  smu(i))
        ibc(i)   = mvidx(lbc_conc,   nbc,   bc_conc(i))
        ireff(i) = mvidx(lsnow_reff, nreff, snow_reff(i))
        iswe(i)  = mvidx(lswe,       nswe,  swe(i))

        if (smu(i) <= lsmu(1)) then
          wsmu(i,2) = 0.0
          wsmu(i,1) = 1.0 - wsmu(i,2)
        else if (smu(i) > lsmu(nsmu)) then
          wsmu(i,2) = 1.0
          wsmu(i,1) = 1.0 - wsmu(i,2)
        else
          wsmu(i,2) = (smu(i) - lsmu(ismu(i))) &
                      / (lsmu(ismu(i) + 1) - lsmu(ismu(i)))
          wsmu(i,1) = 1.0 - wsmu(i,2)
        end if

        if (bc_conc(i) <= lbc_conc(1)) then
          wbc(i,2) = 0.0
          wbc(i,1) = 1.0 - wbc(i,2)
        else if (bc_conc(i) > lbc_conc(nbc)) then
          wbc(i,2) = 1.0
          wbc(i,1) = 1.0 - wbc(i,2)
        else
          wbc(i,2) = (bc_conc(i) - lbc_conc(ibc(i))) &
                    / (lbc_conc(ibc(i) + 1) - lbc_conc(ibc(i)))
          wbc(i,1) = 1.0 - wbc(i,2)
        end if

        if (snow_reff(i) <= lsnow_reff(1)) then
          wreff(i,2) = 0.0
          wreff(i,1) = 1.0 - wreff(i,2)
        else if (snow_reff(i) > lsnow_reff(nreff)) then
          wreff(i,2) = 1.0
          wreff(i,1) = 1.0 - wreff(i,2)
        else
          wreff(i,2) = (snow_reff(i) - lsnow_reff(ireff(i))) &
                      / (lsnow_reff(ireff(i) + 1) &
                      - lsnow_reff(ireff(i)))
          wreff(i,1) = 1.0 - wreff(i,2)
        end if

        if (swe(i) <= lswe(1)) then
          wswe(i,2) = 0.0
          wswe(i,1) = 1.0 - wswe(i,2)
        else if (swe(i) > lswe(nswe)) then
          wswe(i,2) = 1.0
          wswe(i,1) = 1.0 - wswe(i,2)
        else
          wswe(i,2) = (swe(i) - lswe(iswe(i))) &
                      / (lswe(iswe(i) + 1) - lswe(iswe(i)))
          wswe(i,1) = 1.0 - wswe(i,2)
        end if
      end if
    end do ! i

    do ib = 1, nbnd
      do i = il1, il2
        if (c_ind(i) == 1) then

          isalb = mvidx(lsalb,    nsalb, salb(i,ib))

          if (salb(i,ib) <= lsalb(1)) then
            wsalb(2) = 0.0
            wsalb(1) = 1.0 - wsalb(2)
          else if (salb(i,ib) > lsalb(nsalb)) then
            wsalb(2) = 1.0
            wsalb(1) = 1.0 - wsalb(2)
          else
            wsalb(2) = (salb(i,ib) - lsalb(isalb)) &
                      / (lsalb(isalb + 1) - lsalb(isalb))
            wsalb(1) = 1.0 - wsalb(2)
          end if

          trandir(i,ib) = 0.0
          trandif (i,ib) = 0.0

          do iisalb = isalb,isalb + 1
            do iismu = ismu(i),ismu(i) + 1
              do iireff = ireff(i),ireff(i) + 1
                do iiswe = iswe(i), iswe(i) + 1
                  do iibc = ibc(i), ibc(i) + 1

                    wtt = wsmu(i,iismu - ismu(i) + 1) &
                          * wreff(i,iireff - ireff(i) + 1) &
                          * wswe(i,iiswe - iswe(i) + 1) &
                          * wbc(i,iibc - ibc(i) + 1) &
                          * wsalb(iisalb - isalb + 1)

                    trandif (i,ib) = trandif (i,ib) + wtt &
                                    * trandif_lut(iibc,iiswe,iireff,iismu,iisalb,ib)
                    trandir(i,ib) = trandir(i,ib) + wtt &
                                    * trandir_lut(iibc,iiswe,iireff,iismu,iisalb,ib)

                  end do ! iibc
                end do  ! iiswe
              end do     ! iireff
            end do        ! iismu
          end do           ! iisalb

          if (trandif (i,ib) > 1.3 .or. &
              trandif (i,ib) < 0.0) then
            write(6, * ) 'Bad trandif ',i,ib,smu(i),bc_conc(i), &
                        snow_reff(i),swe(i),salb(i,ib),trandif (i,ib)
            write(6, * ) i,ib,ismu(i),ibc(i),ireff(i),iswe(i),isalb
            call errorHandler('SNOW_TRANVAL', - 3)
          end if
          if (trandir(i,ib) > 1.3 .or. &
              trandir(i,ib) < 0.0) then
            write(6, * ) 'Bad trandir ',i,ib,smu(i),bc_conc(i), &
                        snow_reff(i),swe(i),salb(i,ib),trandir(i,ib)
            write(6, * ) i,ib,ismu(i),ibc(i),ireff(i),iswe(i),isalb
            call errorHandler('SNOW_TRANVAL', - 3)
          end if
        else
          trandif (i,ib) = - 999.0
          trandir(i,ib) = - 999.0
        end if
      end do                 ! i
    end do                    ! ib

    return
  end subroutine snow_tranval
  !! @}
  ! ------------------------------------------------------------------------------
  !> \namespace snowAlbedoTransmissMod
  !! Central module for all snow albedo and transmissivity operations
  !!
end module snowAlbedoTransmissMod
