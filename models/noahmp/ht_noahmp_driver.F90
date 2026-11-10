! A single-column driver for Noah-MP v5.2.1 under HydroTuring.
!
! It does what HRLDAS's hrldas.exe does for one land point, with the netCDF I/O
! taken out: land_driver_ini's set-up calls (NoahmpReadNamelist,
! NoahmpIOVarInitDefault, NoahmpReadTable, NoahmpInitMain), then for every row
! land_driver_exe's forcing assignment and the body of NoahmpDriverMain, whose
! calls are repeated here unchanged so that the column type is still in scope
! when NoahmpMain returns. That is the one reason this file exists instead of a
! call to NoahmpDriverMain: the conductivity behind hfg_bottom and the budget
! terms HRLDAS does not copy out are readable only there. Nothing in upstream's
! source is edited. verify_driver.py compares this program with hrldas.exe.
!
! Input, in the working directory: namelist.hrldas (upstream's own example with
! the step and layers set), NoahmpTable.TBL (upstream's, unedited), control.nml
! (below) and forcing.bin, one record per Noah-MP step of FORCING_COLUMNS native
! doubles. Output: output.bin, one record per row of OUTPUT_COLUMNS native doubles.
! A row is steps_per_row Noah-MP steps. Over them a flux is averaged, and a state
! and an instantaneous value are the last step's, in hold_row. At one step per row
! the record is the step's own values.
!
! A negative value in the overrides leaves upstream's table alone.

program ht_noahmp_driver

  use Machine
  use NoahmpVarType
  use NoahmpIOVarType
  use NoahmpIOVarInitMod
  use NoahmpReadNamelistMod
  use NoahmpReadTableMod
  use NoahmpInitMainMod
  use ConfigVarInitMod
  use EnergyVarInitMod
  use ForcingVarInitMod
  use WaterVarInitMod
  use BiochemVarInitMod
  use ConfigVarInTransferMod
  use EnergyVarInTransferMod
  use ForcingVarInTransferMod
  use WaterVarInTransferMod
  use BiochemVarInTransferMod
  use ConfigVarOutTransferMod
  use ForcingVarOutTransferMod
  use EnergyVarOutTransferMod
  use WaterVarOutTransferMod
  use BiochemVarOutTransferMod
  use NoahmpMainMod

  implicit none

  integer, parameter :: FORCING_COLUMNS = 10
  integer, parameter :: OUTPUT_COLUMNS = 49

  type(NoahmpIO_type) :: NoahmpIO

  integer :: nrows = 0, steps_per_row = 1
  integer :: vegtype = 12, soiltype = 8
  integer :: iswater = 17, isurban = 13, isice = 15, islake = 21
  real(8) :: latitude = 40.01d0, longitude = -88.37d0
  real(8) :: tmn = 285.d0, skin = 285.d0
  real(8) :: gvfmax = 96.d0, gvfmin = 1.d0, lai_initial = 2.d0
  real(8) :: soil_t(4) = 285.d0, soil_m(4) = 0.2d0
  real(8) :: swe = 0.d0, snow_depth = 0.d0
  real(8) :: csoil_override = -1.d0, eg_override = -1.d0, ch2op_override = -1.d0
  namelist /control/ nrows, steps_per_row, vegtype, soiltype, iswater, isurban, isice, islake, &
       latitude, longitude, tmn, skin, gvfmax, gvfmin, lai_initial, soil_t, soil_m, &
       swe, snow_depth, csoil_override, eg_override, ch2op_override

  real(8) :: row(FORCING_COLUMNS), out(OUTPUT_COLUMNS), held(OUTPUT_COLUMNS)
  integer :: unit_control, unit_forcing, unit_output, step, k, irow, isub

  open(newunit=unit_control, file='control.nml', status='old', action='read')
  read(unit_control, nml=control)
  close(unit_control)

  ! land_driver_ini, in its order. HRLDAS reads the domain from a setup file;
  ! here it is one point.
  call NoahmpReadNamelist(NoahmpIO)
  NoahmpIO%xstart = 1; NoahmpIO%xend = 1; NoahmpIO%ystart = 1; NoahmpIO%yend = 1
  NoahmpIO%ix = 1; NoahmpIO%jx = 1
  NoahmpIO%iswater = iswater; NoahmpIO%islake = islake
  NoahmpIO%isurban = isurban; NoahmpIO%isice = isice
  NoahmpIO%llanduse = 'MODIFIED_IGBP_MODIS_NOAH'
  NoahmpIO%dx = 1000.0; NoahmpIO%dy = 1000.0
  NoahmpIO%ids = 1; NoahmpIO%ide = 1; NoahmpIO%jds = 1; NoahmpIO%jde = 1
  NoahmpIO%kds = 1; NoahmpIO%kde = 2
  NoahmpIO%its = 1; NoahmpIO%ite = 1; NoahmpIO%jts = 1; NoahmpIO%jte = 1
  NoahmpIO%kts = 1; NoahmpIO%kte = 2
  NoahmpIO%ims = 1; NoahmpIO%ime = 1; NoahmpIO%jms = 1; NoahmpIO%jme = 1
  NoahmpIO%kms = 1; NoahmpIO%kme = 2
  NoahmpIO%ixfull = 1; NoahmpIO%jxfull = 1; NoahmpIO%ixpar = 1; NoahmpIO%jxpar = 1
  NoahmpIO%xstartpar = 1; NoahmpIO%ystartpar = 1

  call NoahmpIOVarInitDefault(NoahmpIO)
  call NoahmpReadTable(NoahmpIO)

  ! The case's capacities, emissivity and heat capacity, in the arrays
  ! NoahmpReadTable filled. The table file itself is not edited.
  if (csoil_override >= 0.d0) NoahmpIO%CSOIL_TABLE = real(csoil_override, kind_noahmp)
  if (eg_override >= 0.d0) NoahmpIO%EG_TABLE(1) = real(eg_override, kind_noahmp)
  if (ch2op_override >= 0.d0) NoahmpIO%CH2OP_TABLE(vegtype) = real(ch2op_override, kind_noahmp)

  NoahmpIO%IVGTYP(1,1) = vegtype
  NoahmpIO%ISLTYP(1,1) = soiltype
  NoahmpIO%TERRAIN(1,1) = 0.0
  NoahmpIO%TMN(1,1) = real(tmn, kind_noahmp)
  NoahmpIO%XLAT(1,1) = real(latitude, kind_noahmp)
  NoahmpIO%XLONG(1,1) = real(longitude, kind_noahmp)
  NoahmpIO%XLAND(1,1) = 1.0
  NoahmpIO%SEAICE(1,1) = 0.0
  NoahmpIO%XICE(1,1) = 0.0
  NoahmpIO%MSFTX(1,1) = 1.0
  NoahmpIO%MSFTY(1,1) = 1.0
  NoahmpIO%CROPTYPE = 0
  NoahmpIO%TD_FRACTION = 0.0
  NoahmpIO%SLOPETYP = 1
  NoahmpIO%DZS = NoahmpIO%SOIL_THICK_INPUT(1:NoahmpIO%NSOIL)
  NoahmpIO%ITIMESTEP = 1
  NoahmpIO%restart_flag = .false.

  ! READINIT_HRLDAS and READVEG_HRLDAS, without the file.
  do k = 1, NoahmpIO%NSOIL
     NoahmpIO%SMOIS(1,k,1) = real(soil_m(k), kind_noahmp)
     NoahmpIO%TSLB(1,k,1) = real(soil_t(k), kind_noahmp)
  enddo
  NoahmpIO%CANWAT(1,1) = 0.0
  NoahmpIO%TSK(1,1) = real(skin, kind_noahmp)
  NoahmpIO%SNOW(1,1) = real(swe, kind_noahmp)
  NoahmpIO%SNOWH(1,1) = real(snow_depth, kind_noahmp)
  NoahmpIO%FNDSNOWH = .true.
  NoahmpIO%VEGFRA(1,1) = real(gvfmax, kind_noahmp)
  NoahmpIO%LAI(1,1) = real(lai_initial, kind_noahmp)
  NoahmpIO%GVFMIN(1,1) = real(gvfmin, kind_noahmp)
  NoahmpIO%GVFMAX(1,1) = real(gvfmax, kind_noahmp)

  call NoahmpInitMain(NoahmpIO)

  open(newunit=unit_forcing, file='forcing.bin', status='old', action='read', &
       access='stream', form='unformatted')
  open(newunit=unit_output, file='output.bin', status='replace', action='write', &
       access='stream', form='unformatted')

  step = 0
  do irow = 1, nrows
     do isub = 1, steps_per_row
        step = step + 1
        read(unit_forcing) row
        call assign_forcing(row, step)
        call column_step(NoahmpIO, out)
        call hold_row(out, isub, held)
     enddo
     write(unit_output) held
  enddo

  close(unit_forcing)
  close(unit_output)

contains

  ! What a row reports when it took several Noah-MP steps. The columns are those
  ! of read_column.
  subroutine hold_row(out, isub, held)
    real(8), intent(in) :: out(OUTPUT_COLUMNS)
    integer, intent(in) :: isub
    real(8), intent(inout) :: held(OUTPUT_COLUMNS)
    integer :: i

    if (isub == 1) then
       held = out
       return
    endif
    do i = 1, OUTPUT_COLUMNS
       select case (i)
       case (6, 7, 8, 10, 11, 37:42, 44:46, 48)
          ! End of the row: the stores and temperatures, and the emitted longwave
          ! and radiative temperature, which are one instant's identity.
          held(i) = out(i)
       case (20, 27)
          ! Upstream's balance errors: the worst step, which a mean would hide.
          if (abs(out(i)) > abs(held(i))) held(i) = out(i)
       case (43)
          held(i) = min(held(i), out(i))
       case default
          held(i) = held(i) + (out(i) - held(i)) / isub
       end select
    enddo
  end subroutine hold_row

  ! land_driver_exe between READFORC_HRLDAS and the call to NoahmpDriverMain.
  ! Columns: SWDOWN, LWDOWN, T2D (K), wind, specific humidity, pressure (Pa),
  ! rain (mm/s), cosine of the zenith angle, JULIAN, year.
  subroutine assign_forcing(row, step)
    real(8), intent(in) :: row(FORCING_COLUMNS)
    integer, intent(in) :: step
    real(kind_noahmp) :: q

    NoahmpIO%SWDOWN(1,1) = real(row(1), kind_noahmp)
    NoahmpIO%GLW(1,1) = real(row(2), kind_noahmp)
    NoahmpIO%T_PHY(1,1,1) = real(row(3), kind_noahmp)
    NoahmpIO%U_PHY(1,1,1) = real(row(4), kind_noahmp)
    NoahmpIO%V_PHY(1,1,1) = 0.0
    q = real(row(5), kind_noahmp)
    NoahmpIO%QV_CURR(1,1,1) = q / (1.0 - q)
    NoahmpIO%P8W(1,1,1) = real(row(6), kind_noahmp)
    NoahmpIO%RAINBL(1,1) = real(row(7), kind_noahmp)
    NoahmpIO%SNOWBL(1,1) = 0.0

    NoahmpIO%P8W(:,2,:) = NoahmpIO%P8W(:,1,:)
    NoahmpIO%T_PHY(:,2,:) = NoahmpIO%T_PHY(:,1,:)
    NoahmpIO%U_PHY(:,2,:) = NoahmpIO%U_PHY(:,1,:)
    NoahmpIO%V_PHY(:,2,:) = NoahmpIO%V_PHY(:,1,:)
    NoahmpIO%QV_CURR(:,2,:) = NoahmpIO%QV_CURR(:,1,:)
    NoahmpIO%RAINBL = NoahmpIO%RAINBL * NoahmpIO%DTBL
    NoahmpIO%SNOWBL = NoahmpIO%SNOWBL * NoahmpIO%DTBL
    NoahmpIO%SR = 0.0
    NoahmpIO%RAINCV = 0.0
    NoahmpIO%RAINNCV = NoahmpIO%RAINBL
    NoahmpIO%RAINSHV = 0.0
    NoahmpIO%SNOWNCV = NoahmpIO%SNOWBL
    NoahmpIO%GRAUPELNCV = 0.0
    NoahmpIO%HAILNCV = 0.0
    NoahmpIO%DZ8W = 2 * NoahmpIO%ZLVL
    NoahmpIO%SWDDIR = NoahmpIO%SWDOWN * 0.7
    NoahmpIO%SWDDIF = NoahmpIO%SWDOWN * 0.3

    NoahmpIO%COSZEN(1,1) = real(row(8), kind_noahmp)
    NoahmpIO%JULIAN = real(row(9), kind_noahmp)
    NoahmpIO%YR = nint(row(10))

    if (step == 1) then
       NoahmpIO%EAHXY = (NoahmpIO%P8W(:,1,:) * NoahmpIO%QV_CURR(:,1,:)) / (0.622 + NoahmpIO%QV_CURR(:,1,:))
       NoahmpIO%TAHXY = NoahmpIO%T_PHY(:,1,:)
       NoahmpIO%CHXY = 0.1
       NoahmpIO%CMXY = 0.1
    endif

    NoahmpIO%MP_RAINC = NoahmpIO%RAINCV
    NoahmpIO%MP_RAINNC = NoahmpIO%RAINNCV
    NoahmpIO%MP_SHCV = NoahmpIO%RAINSHV
    NoahmpIO%MP_SNOW = NoahmpIO%SNOWNCV
    NoahmpIO%MP_GRAUP = NoahmpIO%GRAUPELNCV
    NoahmpIO%MP_HAIL = NoahmpIO%HAILNCV
    NoahmpIO%ITIMESTEP = step
  end subroutine assign_forcing

  ! NoahmpDriverMain for I = J = 1, a land point with no sea ice, no urban
  ! physics and soil processes at every step. `noahmp` is local, as it is there,
  ! so each step starts from what the transfers put into it.
  subroutine column_step(NoahmpIO, out)
    type(NoahmpIO_type), intent(inout) :: NoahmpIO
    real(8), intent(out) :: out(OUTPUT_COLUMNS)
    type(noahmp_type) :: noahmp
    integer :: K

    NoahmpIO%CALCULATE_SOIL = .false.
    NoahmpIO%SOIL_UPDATE_STEPS = nint(NoahmpIO%SOILTSTEP / NoahmpIO%DTBL)
    NoahmpIO%SOIL_UPDATE_STEPS = max(NoahmpIO%SOIL_UPDATE_STEPS, 1)
    if (NoahmpIO%SOIL_UPDATE_STEPS /= 1) stop 'this driver runs soil processes at every step'
    NoahmpIO%ACC_SSOILXY = 0.0
    NoahmpIO%ACC_QINSURXY = 0.0
    NoahmpIO%ACC_QSEVAXY = 0.0
    NoahmpIO%ACC_ETRANIXY = 0.0
    NoahmpIO%ACC_DWATERXY = 0.0
    NoahmpIO%ACC_PRCPXY = 0.0
    NoahmpIO%ACC_ECANXY = 0.0
    NoahmpIO%ACC_ETRANXY = 0.0
    NoahmpIO%ACC_EDIRXY = 0.0
    NoahmpIO%ACC_GLAFLWXY = 0.0
    NoahmpIO%CALCULATE_SOIL = mod(NoahmpIO%ITIMESTEP, NoahmpIO%SOIL_UPDATE_STEPS) == 0

    NoahmpIO%YEARLEN = 365
    if (mod(NoahmpIO%YR, 4) == 0) then
       NoahmpIO%YEARLEN = 366
       if (mod(NoahmpIO%YR, 100) == 0) then
          NoahmpIO%YEARLEN = 365
          if (mod(NoahmpIO%YR, 400) == 0) NoahmpIO%YEARLEN = 366
       endif
    endif

    NoahmpIO%ZSOIL(1) = -NoahmpIO%DZS(1)
    do K = 2, NoahmpIO%NSOIL
       NoahmpIO%ZSOIL(K) = -NoahmpIO%DZS(K) + NoahmpIO%ZSOIL(K-1)
    enddo

    NoahmpIO%J = 1
    NoahmpIO%I = 1

    call ConfigVarInitDefault(noahmp)
    call ConfigVarInTransfer(noahmp, NoahmpIO)
    call ForcingVarInitDefault(noahmp)
    call ForcingVarInTransfer(noahmp, NoahmpIO)
    call EnergyVarInitDefault(noahmp)
    call EnergyVarInTransfer(noahmp, NoahmpIO)
    call WaterVarInitDefault(noahmp)
    call WaterVarInTransfer(noahmp, NoahmpIO)
    call BiochemVarInitDefault(noahmp)
    call BiochemVarInTransfer(noahmp, NoahmpIO)

    if (noahmp%config%domain%VegType == noahmp%config%domain%IndexIcePoint) &
         stop 'this driver does not run a glacier point'
    noahmp%config%domain%IndicatorIceSfc = 0
    call NoahmpMain(noahmp)

    call read_column(noahmp, out)

    call ConfigVarOutTransfer(noahmp, NoahmpIO)
    call ForcingVarOutTransfer(noahmp, NoahmpIO)
    call EnergyVarOutTransfer(noahmp, NoahmpIO)
    call WaterVarOutTransfer(noahmp, NoahmpIO)
    call BiochemVarOutTransfer(noahmp, NoahmpIO)
  end subroutine column_step

  ! The column as NoahmpMain left it, in the order ht_adapter.OUTPUT_NAMES names.
  ! Read before the out-transfers, which multiply nothing on a soil point.
  subroutine read_column(noahmp, out)
    type(noahmp_type), intent(in) :: noahmp
    real(8), intent(out) :: out(OUTPUT_COLUMNS)
    integer :: k
    real(8) :: soil_water, flux_bottom

    soil_water = 0.d0
    do k = 1, noahmp%config%domain%NumSoilLayer
       soil_water = soil_water + dble(noahmp%water%state%SoilMoisture(k)) * &
            dble(noahmp%config%domain%ThicknessSnowSoilLayer(k)) * 1000.d0
    enddo
    ! The Fourier flux of SoilSnowThermalDiffusion through the base of the first
    ! soil layer, with the conductivity the solver used and the temperatures it
    ! returned, so the same flux its implicit equation balanced. Not a residual.
    flux_bottom = 2.d0 * dble(noahmp%energy%state%ThermConductSoilSnow(1)) * &
         (dble(noahmp%energy%state%TemperatureSoilSnow(1)) - dble(noahmp%energy%state%TemperatureSoilSnow(2))) / &
         (-dble(noahmp%config%domain%DepthSoilLayer(2)))

    out = (/ &
         dble(noahmp%forcing%RadSwDownRefHeight), &           !  1
         dble(noahmp%forcing%RadLwDownRefHeight), &           !  2
         dble(noahmp%energy%flux%RadSwAbsSfc), &              !  3
         dble(noahmp%energy%flux%RadSwReflSfc), &             !  4
         dble(noahmp%energy%flux%RadLwNetSfc), &              !  5
         dble(noahmp%energy%flux%RadLwEmitSfc), &             !  6
         dble(noahmp%energy%state%TemperatureRadSfc), &       !  7
         dble(noahmp%energy%state%EmissivitySfc), &           !  8
         dble(noahmp%energy%state%AlbedoSfc), &               !  9
         dble(noahmp%energy%state%TemperatureGrd), &          ! 10
         dble(noahmp%energy%state%TemperatureCanopy), &       ! 11
         dble(noahmp%energy%flux%HeatSensibleSfc), &          ! 12
         dble(noahmp%energy%flux%HeatLatentCanopy), &         ! 13
         dble(noahmp%energy%flux%HeatLatentGrd), &            ! 14
         dble(noahmp%energy%flux%HeatLatentTransp), &         ! 15
         dble(noahmp%energy%flux%HeatLatentIrriEvap), &       ! 16
         dble(noahmp%energy%flux%HeatGroundTot), &            ! 17
         dble(noahmp%energy%flux%HeatPrecipAdvSfc), &         ! 18
         dble(noahmp%energy%flux%HeatCanStorageChg), &        ! 19
         dble(noahmp%energy%state%EnergyBalanceError), &      ! 20
         dble(noahmp%energy%flux%RadSwPenetrateGrd(1)), &     ! 21
         dble(noahmp%water%flux%EvapCanopyNet), &             ! 22
         dble(noahmp%water%flux%Transpiration), &             ! 23
         dble(noahmp%water%flux%EvapGroundNet), &             ! 24
         dble(noahmp%water%flux%RunoffSurface), &             ! 25
         dble(noahmp%water%flux%RunoffSubsurface), &          ! 26
         dble(noahmp%water%state%WaterBalanceError), &        ! 27
         dble(noahmp%water%flux%SnowBotOutflow), &            ! 28
         dble(noahmp%water%flux%RainfallGround), &            ! 29
         dble(noahmp%water%flux%SnowfallGround), &            ! 30
         dble(noahmp%water%flux%SublimSnowSfcIce), &          ! 31
         dble(noahmp%water%flux%FrostSnowSfcIce), &           ! 32
         dble(noahmp%water%flux%SublimCanopyIce), &           ! 33
         dble(noahmp%water%state%PondSfcThinSnwMelt), &        ! 34
         dble(noahmp%water%state%PondSfcThinSnwComb), &        ! 35
         dble(noahmp%water%state%PondSfcThinSnwTrans), &       ! 36
         soil_water, &                                        ! 37
         dble(noahmp%water%state%SnowWaterEquiv), &           ! 38
         dble(noahmp%water%state%SnowDepth), &                ! 39
         dble(noahmp%water%state%CanopyLiqWater), &           ! 40
         dble(noahmp%water%state%CanopyIce), &                ! 41
         dble(noahmp%water%state%WaterStorageAquifer), &      ! 42
         dble(noahmp%config%domain%NumSnowLayerNeg), &        ! 43
         dble(noahmp%energy%state%VegFrac), &                 ! 44
         dble(noahmp%energy%state%LeafAreaIndEff + noahmp%energy%state%StemAreaIndEff), & ! 45
         dble(noahmp%energy%state%TemperatureSoilSnow(1)), &  ! 46
         flux_bottom, &                                       ! 47
         dble(noahmp%energy%state%HeatCapacSoilSnow(1)), &    ! 48
         merge(dble(noahmp%water%flux%RainfallGround), 0.d0, &
               noahmp%config%domain%NumSnowLayerNeg == 0) /)  ! 49
  end subroutine read_column

end program ht_noahmp_driver
