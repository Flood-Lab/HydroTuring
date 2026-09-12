program thermal_component_driver
  use Machine
  use NoahmpVarType
  use SoilThermalPropertyMod, only: SoilThermalProperty
  use SoilSnowThermalDiffusionMod, only: SoilSnowThermalDiffusion
  use SoilSnowTemperatureSolverMod, only: SoilSnowTemperatureSolver
  implicit none
  type(noahmp_type) :: model
  real(kind=kind_noahmp), allocatable :: a(:), b(:), c(:), rhs(:)
  real(kind=kind_noahmp) :: dt, initial_temperature, top_flux, bottom_flux
  real(kind=kind_noahmp) :: start_temperature, mean_bottom_flux, heat_capacity_area, bottom_hour_start
  integer :: hour, substep, steps_per_hour, ns, lower_index
  character(len=32) :: phase
  character(len=32) :: dt_argument

  ! Controlled thermal component experiment, not an HRLDAS/Noah-MP system run.
  ! Native soil moisture is held fixed, with no snow, ice, or water movement.
  ! The first native 0.1 m layer is the evaluated control volume. Deeper layers
  ! supply its physical lower-boundary interaction and are not its heat storage.
  ns = 4
  lower_index = -2
  dt = 300.0_kind_noahmp
  call get_command_argument(1,dt_argument)
  if (len_trim(dt_argument)>0) read(dt_argument,*) dt
  steps_per_hour = nint(3600.0_kind_noahmp/dt)
  if (steps_per_hour < 1) stop 'time step exceeds output interval'
  if (abs(dt*steps_per_hour-3600.0_kind_noahmp)>1.0e-6_kind_noahmp) stop 'time step must divide 3600 s'
  initial_temperature = 285.0_kind_noahmp
  model%config%domain%NumSoilLayer = ns
  model%config%domain%NumSnowLayerMax = 3
  model%config%domain%NumSnowLayerNeg = 0
  model%config%domain%SoilTimeStep = dt
  model%config%nmlist%OptSoilTemperatureBottom = 2
  model%config%nmlist%OptSnowSoilTempTime = 1
  model%forcing%TemperatureSoilBottom = initial_temperature
  model%energy%state%DepthSoilTempBotToSno = -8.0_kind_noahmp
  allocate(model%config%domain%DepthSnowSoilLayer(lower_index:ns))
  model%config%domain%DepthSnowSoilLayer = 0.0_kind_noahmp
  model%config%domain%DepthSnowSoilLayer(1:ns) = &
     [-0.1_kind_noahmp,-0.4_kind_noahmp,-1.0_kind_noahmp,-2.0_kind_noahmp]
  allocate(model%energy%state%TemperatureSoilSnow(lower_index:ns))
  allocate(model%energy%state%ThermConductSoilSnow(lower_index:ns))
  allocate(model%energy%state%HeatCapacSoilSnow(lower_index:ns))
  allocate(model%energy%flux%RadSwPenetrateGrd(lower_index:ns))
  allocate(model%water%param%SoilMoistureSat(ns))
  allocate(model%water%state%SoilMoisture(ns),model%water%state%SoilLiqWater(ns))
  allocate(model%energy%param%SoilQuartzFrac(ns))
  allocate(model%energy%state%HeatCapacVolSoil(ns),model%energy%state%ThermConductSoil(ns))
  allocate(a(lower_index:ns), b(lower_index:ns), c(lower_index:ns), rhs(lower_index:ns))
  model%water%param%SoilMoistureSat = 0.45_kind_noahmp
  model%water%state%SoilMoisture = 0.2_kind_noahmp
  model%water%state%SoilLiqWater = 0.2_kind_noahmp
  model%energy%param%SoilHeatCapacity = 2.0e6_kind_noahmp
  model%energy%param%SoilQuartzFrac = 0.4_kind_noahmp
  model%energy%state%TemperatureSoilSnow = initial_temperature
  model%energy%state%ThermConductSoilSnow = 0.0_kind_noahmp
  model%energy%state%HeatCapacSoilSnow = 0.0_kind_noahmp
  model%energy%flux%RadSwPenetrateGrd = 0.0_kind_noahmp
  call SoilThermalProperty(model)
  model%energy%state%ThermConductSoilSnow(1:ns) = model%energy%state%ThermConductSoil
  model%energy%state%HeatCapacSoilSnow(1:ns) = model%energy%state%HeatCapacVolSoil
  heat_capacity_area = -model%config%domain%DepthSnowSoilLayer(1) * &
     model%energy%state%HeatCapacSoilSnow(1)

  write(*,'(a)') 'hour,phase,g_top_mean_w_m2,g_bottom_mean_w_m2,t_start_k,t_end_k,' // &
     'heat_capacity_area_j_m2_k,g_bottom_hour_start_w_m2,g_bottom_hour_end_w_m2'
  do hour = 0, 95
     top_flux = 0.0_kind_noahmp
     phase = 'spinup'
     if (hour >= 24 .and. hour < 36) then
        phase = 'heating'
        top_flux = 80.0_kind_noahmp
     else if (hour >= 36) then
        phase = 'recovery'
     endif
     start_temperature = model%energy%state%TemperatureSoilSnow(1)
     bottom_hour_start = 2.0_kind_noahmp * model%energy%state%ThermConductSoilSnow(1) * &
        (model%energy%state%TemperatureSoilSnow(1)-model%energy%state%TemperatureSoilSnow(2)) / &
        (-model%config%domain%DepthSnowSoilLayer(2))
     mean_bottom_flux = 0.0_kind_noahmp
     do substep = 1, steps_per_hour
        model%energy%flux%HeatGroundTotMean = top_flux
        call SoilSnowThermalDiffusion(model,a,b,c,rhs)
        call SoilSnowTemperatureSolver(model,dt,a,b,c,rhs,bottom_flux)
        mean_bottom_flux = mean_bottom_flux + bottom_flux / real(steps_per_hour,kind_noahmp)
     enddo
     write(*,'(i0,a,a,7(a,es24.16e3))') hour,',',trim(phase),',',top_flux, &
        ',',mean_bottom_flux,',',start_temperature, &
        ',',model%energy%state%TemperatureSoilSnow(1),',',heat_capacity_area, &
        ',',bottom_hour_start,',',bottom_flux
  enddo
end program thermal_component_driver
