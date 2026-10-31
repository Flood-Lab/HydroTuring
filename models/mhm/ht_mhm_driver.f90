!> HydroTuring's driver around mHM's own library.
!!
!! This is upstream's src/mHM/mhm_driver.f90 with mhm_eval's loop written out, in
!! mhm_eval's order and through the same mhm_interface_run_* calls, and one thing
!! added: after each time step it reads the states and fluxes out of upstream's
!! global arrays, and every `steps_per_row` steps it writes one record of them over
!! the domain to `ht_rows.bin`. Upstream's netCDF writer can only give a day's
!! states as the day's mean, or an hourly table too large for a run's tmpfs, and
!! the contract needs a store's value at the end of the row. Nothing is computed
!! here that upstream does not compute, apart from weighting cells by their area.
!!
!! The steps per row come from the environment variable HT_STEPS_PER_ROW. Arguments
!! are upstream's own, parsed by upstream's parser.
program ht_mhm_driver
  use mo_kind, only: dp, i4
  use mo_common_mpi_tools, only: mpi_tools_init, mpi_tools_finalize
  use mo_mhm_cli, only: parse_command_line
  use mo_mhm_interface, only: mhm_interface_init, mhm_interface_finalize
  use mo_mhm_interface_run, only: mhm_interface_run_prepare, mhm_interface_run_prepare_domain, &
    mhm_interface_run_do_time_step, mhm_interface_run_write_output, mhm_interface_run_finished, &
    mhm_interface_run_finalize_domain, mhm_interface_run_finalize
  use mo_common_variables, only: level1, global_parameters
  use mo_global_variables, only: L1_inter, L1_snowPack, L1_soilMoist, L1_sealSTW, L1_unsatSTW, L1_satSTW, &
    L1_aETSoil, L1_aETCanopy, L1_aETSealed, L1_total_runoff, L1_preEffect, L1_melt, L1_pet_calc
  use mo_mpr_global_variables, only: L1_fSealed, L1_soilMoistSat, L1_karstLoss, L1_fAsp
  use mo_mrm_global_variables, only: mRM_runoff, L11_TSrout

  implicit none

  integer, parameter :: n_fields = 11
  integer(i4) :: steps_per_row, step, s, nrow, unit_rows, unit_q, unit_meta, h, ios
  logical :: finished
  character(len = 32) :: text
  real(dp) :: record(n_fields), accumulated(6)
  real(dp), allocatable :: w(:), fs(:), fn(:), outlet(:)

  call get_environment_variable("HT_STEPS_PER_ROW", text)
  read(text, *, iostat = ios) steps_per_row
  if (ios /= 0 .or. steps_per_row < 1) steps_per_row = 1

  call mpi_tools_init()
  call parse_command_line()
  call mhm_interface_init()

  call mhm_interface_run_prepare(global_parameters(:, 3))
  call mhm_interface_run_prepare_domain(1)

  w = level1(1)%CellArea / sum(level1(1)%CellArea)
  fs = L1_fSealed(:, 1, 1)
  fn = 1.0_dp - fs

  open(newunit = unit_rows, file = "ht_rows.bin", access = "stream", form = "unformatted", status = "replace")

  accumulated = 0.0_dp
  step = 0
  finished = .false.
  do while (.not. finished)
    call mhm_interface_run_do_time_step()
    call mhm_interface_run_write_output()
    call mhm_interface_run_finished(finished)
    step = step + 1

    ! Fluxes as mHM_updateDataset composes them, in mm per step, summed over the row.
    accumulated(1) = accumulated(1) + sum(w * (sum(L1_aETSoil, dim = 2) * fn + L1_aETCanopy + L1_aETSealed * fs))
    accumulated(2) = accumulated(2) + sum(w * L1_total_runoff)
    accumulated(3) = accumulated(3) + sum(w * L1_preEffect)
    accumulated(4) = accumulated(4) + sum(w * L1_melt)
    accumulated(5) = accumulated(5) + sum(w * L1_pet_calc)

    if (mod(step, steps_per_row) == 0) then
      ! Stores at the end of the row, each per unit of the share of the cell it covers.
      record(1) = sum(w * sum(L1_soilMoist, dim = 2) * fn)
      record(2) = sum(w * L1_sealSTW * fs)
      record(3) = sum(w * L1_unsatSTW * fn)
      record(4) = sum(w * L1_satSTW * fn)
      record(5) = sum(w * L1_snowPack)
      record(6) = sum(w * L1_inter)
      record(7:11) = accumulated(1:5)
      write(unit_rows) record
      accumulated = 0.0_dp
    end if
  end do
  close(unit_rows)

  call mhm_interface_run_finalize_domain()
  call mhm_interface_run_finalize()

  ! mRM's discharge at the gauge, m3/s at each step, averaged over the row. Only
  ! mhm_interface_finalize deallocates it.
  nrow = step / steps_per_row
  allocate(outlet(nrow))
  do s = 1, nrow
    outlet(s) = sum(mRM_runoff((s - 1) * steps_per_row + 1 : s * steps_per_row, 1)) / real(steps_per_row, dp)
  end do
  open(newunit = unit_q, file = "ht_outlet.bin", access = "stream", form = "unformatted", status = "replace")
  write(unit_q) outlet
  close(unit_q)

  open(newunit = unit_meta, file = "ht_meta.txt", status = "replace")
  write(unit_meta, *) "cells", size(w)
  write(unit_meta, *) "area_km2", sum(level1(1)%CellArea) * 1.0e-6_dp
  write(unit_meta, *) "sealed_share", sum(w * fs)
  do h = 1, size(L1_soilMoistSat, 2)
    write(unit_meta, *) "saturated_soil_mm_h", h, sum(w * L1_soilMoistSat(:, h, 1) * fn)
  end do
  write(unit_meta, *) "aspect_factor_mean", sum(w * L1_fAsp(:, 1, 1))
  write(unit_meta, *) "karst_loss_min", minval(L1_karstLoss(:, 1, 1))
  write(unit_meta, *) "routing_step_s", L11_TSrout(1)
  close(unit_meta)

  call mhm_interface_finalize()
  call mpi_tools_finalize()
end program ht_mhm_driver
