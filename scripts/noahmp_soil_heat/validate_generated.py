"""Run complete HRLDAS/Noah-MP on exported soil-heat probe cases.

Modes: export the current harness cases, execute the official native binaries,
and score their outputs with the current registered criterion. No native
physics is replaced and the native timestep equals the hourly case timestep.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np
import pandas as pd

PROBE_ID = 'energy/soil-heat-storage-consistency'


def export_cases(destination: Path) -> None:
    from hydroturing.harness import build_case
    from hydroturing.registry import find_probe
    from hydroturing.seeds import gate_seeds

    probe = find_probe(PROBE_ID)
    destination.mkdir(parents=True, exist_ok=True)
    seeds = gate_seeds(probe.id, probe.n_seeds)
    for seed in seeds:
        case = build_case(probe, seed)
        required = {'rsds','rlds','sfcWind','huss','ps'}
        if not required.issubset(case.forcing.columns):
            raise ValueError('The current generator does not expose the atmospheric forcing')
        folder = destination / str(seed)
        folder.mkdir()
        case.forcing[[c for c in case.forcing if not c.startswith('_')]].to_csv(folder / 'forcing.csv', index=False)
        (folder / 'static.json').write_text(json.dumps(case.static,indent=2)+'\n')
        (folder / 'case.json').write_text(json.dumps({
            'probe_id':case.probe_id, 'seed':seed, 'spinup_steps':case.spinup_steps,
            'timestep':case.timestep,
            'phases':case.forcing['_phase'].tolist(),
            'criteria':[{'name':c.name,'params':dict(c.params)} for c in probe.criteria],
        },indent=2)+'\n')
    (destination / 'seeds.json').write_text(json.dumps(seeds)+'\n')
    print('Exported current harness cases:', seeds)


def namelist_setting(text: str, name: str, value: str) -> str:
    result, count = re.subn(rf'(?m)^\s*{re.escape(name)}\s*=.*$',f' {name} = {value}',text)
    if count != 1:
        raise ValueError(f'Expected one namelist setting: {name}')
    return result


def prepare_case(inputs: Path, destination: Path, native_case: Path) -> pd.DataFrame:
    forcing = pd.read_csv(inputs / 'forcing.csv')
    static = json.loads((inputs / 'static.json').read_text())
    meta = json.loads((inputs / 'case.json').read_text())
    if meta['timestep'] != 'PT1H' or len(forcing) != 96:
        raise ValueError('This reproduction expects the complete 96-hour PT1H case')
    if static['soil_porosity'] != 0.464:
        raise ValueError('The declared porosity must match native STAS category 8')
    if not np.all(forcing.pr == 0) or not np.all(forcing.huss == 0):
        raise ValueError('This reproduction requires the specified dry, rain-free case')
    destination.mkdir(parents=True)
    forcdir = destination / 'forcing'
    forcdir.mkdir()
    times = pd.to_datetime(forcing.time)
    if not np.all(np.diff(times.astype('int64')) == 3600*10**9):
        raise ValueError('Case timestamps must be consecutive hours')
    depth = static['soil_layer_depth_m']
    thickness = np.array([depth,.3,.6,1.0])
    nodes = np.cumsum(thickness)-thickness/2
    temperature = static['soil_temperature_initial']
    header = (native_case / 'forcing/bondville.dat').read_text().splitlines()[:54]
    settings = {
        'deep_soil_temperature':str(temperature),
        'skin_temperature':str(temperature),
        'soil_temperature':', '.join([str(temperature)]*4),
        'soil_moisture':', '.join([str(static['soil_water_content_initial'])]*4),
        'soil_layer_thickness':', '.join(map(str,thickness)),
        'soil_layer_nodes':', '.join(map(str,nodes)),
    }
    for index,line in enumerate(header):
        name = line.split('=',1)[0].strip()
        if name in settings:
            header[index] = f' {name} = {settings[name]}'
    header[50:54] = [
        '! Exported HydroTuring soil-heat-storage-consistency case',
        '! Native model and forcing timestep are both exactly one hour',
        '! Row i forces the integration ending at case.time[i] + 1 hour',
        '! No interpolation or internal timestep substitution by this adapter',
    ]
    rows=[]
    # HRLDAS reads the endpoint forcing for each native integration interval.
    for index in range(-1,len(forcing)):
        row = forcing.iloc[max(index,0)]
        timestamp = times.iloc[0] if index == -1 else times.iloc[index]+pd.Timedelta(hours=1)
        rows.append(f'{timestamp:%Y %m %d %H %M} {row.sfcWind:.12g} '
                    f'{row.tas+273.15:.12g} {row.huss:.12g} {row.ps:.12g} '
                    f'{row.rsds:.12g} {row.rlds:.12g} {row.pr/86400:.12g}')
    (forcdir/'bondville.dat').write_text('\n'.join(header+rows)+'\n')
    with (destination/'create_forcing.log').open('w') as log:
        subprocess.run([str(native_case/'create_point_data.exe')],cwd=forcdir,
                       stdout=log,stderr=subprocess.STDOUT,check=True)
    # The official text converter writes minute suffixes; the hourly HRLDAS
    # filename convention omits the two zero minute digits.
    for path in forcdir.glob('????????????.LDASIN_DOMAIN1'):
        path.rename(path.with_name(path.name[:10]+path.name[12:]))
    start=times.iloc[0]
    for variant in ('unmodified','diagnostic'):
        run=destination/variant
        (run/'output').mkdir(parents=True)
        shutil.copy2(native_case/variant/'hrldas.exe',run/'hrldas.exe')
        table=(native_case/variant/'NoahmpTable.TBL').read_text()
        table=namelist_setting(table,'CSOIL_DATA',str(static['soil_solid_heat_capacity']))
        (run/'NoahmpTable.TBL').write_text(table)
        nml=(native_case/variant/'namelist.hrldas').read_text()
        settings={
            'HRLDAS_SETUP_FILE':f'"{forcdir/"hrldas_setup_single_point.nc"}"',
            'INDIR':f'"{forcdir}"','OUTDIR':f'"{run/"output"}"',
            'START_YEAR':str(start.year),'START_MONTH':str(start.month),
            'START_DAY':str(start.day),'START_HOUR':str(start.hour),
            'START_MIN':str(start.minute),'KDAY':'4',
            'FORCING_TIMESTEP':'3600','NOAH_TIMESTEP':'3600',
            'SOIL_TIMESTEP':'3600','OUTPUT_TIMESTEP':'3600',
        }
        settings.update({f'soil_thick_input({i+1})':str(d) for i,d in enumerate(thickness)})
        for key,value in settings.items():
            nml=namelist_setting(nml,key,value)
        (run/'namelist.hrldas').write_text(nml)
    return forcing


def run_cases(inputs: Path, destination: Path, native_case: Path) -> None:
    from netCDF4 import Dataset

    destination.mkdir(parents=True,exist_ok=True)
    seeds=sorted(int(path.name) for path in inputs.iterdir() if path.is_dir() and path.name.isdigit())
    provenance=[]
    for seed in seeds:
        folder=destination/str(seed)
        forcing=prepare_case(inputs/str(seed),folder,native_case)
        elapsed={}
        for variant in ('unmodified','diagnostic'):
            run=folder/variant
            started=time.monotonic()
            with (run/'run.log').open('w') as log:
                subprocess.run([str(run/'hrldas.exe')],cwd=run,stdout=log,
                               stderr=subprocess.STDOUT,check=True)
            elapsed[variant]=time.monotonic()-started
        left_path=next((folder/'unmodified/output').glob('*.LDASOUT_DOMAIN1'))
        right_path=folder/'diagnostic/output'/left_path.name
        differences=[]
        with Dataset(left_path) as left,Dataset(right_path) as right:
            left.set_auto_mask(False)
            right.set_auto_mask(False)
            if set(left.variables) != set(right.variables):
                raise ValueError('Instrumentation changed native output variables')
            for name in left.variables:
                a,b=left.variables[name][:],right.variables[name][:]
                equal=np.array_equal(a,b,equal_nan=True) if a.dtype.kind in 'fc' else np.array_equal(a,b)
                if not equal:
                    differences.append(name)
            if differences or len(left.dimensions['Time']) != len(forcing)+1:
                raise ValueError(f'Native output comparison failed: {differences}')
            forcing_check={}
            for name,column in [('SWFORC','rsds'),('LWFORC','rlds')]:
                consumed=np.asarray(right.variables[name][1:]).reshape(len(forcing),-1)[:,0]
                expected=forcing[column].to_numpy()
                error=float(np.max(np.abs(consumed-expected)))
                if not np.allclose(consumed,expected,rtol=2e-7,atol=2e-5):
                    raise ValueError(f'{seed}: native {name} does not match input intervals; max error {error}')
                forcing_check[name]={'max_abs_float32_error':error,'all_intervals_match':True}
            source_record={'seed':seed,'native_dt_seconds':3600,'run_seconds':elapsed,
                           'native_variables_compared':len(left.variables),'native_records':len(left.dimensions['Time']),
                           'native_changed_variables':differences,'forcing_check':forcing_check}
        (folder/'provenance.json').write_text(json.dumps(source_record,indent=2)+'\n')
        provenance.append(source_record)
        print(json.dumps(source_record),flush=True)
    (destination/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')


def score_cases(inputs: Path, destination: Path) -> None:
    from hydroturing.criteria import get
    from hydroturing.protocol import Case,RunResult
    from hydroturing.registry import find_probe

    reports=[]
    probe=find_probe(PROBE_ID)
    seeds=sorted(int(path.name) for path in inputs.iterdir() if path.is_dir() and path.name.isdigit())
    for seed in seeds:
        source=inputs/str(seed)
        folder=destination/str(seed)
        meta=json.loads((source/'case.json').read_text())
        static=json.loads((source/'static.json').read_text())
        forcing=pd.read_csv(source/'forcing.csv')
        forcing['_phase']=meta['phases']
        case=Case(probe_id=meta['probe_id'],seed=seed,forcing=forcing,static=static,
                  spinup_steps=meta['spinup_steps'],timestep=meta['timestep'])
        thermal=pd.read_csv(folder/'diagnostic/native_thermal_steps.csv')
        column=pd.read_csv(folder/'diagnostic/native_column_steps.csv')
        if len(thermal)!=len(forcing) or len(column)!=len(forcing):
            raise ValueError('Expected exactly one native thermal step per case input row')
        for data in (thermal,column):
            if not np.isfinite(data.to_numpy(dtype=float)).all():
                raise ValueError('Non-finite native diagnostics')
            if not np.array_equal(data.step,np.arange(1,len(forcing)+1)) or not np.all(data.dt_s==3600):
                raise ValueError('Native steps do not match the one-hour case')
        if not np.all(column.soil_step_executed==1):
            raise ValueError('Soil integration did not execute at every case interval')
        capacity=static['soil_heat_capacity_areal']
        native_capacity=thermal.heat_capacity_area_j_m2_k.to_numpy()
        delta=column.t1_final_k.to_numpy()-column.t1_start_k.to_numpy()
        capacity_error=(capacity-native_capacity)*delta/3600
        step_contract={
            'max_abs_layer_depth_difference_m':float(np.max(abs(-thermal.z1_m-static['soil_layer_depth_m']))),
            'capacity_from_case_static_j_m2_k':capacity,
            'native_capacity_range_j_m2_k':[float(native_capacity.min()),float(native_capacity.max())],
            'native_capacity_relative_range':float(np.ptp(native_capacity)/capacity),
            'max_abs_capacity_approximation_error_w_m2':float(np.max(abs(capacity_error))),
            'soil_moisture_range':[float(column.soil_moisture_1_end.min()),float(column.soil_moisture_1_end.max())],
            'max_abs_latent_heat_w_m2':float(column.heat_latent_ground_w_m2.abs().max()),
            'max_abs_precip_heat_w_m2':float(column.heat_precip_advected_w_m2.abs().max()),
            'max_abs_penetrating_radiation_w_m2':float(thermal.penetrated_sw_w_m2.abs().max()),
            'max_snow_swe_mm':float(column.snow_swe_mm.max()),
            'max_soil_ice_fraction':float(column.soil_ice_max_fraction.max()),
            'max_vegetation_fraction':float(column.veg_fraction.max()),
            'max_abs_post_thermal_temperature_change_k':float((column.t1_final_k-thermal.t1_thermal_end_k).abs().max()),
            'minimum_soil_temperature_k':float(column.tsoil_min_k.min()),
            'max_native_energy_balance_error_w_m2':float(column.energy_balance_error_w_m2.abs().max()),
        }
        base=pd.DataFrame({'time':forcing.time,'hfg':thermal.g_top_w_m2,
                           'hfg_bottom':thermal.g_bottom_w_m2,'tsoil_layer':column.t1_final_k})
        base.to_csv(folder/'native_hydroturing_output.csv',index=False)
        results={}
        for name,scale in [('native',1.),('frozen_temperature',0.),('half_temperature_change',.5)]:
            table=base.copy()
            table.tsoil_layer=static['soil_temperature_initial']+scale*(table.tsoil_layer-static['soil_temperature_initial'])
            results[name]=[asdict(get(c.name)(RunResult(case,table,{},0.),probe,dict(c.params))) for c in probe.criteria]
        phase_error={}
        for label in forcing._phase.unique():
            mask=forcing._phase.eq(label).to_numpy()
            phase_error[label]={'mean_abs_w_m2':float(np.mean(abs(capacity_error[mask]))),
                                'max_abs_w_m2':float(np.max(abs(capacity_error[mask])))}
        report={'seed':seed,'probe_id':case.probe_id,'case_timestep':case.timestep,
                'native_timestep_seconds':3600,'case_contract_checks':step_contract,
                'phase_capacity_approximation_error':phase_error,'criterion_results':results}
        (folder/'criterion_validation.json').write_text(json.dumps(report,indent=2)+'\n')
        reports.append(report)
        print(seed,{name:[r['status'] for r in values] for name,values in results.items()},flush=True)
    (destination/'criterion_validation.json').write_text(json.dumps(reports,indent=2)+'\n')
    unexpected=[]
    for report in reports:
        for name,results in report['criterion_results'].items():
            expected='pass' if name=='native' else 'fail'
            if any(result['status']!=expected for result in results):
                unexpected.append(f"{report['seed']} {name}: expected {expected}")
    if unexpected:
        raise ValueError('Unexpected criterion result: '+'; '.join(unexpected))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['export','run','score'])
    parser.add_argument('--cases',type=Path,required=True)
    parser.add_argument('--workdir',type=Path)
    parser.add_argument('--native-case',type=Path)
    args=parser.parse_args()
    if args.mode=='export':
        export_cases(args.cases.resolve())
    elif args.mode=='run':
        if args.workdir is None or args.native_case is None:
            parser.error('run requires --workdir and --native-case')
        run_cases(args.cases.resolve(),args.workdir.resolve(),args.native_case.resolve())
    else:
        if args.workdir is None:
            parser.error('score requires --workdir')
        score_cases(args.cases.resolve(),args.workdir.resolve())


if __name__=='__main__':
    main()
