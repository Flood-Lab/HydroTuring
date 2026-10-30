% Run MARRMoT's m_18_simhyd_7p_3s, unmodified, on one case.
%
%   octave-cli --norc --no-history simhyd_run.m <in.bin> <out.bin> <delta_t> <seed> <theta,...> <S0,...>
%
% in.bin holds three doubles per row: P and Ep in mm/d, which is what the numeric form
% of input_climate takes, and T, which this model's equations never read. out.bin holds
% one row per step of what get_output returns: the fluxes in mm per step, the stores in
% mm, and the solver's step data.

args = argv();
in_path = args{1};
out_path = args{2};
delta_t = str2double(args{3});
seed = str2double(args{4});
theta = str2double(strsplit(args{5}, ','));
s0 = str2double(strsplit(args{6}, ','));

root = getenv('MARRMOT_ROOT');
addpath(genpath(fullfile(root, 'Functions')));
addpath(genpath(fullfile(root, 'Models')));

% The solver's fifth and later restarts draw from randn (MARRMoT_model.rerunSolver).
% Seeded from the request so paired variants of one seed draw the same numbers.
rand('twister', seed);
randn('state', seed);

fid = fopen(in_path, 'rb');
climate = fread(fid, [3, Inf], 'double')';
fclose(fid);

% The settings workflow_example_1.m gives, which it notes are the defaults.
solver_opts = struct();
solver_opts.resnorm_tolerance = 0.1;
solver_opts.resnorm_maxiter = 6;

m = feval('m_18_simhyd_7p_3s');
m.delta_t = delta_t;
m.theta = theta;
m.S0 = s0;
m.solver_opts = solver_opts;
m.input_climate = climate;

[flux_out, flux_in, stores, balance, solver] = m.get_output();

table = [flux_out.Ea(:), flux_out.Q(:), ...
         flux_in.Ei(:), flux_in.EXC(:), flux_in.INF(:), flux_in.INT(:), flux_in.REC(:), ...
         flux_in.Et(:), flux_in.GWF(:), flux_in.BAS(:), flux_in.SRUN(:), flux_in.Qt(:), ...
         stores.S1(:), stores.S2(:), stores.S3(:), ...
         solver.resnorm(:), double(solver.solver(:)), solver.iter(:)];

fid = fopen(out_path, 'wb');
fwrite(fid, table', 'double');
fclose(fid);

printf('MARRMOT_WATER_BALANCE_MM %.17g\n', balance);
