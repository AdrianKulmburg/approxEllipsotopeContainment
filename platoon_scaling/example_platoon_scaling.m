% This script demonstrates how to compute a terminal region for the
% platoon benchmark using the containmentLinSys algorithm, for different
% numbers of vehicles.
clear;
nVehicles = 2;

% The next line does nothing, but is a reminder that for the proper timing
% of this code, one needs to change line 29 in
% cora/global/macros/CHECKS_ENABLED.m
% to false
CHECKS_ENABLED = false;


compute_zonoLinSys = true; % Deactivate this if you don't want to compute the results from Gruber et al.
compute_volume = true; % Deactivate this if you don't want to compute the volume of the safe sets (recommended for higher dimensions)

yalmip('clear');
rng(1234);
N = (nVehicles+1)*10;
Nsimulation = 100;

%% Params
% init
Param = param_platoon_N(nVehicles);

benchmark = ['platoon_' num2str(nVehicles)];


%% Options
% general options
Opts.N = N;

Opts.timeStep = 0.1;

Opts.taylorOrder = 5;

algorithm = 'containmentLinSys';

% Testing 'feedback' methods
Opts.controlMethod = 'feedback';

Opts.terminalRegionType = 'ellipsoid';
Opts.genMethod = 'provided';
Opts.G = speye(dim(Param.X));
t_ell_feedback = tic;
T_ell_feedback = computeTerminalRegion(benchmark, algorithm, Param, Opts);
t_ell_feedback = toc(t_ell_feedback);
disp("Time elapsed for containmentLinSys with ellipsoids using the feedback method: " + num2str(t_ell_feedback))

yalmip('clear')
Opts.terminalRegionType = 'zonotope';
Opts.genMethod = 'spherical';
Opts.nGenerators = nVehicles*5;
t_zono_feedback = tic;
T_zono_feedback = computeTerminalRegion(benchmark, algorithm, Param, Opts);
t_zono_feedback = toc(t_zono_feedback);
disp("Time elapsed for containmentLinSys with zonotopes using the feedback method: " + num2str(t_zono_feedback))


if compute_zonoLinSys
    % Comparison with zonoLinSys
    yalmip('clear')
    clear Opts
    Opts.timeStep = 0.1;                            % time step size
    Opts.N = N;                                     % number of time steps
    Opts.costFun = 'geomean';                       % cost fun. opt. problem
    algorithm = 'zonoLinSys';
    
    Opts.xEq = [zeros([2 1]); repmat([10; 0], nVehicles-1, 1)]; % equilibrium point
    Opts.uEq = zeros([nVehicles 1]);                    % control input eq. point
    
    % compute terminal region
    t_zonoLinSys = tic;
    T_zonoLinSys = computeTerminalRegion(benchmark, algorithm, Param, Opts);
    t_zonoLinSys = toc(t_zonoLinSys);
    disp("Time elapsed for zonoLinSys: " + num2str(t_zonoLinSys))
end

if compute_volume
    V_ell = volume(T_ell_feedback.set)
    V_zono = volume(T_zono_feedback.set)
    if compute_zonoLinSys
        V_zonoLinSys = volume(T_zonoLinSys.set)
    end
end

yalmip('clear')
% Create simulations
t_simulations_ell_feedback = tic;
simulations_ell_feedback = T_ell_feedback.simulateRandom(Nsimulation, 'extreme');
t_simulations_ell_feedback = toc(t_simulations_ell_feedback);
t_simulations_zono_feedback = tic;
simulations_zono_feedback = T_zono_feedback.simulateRandom(Nsimulation, 'extreme');
t_simulations_zono_feedback = toc(t_simulations_zono_feedback);

if compute_zonoLinSys
    % Create simulations
    
    t_simulations_zonoLinSys = tic;
    [res,t,x,u] = T_zonoLinSys.simulateRandom(Opts.N * Opts.timeStep,Nsimulation);
    t_simulations_zonoLinSys = toc(t_simulations_zonoLinSys);
    
    % Converting to our format
    simulations_zonoLinSys = cell([1 Nsimulation]);
    for i=1:Nsimulation
        simulations_zonoLinSys{i}.x = x{i}';
        simulations_zonoLinSys{i}.y = x{i}';
        simulations_zonoLinSys{i}.yexact = x{i}';
        simulations_zonoLinSys{i}.t = t{i};
        % The rest is not necessary for our purpose
    end
end

yalmip('clear')

[averageRuntime_ell_feedback, stdErrorRuntime_ell_feedback] = re_computeInputs(T_ell_feedback, simulations_ell_feedback);
[averageRuntime_zono_feedback, stdErrorRuntime_zono_feedback] = re_computeInputs(T_zono_feedback, simulations_zono_feedback);
if compute_zonoLinSys
    [averageRuntime_zonoLinSys, stdErrorRuntime_zonoLinSys] = re_computeInputs_zonoLinSys(T_zonoLinSys, simulations_zonoLinSys,Opts.timeStep,Opts.N);
end

save("platoon_"+num2str(nVehicles)+".mat")