% This script demonstrates how to compute a terminal region for the
% platoon benchmark using the containmentLinSys algorithm

% Careful! The following code requires about 6 hours to be executed, due to
% the zonoLinSys taking a long time


clear;
yalmip('clear');
rng(1234);
N = 50;
Nsimulation = 100;
tol_verification = 1e-3;

%% Params
% init
Param = param_platoon();

benchmark = 'platoon';

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
Opts.nGenerators = 20;
t_zono_feedback = tic;
T_zono_feedback = computeTerminalRegion(benchmark, algorithm, Param, Opts);
t_zono_feedback = toc(t_zono_feedback);
disp("Time elapsed for containmentLinSys with zonotopes using the feedback method: " + num2str(t_zono_feedback))

% Comparison with zonoLinSys
yalmip('clear')
clear Opts
Opts.timeStep = 0.1;                            % time step size
Opts.N = N;                                     % number of time steps
Opts.costFun = 'geomean';                       % cost fun. opt. problem
algorithm = 'zonoLinSys';

Opts.xEq = [0; 0; 10; 0; 10; 0; 10; 0];     % equilibrium point
Opts.uEq = [0; 0; 0; 0];                    % control input eq. point

% compute terminal region
t_zonoLinSys = tic;
T_zonoLinSys = computeTerminalRegion('platoon', algorithm, Param, Opts);
t_zonoLinSys = toc(t_zonoLinSys);
disp("Time elapsed for zonoLinSys: " + num2str(t_zonoLinSys))




yalmip('clear')
% Create simulations
t_simulations_ell_feedback = tic;
simulations_ell_feedback = T_ell_feedback.simulateRandom(Nsimulation, 'extreme');
t_simulations_ell_feedback = toc(t_simulations_ell_feedback);
t_simulations_zono_feedback = tic;
simulations_zono_feedback = T_zono_feedback.simulateRandom(Nsimulation, 'extreme');
t_simulations_zono_feedback = toc(t_simulations_zono_feedback);

disp("Ellipsoid Feedback Verification...")
T_ell_feedback.verifyTrajectory(simulations_ell_feedback,tol_verification);
disp("Zonotope Feedback Verification...")
T_zono_feedback.verifyTrajectory(simulations_zono_feedback,tol_verification);

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

yalmip('clear')

averageRuntime_ell_feedback = re_computeInputs(T_ell_feedback, simulations_ell_feedback);
averageRuntime_zono_feedback = re_computeInputs(T_zono_feedback, simulations_zono_feedback);
averageRuntime_zonoLinSys = re_computeInputs_zonoLinSys(T_zonoLinSys, simulations_zonoLinSys,Opts.timeStep,Opts.N);

%% Creating Plot

% Color palette for people with colorblindness. See
% T. B. Plante, M. Cushman, "Choosing color palettes for scientific
% figures", 2020
RPTH_blue = [0, 92, 171]./255;
RPTH_red = [227, 27, 35]./255;
RPTH_yellow = [255, 195, 37]./255;

figure; hold on;
title("Terminal region for platoon with 4 vehicles");

pell_feedback = plot(T_ell_feedback.set,[1 2],'Color', RPTH_blue,'DisplayName','Ellipsoid approach');
pzono_feedback = plot(T_zono_feedback.set,[1 2],'Color',RPTH_red,'DisplayName','Zonotope approach');

pzonoLinSys = plot(T_zonoLinSys.set, [1 2], 'Color', RPTH_yellow, 'DisplayName','Approach from \cite{gruber}');

legend();
xlabel("$x_1$", 'Interpreter', 'latex')
ylabel("$x_2$", 'Interpreter', 'latex')

save large_platoon.mat
matlab2tikz('large_platoon.tex')