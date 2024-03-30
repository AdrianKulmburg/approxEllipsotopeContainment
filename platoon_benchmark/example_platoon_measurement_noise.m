% This script demonstrates how to compute a terminal region for the
% platoon benchmark using the containmentLinSys algorithm, with measurement
% noise

clear;
yalmip('clear');
rng(1234);
N = 50;
Nsimulation = 10;
tol_verification = 1e-3;
nVehicles = 2;

%% Params
% init
Param = param_platoon_N(nVehicles);

benchmark = ['platoon_' num2str(nVehicles)];

% Setting up custom E, F, and V
% Let's start with E
Param.D = [];
for i=1:nVehicles
    e2i = sparse(2*nVehicles,1);
    e2i(2*i) = 1;
    Param.D = [Param.D e2i];
end
% We continue with F; for our case, we just choose F = E
Param.D = 0.01*Param.D;
Param.F = eye(2*nVehicles);
% And finally, we choose a simple V
Param.V = 0.01*interval(-ones([2*nVehicles 1]), ones([2*nVehicles 1]));



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

if volume(T_ell_feedback.set) < 1e-5
    error("Too much")
end

yalmip('clear')
Opts.terminalRegionType = 'zonotope';
Opts.genMethod = 'spherical';
Opts.nGenerators = 20;
t_zono_feedback = tic;
T_zono_feedback = computeTerminalRegion(benchmark, algorithm, Param, Opts);
t_zono_feedback = toc(t_zono_feedback);
disp("Time elapsed for containmentLinSys with zonotopes using the feedback method: " + num2str(t_zono_feedback))

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

yalmip('clear')

averageRuntime_ell_feedback = re_computeInputs(T_ell_feedback, simulations_ell_feedback);
averageRuntime_zono_feedback = re_computeInputs(T_zono_feedback, simulations_zono_feedback);


%% Creating Plot

% Color palette for people with colorblindness. See
% T. B. Plante, M. Cushman, "Choosing color palettes for scientific
% figures", 2020
RPTH_blue = [0, 92, 171]./255;
RPTH_red = [227, 27, 35]./255;
RPTH_yellow = [255, 195, 37]./255;

figure;
sgtitle(sprintf("Safe terminal region for platoon benchmark with 2 vehicles\nand measurement noise"));

subplot(1,2,1)
hold on
title("Ellipsoid approach")

pell_feedback = plot(T_ell_feedback.set,[1 2],'Color', RPTH_blue);
pzono_feedback = plot(T_zono_feedback.set,[1 2],'--','Color', RPTH_red);

for i = 1:Nsimulation
    y = simulations_ell_feedback{i}.y;
    y1 = y(1,:);
    y2 = y(2,:);
    plot(y1, y2, 'k')
end
axis square
xlabel("$x_1$", 'Interpreter', 'latex')
ylabel("$x_2$", 'Interpreter', 'latex')

subplot(1,2,2)
hold on
title("Zonotope approach")

pell_feedback = plot(T_ell_feedback.set,[1 2],'--','Color', RPTH_blue);
pzono_feedback = plot(T_zono_feedback.set,[1 2],'Color', RPTH_red);

for i = 1:Nsimulation
    y = simulations_zono_feedback{i}.y;
    y1 = y(1,:);
    y2 = y(2,:);
    plot(y1, y2, 'k')
end

axis square
xlabel("$x_1$", 'Interpreter', 'latex')
ylabel("$x_2$", 'Interpreter', 'latex')

save platoon_measurement_noise.mat
matlab2tikz('platoon_measurement_noise.tex')