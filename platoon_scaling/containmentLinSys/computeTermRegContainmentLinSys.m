function T = computeTermRegContainmentLinSys(benchmark, Param, Opts)
% COMPUTETERMREGCONTAINMENTLINSYS - compute terminal region with
%       linear system approach, using only containment
%
% Syntax:
%       T = COMPUTETERMREGCONTAINMENTLINSYS(benchmark, Param, Opts)
%
% Description:
%       This function computes a terminal region with the approach in [1].
%
% Input Arguments:
%       -benchmark:  name of the considered benchmark model (see
%                    "aroc/benchmarks/dynamics/...")
%
%       -Param:      a structure containing the benchmark parameters
%           -.U:     set of admissible control inputs
%           -.Ustart:uncertainty of the first control input
%           -.W:     set of uncertain disturbances
%           -.V:     set of measurement noise
%           -.X:     set of state constraints
%           -.Y:     set of output constraints
%           -.Yexact:set of error-free output constraints
%           -.C, -.D, -.F, -.gamma:
%                    defines the output equation
%                    y = Cx + Du + Fv + gamma
%                    (Currently, C must be invertible)
%           U, X, and Y need to be polytopes, or will be transformed to 
%           polytopes within the algorithm (i.e., intervals or zonotopes
%           can also be used here, but will be transformed to polytopes).
%           W, V, and Ustart on the other hand need to be intervals or
%           zonotopes.
%
%       - Opts: a structure containing the algorithm settings
%           -.N:                  number of time steps 
%                                 [{10} / positive integer]
%           -.timeStep:           sampling time [{0.1} / positive scalar]
%           -.terminalRegionType: set type of the terminal region
%                                 [{'zonotope'} or ellipsoid]
%           -.taylorOrder:        Taylor order for certain approximations
%                                 needed during the reachability analysis
%                                 [{5} / positive integer]
%           -.genMethod:          method for computing the fixed directions
%                                 [{'spherical'} or 'provided']
%           -.nGenerators:        if .genMethod='spherical', .nGenerators
%                                 should be the number of generators for
%                                 the fixed generator matrix
%                                 [{dim(X)}/ positive integer]
%           -.G:                  if .genMethod='provided', .matrix should
%                                 be the desired matrix for the fixed
%                                 generators
%           -.K:                  feedback matrix for the terminal
%                                 controller
%           -.Q:                  state weighting matrix for the LQR
%                                 approach applied to determine the
%                                 terminal controller
%                                 [{eye(nx)} / positive definite
%                                 square matrix]
%           -.R:                  input weighting matrix for the LQR
%                                 approach applied to determine the
%                                 terminal controller
%                                 [{eye(nu)} / positive definite
%                                 square matrix]
%           -.controlMethod:      Method for generating the controller,
%                                 i.e., [{'feedback'} or 'iterative']
%
%
% Output Arguments:
%       -T:     object of class termRegContainmentLinSys
%
% See Also:
%       terminalRegion, termRegLinSysApproach
%
% References:
%       * *[1] "Approximability of the Containment Problem for Zonotopes
%              and Ellipsotopes", A. Kulmburg et al., submited to TAC
%              in 2024
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch aroc">AROC</a>, a Toolbox for Automatic Reachset-
% Optimal Controller Synthesis developed at the Chair of Robotics,
% Artificial Intelligence and Embedded Systems,
% Technische Universitaet Muenchen.
%
% For updates and further information please visit <a href="http://aroc.in.tum.de">aroc.in.tum.de</a>
%
% More Toolbox Info by searching <a href="matlab:docsearch aroc">AROC</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Adrian Kulmburg
% Website:      <a href="http://aroc.in.tum.de">aroc.in.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Copyright (c) 2023 Chair of Robotics, Artificial Intelligence and
%               Embedded Systems, TU Muenchen
%------------------------------------------------------------------

%% Emptying Yalmip variables, since this seems to change quite a lot
yalmip('clear')

%% Pre-Processing
tic;
optsInternal = cLS_setOptsInternal(benchmark,Opts,Param);
Opts = optsInternal.Opts;
Param = optsInternal.Param;
t_preprocessing = toc;



%% Setting up reachability
tic;

Rtp = cell([1 Opts.N+1]);
Rti = cell([1 Opts.N]);
ctrl_variable = cell([1 Opts.N]);
utotal = cell([1 Opts.N+1]);
additional_variables = cell([1 Opts.N]);

s = sdpvar(optsInternal.size_s, 1, 'full');
cost = -geomean(s');

constraints = [s >= sparse(optsInternal.size_s, 1)];

X0 = optsInternal.preH * diag(s) * optsInternal.postH;
c = sdpvar(optsInternal.stdDynDims.nx, 1, 'full');

% Preparing everything for the initial step
Rtp{1}.X.G = X0;
Rtp{1}.X.c = c;
% We now set up Yexact and Y, knowing that the input acts
% 'independently' in the first step
utotal{1}.G = [sparse(optsInternal.stdDynDims.nu, size(optsInternal.postH,2)) generators(optsInternal.stdConstraints.Ustart)];
utotal{1}.c = center(optsInternal.stdConstraints.Ustart);

Rtp{1}.Yexact.G = [optsInternal.stdSys.C*X0 sparse(optsInternal.stdDynDims.ny, size(generators(optsInternal.stdConstraints.Ustart),2))] + optsInternal.stdSys.D*utotal{1}.G;
Rtp{1}.Yexact.c = optsInternal.stdSys.C*c + optsInternal.stdSys.D*utotal{1}.c + optsInternal.stdSys.gamma;
Rtp{1}.Y.G = [Rtp{1}.Yexact.G optsInternal.stdSys.F];
Rtp{1}.Y.c = Rtp{1}.Yexact.c;
Rtp{1}.YnoInput.G = [optsInternal.stdSys.C*X0 sparse(optsInternal.stdDynDims.ny, size(generators(optsInternal.stdConstraints.Ustart),2)) optsInternal.stdSys.F];
Rtp{1}.YnoInput.c = optsInternal.stdSys.C*c + optsInternal.stdSys.gamma;
Rtp{1}.YnoInputNoGamma.G = [optsInternal.stdSys.C*X0 sparse(optsInternal.stdDynDims.ny, size(generators(optsInternal.stdConstraints.Ustart),2)) optsInternal.stdSys.F];
Rtp{1}.YnoInputNoGamma.c = optsInternal.stdSys.C*c;

for k=1:Opts.N
    [Rtp_k, Rti_k, ctrl_variable_k, utotal_k, additional_variables_k, additional_constraints_k] = cLS_oneStepReachability(Rtp{k}, utotal{k}, optsInternal);
    Rtp{k+1} = Rtp_k;
    Rti{k} = Rti_k;
    ctrl_variable{k} = ctrl_variable_k;
    utotal{k+1} = utotal_k;
    additional_variables{k} = additional_variables_k;

    constraints = [constraints additional_constraints_k];
end
    
t_reachability = toc;

%% Setting up state, input, and containment constraints

tic;
[variables_sic, constraints_sic] = cLS_stateInputContainment(optsInternal, Opts.terminalRegionType, Rtp, Rti, utotal, optsInternal.preH, s);
constraints = [constraints constraints_sic];
t_state_input_containment_constraints = toc;

%% Solving the optimization problems

tic;

yalmipDiagnostics = optimize(constraints,cost,optsInternal.yalmipOptions);

if yalmipDiagnostics.problem ~= 0
    warning("The solver had problems computing the result of the optimization problem. It is thus very likely that the result will be false. Here is some more detailed information:"+sprintf('\n')+ yalmipDiagnostics.info)
end

% Evaluate all results
s = value(s);
c = value(c);
for i=1:(Opts.N+1)
    Rtp{i}.X.G = value(Rtp{i}.X.G);
    Rtp{i}.X.c = value(Rtp{i}.X.c);
    Rtp{i}.Yexact.G = value(Rtp{i}.Yexact.G);
    Rtp{i}.Yexact.c = value(Rtp{i}.Yexact.c);
    Rtp{i}.Y.G = value(Rtp{i}.Y.G);
    Rtp{i}.Y.c = value(Rtp{i}.Y.c);
    Rtp{i}.YnoInput.G = value(Rtp{i}.YnoInput.G);
    Rtp{i}.YnoInput.c = value(Rtp{i}.YnoInput.c);
    Rtp{i}.YnoInputNoGamma.G = value(Rtp{i}.YnoInputNoGamma.G);
    Rtp{i}.YnoInputNoGamma.c = value(Rtp{i}.YnoInputNoGamma.c);
end

for i=1:Opts.N
    Rti{i}.X.G = value(Rti{i}.X.G);
    Rti{i}.X.c = value(Rti{i}.X.c);
    Rti{i}.Yexact.G = value(Rti{i}.Yexact.G);
    Rti{i}.Yexact.c = value(Rti{i}.Yexact.c);
    Rti{i}.Y.G = value(Rti{i}.Y.G);
    Rti{i}.Y.c = value(Rti{i}.Y.c);

    ctrl_variable{i}.U_alpha = value(ctrl_variable{i}.U_alpha);
    ctrl_variable{i}.U_beta = value(ctrl_variable{i}.U_beta);
    ctrl_variable{i}.u_c = value(ctrl_variable{i}.u_c);
end


R = zonotope(Rtp{end}.Y.c, Rtp{end}.Y.G);

if strcmp(Opts.terminalRegionType, 'ellipsoid')
    RCI = (optsInternal.stdSys.C * optsInternal.preH * diag(s)) * ellipsoid(speye(optsInternal.size_s)) + optsInternal.stdSys.C*c + optsInternal.stdSys.gamma;
    RCI_zono = zonotope(optsInternal.stdSys.C*c + optsInternal.stdSys.gamma, optsInternal.stdSys.C * optsInternal.preH * diag(s) * optsInternal.postH);

    RCIState = (optsInternal.preH * diag(s)) * ellipsoid(speye(optsInternal.size_s)) + c;
    RCIState_zono = zonotope(c, optsInternal.preH * diag(s) * optsInternal.postH);
else
    RCI = zonotope(optsInternal.stdSys.C*c + optsInternal.stdSys.gamma, optsInternal.stdSys.C * optsInternal.preH * diag(s) * optsInternal.postH);
    RCI_zono = RCI;

    RCIState = zonotope(c, optsInternal.preH * diag(s) * optsInternal.postH);
    RCIState_zono = RCIState;
end

% We transform all our results back to the original configuration
[Rtp, Rti, ctrl_variable, R, RCI, RCI_zono, RCIState, RCIState_zono] = cLS_inverseStdDynSys(Rtp, Rti, ctrl_variable, R, RCI, RCI_zono, RCIState, RCIState_zono, optsInternal);

t_solver = toc;

%% Construct the terminal region object

debug.s = s;
debug.cost = cost;
debug.Rtp = Rtp;
debug.Rti = Rti;
debug.ctrl_variable = ctrl_variable;
debug.utotal = utotal;
debug.additional_variables = additional_variables;
debug.constraints = constraints;
debug.optsInternal = optsInternal;
debug.yalmipDiagnostics = yalmipDiagnostics;


runtimes.preprocessing = t_preprocessing;
runtimes.reachability = t_reachability;
runtimes.state_input_containment_constraints = t_state_input_containment_constraints;
runtimes.solver = t_solver;
runtimes.total = t_preprocessing + t_reachability + t_state_input_containment_constraints + t_solver;

T = termRegContainmentLinSys(optsInternal.dynamics, RCI, RCI_zono, RCIState, RCIState_zono, Rtp, R, s, ctrl_variable, debug, runtimes, optsInternal);

end