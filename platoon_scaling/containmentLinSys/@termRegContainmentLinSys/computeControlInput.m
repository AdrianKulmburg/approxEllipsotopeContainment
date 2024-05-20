function u = computeControlInput(T, y, currentInput, nTimeStep, varargin)
% COMPUTECONTROLINPUT - computes the control inputs one has to choose so
%                       that an output y starting in the terminal region
%                       ends up again in the terminal region
%
% Syntax:
%        u = computeControlInput(T, y, currentInput, nTimeStep)
%        u = computeControlInput(T, y, currentInput, nTimeStep, beta)
%
% Input Arguments:
%       -y:         A point that lies in the terminal region of the output
%                   (and for which the state lies in the state terminal
%                   region) if nTimeStep=1, and lies in the reachable set
%                   of the terminal region at time
%                   (nTimeStep-1) * T.Opts.timeStep
%                   (class: double)
%
%       -nTimeStep: Index of the current time step (class: double, must be
%                   a positive integer)
%
%       -beta:      Not necessary for the iterative method, but mandatory
%                   for the feedback method: Initial parametrization of the
%                   point x_0 that led to y; in other words, if y is the
%                   current trajectory that started with y0, then beta is
%                   given by
%                   y0 = G_0 * beta + c_0
%                   where G_0 and c_0 are the generator matrix and center
%                   of the initial set X0
%
% Output Arguments:
%       -u:       Inputs that have to be applied at time
%                 nTimeStep * T.Opts.timeStep, in order to continue
%                 steering the system back to the terminal region
%
% Description:
%       Computes the control inputs one has to choose so that an output y
%       starting in the terminal region ends up again in the terminal
%       region.
%
%
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
% Copyright (c) 2024 Chair of Robotics, Artificial Intelligence and
%               Embedded Systems, TU Muenchen
%------------------------------------------------------------------

global CHECKS_ENABLED
if CHECKS_ENABLED

    inputArgsCheck({{y,'att','numeric',{'vector','nonnan','finite'}}});
    equalDimCheck(T.set, y)
    inputArgsCheck({{nTimeStep,'att','numeric',{'scalar','nonnan','finite','positive','integer'}}});

    if nTimeStep > T.Opts.N
        error('nTimeStep must be a positive integer between 1 and T.Opts.N.')
    end

    if length(varargin) > 2
        error("Too many input arguments.")
    elseif length(varargin) == 1
        error("You specified alpha, but not beta. The algorithm needs both.")
    end
end


if isempty(varargin)
    error("For the feedback method, the parameters alpha and beta have to be specified.")
end
alpha = varargin{1};
%beta = varargin{2};

termReg = T.set;
errorSet = T.V;

if CHECKS_ENABLED
    inputArgsCheck({{alpha,'att','numeric',{'vector','nonnan','finite'}}});
    %inputArgsCheck({{beta,'att','numeric',{'vector','nonnan','finite'}}});
    
    if strcmp(T.Opts.terminalRegionType, 'zonotope') && size(generators(termReg),2) ~= length(alpha)
        error("Dimension mismatch between alpha and the number of generators for the (error-free) initial set.")
    elseif strcmp(T.Opts.terminalRegionType, 'ellipsoid') && dim(termReg) ~= length(alpha)
        error("Dimension mismatch between alpha and the number of generators for the (error-free) initial set.")
    %elseif size(generators(errorSet),2) ~= length(beta)
    %    error("Dimension mismatch between beta and the number of generators for the measurements error set.")
    end
end


K_mod = T.feedback_speedup.K_mod;

U_alpha = T.termRegCtrl{nTimeStep}.U_alpha;
%U_beta = T.termRegCtrl{nTimeStep}.U_beta;
u_c = T.termRegCtrl{nTimeStep}.u_c;

% We can now compute the control input quite easily
%u = K_mod * (y - T.optsInternal.Param.D*currentInput) + U_alpha * alpha + U_beta * beta + u_c;
u = K_mod * y + U_alpha * alpha + u_c;





end