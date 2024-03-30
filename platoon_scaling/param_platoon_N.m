function Param = param_platoon_N(N)
% PARAM_PLATOON - parameters for the platoon benchmark
%
% Syntax:
%       Param = PARAM_PLATOON()
%
% Description:
%       Parameters for the platton benchmark. The parameters include input
%       constraints, disturbances as well as the parameters of the motion
%       primitive like initial set, goal state, etc..
%
% Output Arguments:
%
%       -Param:             a structure containing following options
%
%           -.R0:           initial set of states (class: interval)
%           -.xf:           goal state
%           -.tFinal:       final time after which the goal state should be
%                           reached
%           -.U:            set of admissible control inputs (class:
%                           interval)
%           -.W:            set of uncertain disturbances (class: interval 
%                           or zonotope)
%           -.V:            set of measurement errors (class: interval or
%                           zonotope)
%           -.X:            set of state constraints (class: mptPolytope)
%
% See Also:
%       car, param_car_driveStraight
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
% Authors:      Niklas Kochdumper
% Website:      <a href="http://aroc.in.tum.de">aroc.in.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Copyright (c) 2019 Chair of Robotics, Artificial Intelligence and
%               Embedded Systems, TU Muenchen
%------------------------------------------------------------------
    
    
    % set of admissible control inputs
    width = 10 .* ones([N 1]);
    Param.U = interval(-width,width);
    
    % set of uncertain disturbances
    width = ones([N 1]);
    Param.W = interval(-width,width);

    % set of state constraints
    A = [zeros([N-1 2]) kron(eye(N-1),[-1 0])];
    b = zeros([N-1 1]);

    Param.X = polytope(A,b);

end