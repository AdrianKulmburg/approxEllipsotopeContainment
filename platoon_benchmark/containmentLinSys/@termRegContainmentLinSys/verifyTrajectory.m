function verifyTrajectory(T, trajectories, varargin)
% VERIFYTRAJECTORY - given a set of trajectories, checks whether
%                      they all fulfill the state constraints, and
%                      start and end inside the terminal region
%
% Syntax:
%       verifyTrajectory(T, trajectories)
%       verifyTrajectory(T, trajectories, tol)
%
% Description:
%       Given a set of trajectories (for example given as the output of
%       simulateRandom), checks whether they all fulfill the
%       state/input/output constraints, and start and end inside the
%       terminal region.
%
% Input Arguments:
%
%       -trajectories: cell array of structs as generated by simulate.m
%
%       -tol:          tolerance for containment
%                      (class: non-negative double)
%
% Output Arguments:
%       ---
%
% See Also:
%       simulate, simulateRandom
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

if isempty(varargin)
    tol = 1e-5;
elseif length(varargin) == 1
    tol = varargin{1};
elseif length(varargin) > 1
    error("Too many input arguments.")
end
inputArgsCheck({{tol,'att','numeric',{'scalar','nonnegative','nonnan','finite'}}});


inputArgsCheck({{trajectories,'att',{'cell','struct'}}});
if isa(trajectories, 'struct')
    trajectories = {trajectories};
end

N = length(trajectories);
fields = {'x','y','yexact','u'};

for i=1:N
    trajectory = trajectories{i};
    for i_f = 1:length(fields)
        if ~isfield(trajectory, fields{i_f})
            warning("For trajectory number " + num2str(i) + ", the field " + fields{i_f} + " is missing.")
        end
    end

    L = size(trajectory.x, 2);
    for i_f = 1:length(fields)
        str = ['currentField = trajectory.', fields{i_f}, ';'];
        eval(str);
        if size(currentField,2) ~= L
            warning("For trajectory number " + num2str(i) + ", the field x does not have the same number of time steps as the field " + fields{i_f} + ".")
        end
    end
end


termReg = T.set;
termRegState = T.setState;

measurement_errors = T.optsInternal.Param.F * T.optsInternal.Param.V;

X = T.X;
Y = T.Y;
Yexact = T.Yexact;

U = T.U;
Ustart = T.Ustart;

for i = 1:N
    trajectory = trajectories{i};
    x = trajectory.x;
    y = trajectory.y;
    yexact = trajectory.yexact;
    u = trajectory.u;

    % Check that each trajectory starts and ends in the terminal region
    if ~termRegState.contains(x(:,1),'exact',tol)
        warning("For trajectory number " + num2str(i) + ", the initial value of x is not contained in the (state) terminal region.")
    elseif ~termRegState.contains(x(:,end),'exact',tol)
        warning("For trajectory number " + num2str(i) + ", the final value of x is not contained in the (state) terminal region.")
    end
    if ~aux_pointContainment_sumEllipsotopeZonotope(y(:,1), termReg + T.optsInternal.Param.D*u(:,1), measurement_errors, tol)
        warning("For trajectory number " + num2str(i) + ", the initial value of y is not contained in the terminal region.")
    elseif ~contains(termReg+T.optsInternal.Param.D*u(:,end), y(:,end),'exact',tol)
        warning("For trajectory number " + num2str(i) + ", the final value of y is not contained in the terminal region.")
    end

    % We check the initial value ustart
    if ~Ustart.contains(u(:,1),'exact',tol)
        warning("For trajectory number " + num2str(i) + ", the initial value of u (i.e., ustart) is not contained in the set Param.Ustart.")
    elseif ~U.contains(u(:,1),'exact',tol)
        warning("For trajectory number " + num2str(i) + ", the initial value of u (i.e., ustart) does not fulfill the input constraints in Param.U. This is tolerated, but please make sure that this is intended.")
    end

    % We check the final value of the input
    if ~Ustart.contains(u(:,end),'exact',tol)
        warning("For trajectory number " + num2str(i) + ", the final value of u is not contained in the set Param.Ustart.")
    end

    timeSteps = size(x, 2);
    for t = 1:timeSteps
        % Verify state constraints
        if ~X.contains(x(:,t),'exact',tol)
            warning("For trajectory number " + num2str(i) + ", and time step " + num2str(t) + ", the state x does not satisfy the constraints Param.X.")
        end
        if ~Y.contains(y(:,t),'exact',tol)
            warning("For trajectory number " + num2str(obj.i) + ", and time step " + num2str(t) + ", the output y does not satisfy the constraints Param.Y.")
        end
        if ~Yexact.contains(yexact(:,t),'exact',tol)
            warning("For trajectory number " + num2str(i) + ", and time step " + num2str(t) + ", the (exact) output yexact does not satisfy the constraints Param.Yexact.")
        end
        if ~U.contains(u(:,t),'exact',tol)
            warning("For trajectory number " + num2str(i) + ", and time step " + num2str(t) + ", the output u does not satisfy the constraints Param.U.")
        end
    end
end
end

function res = aux_pointContainment_sumEllipsotopeZonotope(p, E, Z, tol)
    if isa(E, 'zonotope')
        res = contains(E+Z, p, 'exact', tol);
    else
        c_comb = center(E) + center(Z);
        G_comb = [generators(E) generators(Z)];

        m = size(G_comb,2);
        m_Z = size(generators(Z),2);
        m_E = m-m_Z;

        H = 2.*blkdiag(speye(m_E), sparse(m_Z,m_Z));
        f = sparse(m,1);

        A = [sparse(m_Z,m_E) speye(m_Z); sparse(m_Z, m_E) -speye(m_Z)];
        b = ones([2*m_Z 1]);

        Aeq = G_comb;
        beq = sparse(p - c_comb);


        persistent options_quadprog
        if isempty(options_quadprog)
            options_quadprog = optimoptions('quadprog', 'Display', 'none');
        end
        [preimage, fval, exitflag] = quadprog(H,f,A,b,Aeq,beq,[],[],options_quadprog);

        tol = 1e-5;

        % If the problem is not feasible, y has been chosen incorrectly
        if exitflag ~= 1 || fval > 1 + tol
            res = false;
        else
            res = true;
        end

    end
end