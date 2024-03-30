classdef termRegContainmentLinSys < terminalRegion
% TERMREGCONTAINMENTLINSYS - terminal region object computed with the
%                   containment approach proposed in [1]
%
% Syntax:
%       obj = TERMREGCONTAINMENTLINSYS(dynamics, termReg, termReg_zono, termRegState, termRegState_zono, Rtp, R, s, termRegCtrl, debug, runtimes, optsInternal)
%
% Description:
%       This class represents a terminal region computed with the
%       containment approach in [1].
%
% Input Arguments:
%       -dynamics:      function handle to the dynamic function of the open
%                       loop system
%
%       - termReg:      terminal region (class: zonotope or ellipsoid)
%
%       - termReg_zono: zonotopic terminal region
%
%       - termRegState: terminal region of the state (class: zonotope or
%                       ellipsoid)
%
%       - termRegState_zono: zonotopic terminal region of the state
%
%       - Rtp:          over-approximations of the reachable set at time
%                       steps i * Opts.timeStep, for i=0,...,Opts.N
%
%       - R:            over-approximation of the reachable set of termReg,
%                       assuming each point is steered using the control
%                       inputs from termRegCtrl
%
%       - s:            value of the initial stretching factors
%
%       - termRegCtrl:  set-based safety-presevering controller
%
%
%       - debug:        struct containing internal variables, for debug
%                       purposes
%
%       - runtimes:     struct containing information about the runtimes of
%                       the different parts of the algorithm
%
%
%       - optsInternal: a structure containing internal algorithm settings
%                       See ./src/cLS_setOptsInternal
%
%
%
% Output Arguments:
%       -obj:      resulting object of class termRegContainmentLinSys
%
% See Also:
%       terminalRegion, computeTermRegNonlinSysLinApproach,
%
% References:
%       * *[1] "Approximability of the Containment Problem for Zonotopes
%              and Ellipsotopes", A. Kulmburg et al., submited to TAC
%              in 2024
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


properties (SetAccess = private, GetAccess = public)
    set_zono = [];
    setState = [];
    setState_zono = [];
    Rtp = [];
    R = [];
    s = [];
    termRegCtrl = [];
    debug = [];
    runtimes = [];
    Opts = [];
    optsInternal = [];

    % Additional values for Param
    Y = [];
    Yexact = [];
    Ustart = [];

    C = [];
    D = [];
    F = [];
    gamma = [];


    feedback_speedup = []; % In the case of the feedback algorithm,
    % we can speed up calculations for computing the control input,
    % by pre-computing a few things
end

methods
    function obj = termRegContainmentLinSys(dynamics, termReg, termReg_zono, termRegState, termRegState_zono, Rtp, R, s, termRegCtrl, debug, runtimes, optsInternal)

        % call superclass constructor
        obj = obj@terminalRegion(dynamics,termReg,optsInternal.Param);

        % store parameters
        obj.set_zono = termReg_zono;
        obj.setState = termRegState;
        obj.setState_zono = termRegState_zono;
        obj.Rtp = Rtp;
        obj.R = R;
        obj.s = s;
        obj.termRegCtrl = termRegCtrl;
        obj.debug = debug;
        obj.runtimes = runtimes;
        obj.Opts = optsInternal.Opts;
        obj.optsInternal = optsInternal;

        obj.Y = optsInternal.Param.Y;
        obj.Yexact = optsInternal.Param.Yexact;
        obj.Ustart = optsInternal.Param.Ustart;

        obj.C = optsInternal.Param.C;
        obj.D = optsInternal.Param.D;
        obj.F = optsInternal.Param.F;
        obj.gamma = optsInternal.Param.gamma;

        if strcmp(optsInternal.controlMethod, 'feedback')

            if strcmp(optsInternal.terminalRegionType, 'ellipsoid')
                % Create combined system of terminal region and error
                % set
                G_ell = optsInternal.stdSys.svd_U * optsInternal.stdSys.C * optsInternal.preH * diag(s);
                G_error = optsInternal.Param.F*generators(optsInternal.Param.V);
                c_ell = termRegState.q;
                c_error = optsInternal.Param.F*center(optsInternal.Param.V);

                obj.feedback_speedup.comb_G = [G_ell G_error];
                obj.feedback_speedup.comb_c = optsInternal.Param.C*c_ell + c_error + optsInternal.Param.gamma;

                obj.feedback_speedup.is_bijective = (size(obj.feedback_speedup.comb_G,1) == size(obj.feedback_speedup.comb_G,2)) & (rank(obj.feedback_speedup.comb_G) == size(obj.feedback_speedup.comb_G,1));

                if obj.feedback_speedup.is_bijective
                    % So, if the combined matrix is bijective, we can
                    % do some cool stuff very quickly

                    if cond(obj.feedback_speedup.comb_G) > 1e10
                        warning("The projection matrix for the terminal region has a bad condition number. It is likely that when generating points using the simulate function (or, more generally, when computing the control input), the generated trajectories might have numerical errors.")
                    end

                    obj.feedback_speedup.multFactor = inv(obj.feedback_speedup.comb_G);
                    obj.feedback_speedup.constFactor = - obj.feedback_speedup.multFactor * obj.feedback_speedup.comb_c;
                end
            elseif strcmp(optsInternal.terminalRegionType, 'zonotope')
                obj.feedback_speedup.c_comb = optsInternal.Param.C*center(termRegState) + optsInternal.Param.F*center(optsInternal.Param.V) + optsInternal.Param.gamma;
                G_state = generators(termReg);
                G_error = optsInternal.Param.F*generators(optsInternal.Param.V);
                obj.feedback_speedup.G_comb = [G_state G_error];
            end

            % Finally, we can speed up the computation of the gain
            % matrix multiplication

            K_mod = optsInternal.K * optsInternal.stdSys.svd_U';

            % Saving
            obj.feedback_speedup.K_mod = K_mod;

        end
    end
end
end

