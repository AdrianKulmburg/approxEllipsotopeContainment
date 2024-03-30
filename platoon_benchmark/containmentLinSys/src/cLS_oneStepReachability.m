function [Rtp, Rti, ctrl_variable, utotal, additional_variables, additional_constraints] = cLS_oneStepReachability(Rtp_previous, utotal_previous, optsInternal)
% CLS_ONESTEPREACHABILITY - Computes one step of linear reachability
%       analysis, as used by the containmentLinSys algorithm
%
% Syntax:
%       [Rtp, Rti, ctrl_variable, utotal, additional_variables, additional_constraints] = cLS_oneStepReachability(Rtp_previous, utotal_previous, optsInternal)
%
% Description:
%       Computes one step of linear reachability analysis, either using the
%       iterative or feedback approach for the controller. This entire
%       algorithm is set up such that the reachability is computed
%       symbolically, which is why we couldn't use CORA directly here.
%
% Input Arguments:
%       -Rtp_previous: Struct containing the symbolic variables describing
%                      the time-point solution from the previous step
%       -utotal_previous: Struct containing the symbolic variables
%                      describing the input from the previous step
%
%       -optsInternal: Struct containing internal parameters
%
% Output Arguments:
%       -Rtp:          Struct containing the symbolic variables describing
%                      the time-point solution
%       -Rti:          Struct containing the symbolic variables describing
%                      the time-interval solution
%
%       - ctrl_variable: Struct containing the symbolic variables necessary
%                      to describe the controller
%
%       - utotal:      Struct containing the symbolic variables describing
%                      the input 
%
%       -additional_variables, additional_constraints:
%                      These are the auxilliary variables and constraints
%                      that are needed to build the optimization problem,
%                      but don't necessarily have an interpretation.
%
% See Also:
%       ---
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
% Copyright (c) 2023 Chair of Robotics, Artificial Intelligence and
%               Embedded Systems, TU Muenchen
%------------------------------------------------------------------

    % More convenient
    reachability = optsInternal.reachability;
    sys = optsInternal.stdSys;
    dims = optsInternal.stdDynDims;


    if strcmp(optsInternal.controlMethod, 'iterative')
        % Initialize the variable for the input
        ctrl_variable.U_alpha = sdpvar(dims.nu, size(Rtp_previous.Y.G,2) - dims.nv, 'full');
        ctrl_variable.U_beta = sdpvar(dims.nu, dims.nv, 'full');
        ctrl_variable.u_c = sdpvar(dims.nu, 1, 'full');

        additional_variables.U_alpha = ctrl_variable.U_alpha;
        additional_variables.U_beta = ctrl_variable.U_beta;
        additional_variables.u_c = ctrl_variable.u_c;

        utotal.G = [ctrl_variable.U_alpha ctrl_variable.U_beta];
        utotal.G = [utotal.G sparse(dims.nu, dims.nv)];
        utotal.c = ctrl_variable.u_c;

        psi = sys.C*reachability.eAminusdt*reachability.Phi*(sys.B*ctrl_variable.u_c+sys.chi);
        Psi_alpha = sys.C*reachability.eAminusdt*reachability.Phi*sys.B*ctrl_variable.U_alpha;
        Psi_beta = sys.C*reachability.eAminusdt*reachability.Phi*sys.B*ctrl_variable.U_beta;


        m = size(Rtp_previous.Y.G,2) + dims.nv;
        Xhom.G = sdpvar(dims.nx, m, 'full');
        Xhom.c = sdpvar(dims.nx, 1, 'full');

        additional_variables.Xhom_G = Xhom.G;
        additional_variables.Xhom_c = Xhom.c;

        m_utotal_previous = size(utotal_previous.G,2);

        
         
        pre_constraints = sys.C * reachability.eAminusdt * Xhom.G == [Rtp_previous.Y.G+[Psi_alpha Psi_beta] sys.F] - [sys.D*utotal_previous.G sparse(dims.ny, m-m_utotal_previous)];
        pre_constraints = [pre_constraints sys.C * reachability.eAminusdt * Xhom.c == Rtp_previous.Y.c + psi - sys.D*utotal_previous.c - sys.gamma];


        Rtp.X.G = [Xhom.G reachability.P_dist];
        Rtp.X.c = Xhom.c;

        Rtp.Yexact.G = sys.C*Rtp.X.G + sys.D*[utotal.G sparse(dims.nu, reachability.m_P_dist)];
        Rtp.Yexact.c = sys.C*Rtp.X.c + sys.D*utotal.c + sys.gamma;

        % Extending
        utotal_hom.G = utotal.G;
        utotal.G = [utotal.G sparse(dims.nu, reachability.m_P_dist)];

        Rtp.YnoInput.G = [sys.C*Rtp.X.G sys.F];
        Rtp.YnoInput.c = sys.C*Rtp.X.c + sys.gamma;

        Rtp.YnoInputNoGamma.G = [sys.C*Rtp.X.G sys.F];
        Rtp.YnoInputNoGamma.c = sys.C*Rtp.X.c;

        Rtp.Y.G = [Rtp.Yexact.G sys.F];
        Rtp.Y.c = Rtp.Yexact.c;

        % Extending again
        utotal.G = [utotal.G sparse(dims.nu, dims.nv)];
    else
        ctrl_variable.U_alpha = sdpvar(dims.nu, size(optsInternal.postH,1), 'full');
        ctrl_variable.U_beta = sdpvar(dims.nu, dims.nv, 'full');
        ctrl_variable.u_c = sdpvar(dims.nu, 1, 'full');

        additional_variables.U_alpha = ctrl_variable.U_alpha;
        additional_variables.U_beta = ctrl_variable.U_beta;
        additional_variables.u_c = ctrl_variable.u_c;

        m_y_previous = size(Rtp_previous.Y.G,2);

        utotal.G = optsInternal.K*Rtp_previous.Y.G + [ctrl_variable.U_alpha*optsInternal.postH ctrl_variable.U_beta sparse(dims.nu, m_y_previous-size(optsInternal.postH,2)-dims.nv)];
        utotal.G = [utotal.G sparse(dims.nu, dims.nv)];
        utotal.c = optsInternal.K*Rtp_previous.Y.c + ctrl_variable.u_c;

        psi = sys.C*reachability.eAminusdt*reachability.Phi*(sys.B*ctrl_variable.u_c+sys.B*optsInternal.K*Rtp_previous.Y.c+sys.chi);
        Psi_alpha = sys.C*reachability.eAminusdt*reachability.Phi*sys.B*ctrl_variable.U_alpha*optsInternal.postH;
        Psi_beta = sys.C*reachability.eAminusdt*reachability.Phi*sys.B*ctrl_variable.U_beta;

        m = size(Rtp_previous.Y.G,2) + dims.nv;
        Xhom.G = sdpvar(dims.nx, m, 'full');
        Xhom.c = sdpvar(dims.nx, 1, 'full');

        additional_variables.Xhom_G = Xhom.G;
        additional_variables.Xhom_c = Xhom.c;

        m_utotal_previous = size(utotal_previous.G,2);

        pre_constraints = sys.C*reachability.eAminusdt*Xhom.G == [(speye(dims.ny)+sys.C*reachability.eAminusdt*reachability.Phi*sys.B*optsInternal.K)*Rtp_previous.Y.G sys.F] + [Psi_alpha Psi_beta sparse(dims.ny,m-size(optsInternal.postH,2)-dims.nv)] - [sys.D*utotal_previous.G sparse(dims.ny, m-m_utotal_previous)];
        pre_constraints = [pre_constraints sys.C * reachability.eAminusdt * Xhom.c == Rtp_previous.Y.c + psi - sys.D*utotal_previous.c - sys.gamma];

        Rtp.X.G = [Xhom.G reachability.P_dist];
        Rtp.X.c = Xhom.c;

        Rtp.YnoInput.G = [sys.C*Rtp.X.G sys.F];
        Rtp.YnoInput.c = sys.C*Rtp.X.c + sys.gamma;

        Rtp.YnoInputNoGamma.G = [sys.C*Rtp.X.G sys.F];
        Rtp.YnoInputNoGamma.c = sys.C*Rtp.X.c;

        Rtp.Yexact.G = sys.C*Rtp.X.G + sys.D*[utotal.G sparse(dims.nu, reachability.m_P_dist)];
        Rtp.Yexact.c = sys.C*Rtp.X.c + sys.D*utotal.c + sys.gamma;

        % Extending
        utotal_hom.G = utotal.G;
        utotal.G = [utotal.G sparse(dims.nu, reachability.m_P_dist)];

        Rtp.Y.G = [Rtp.Yexact.G sys.F];
        Rtp.Y.c = Rtp.Yexact.c;

        % Extending again
        utotal.G = [utotal.G sparse(dims.nu, dims.nv)];

    end
    % Let us extract the homogeneous parts:
    m_X_previous = size(Rtp_previous.X.G, 2);
    H_0 = [Rtp_previous.X.G sparse(dims.nx, m-m_X_previous)]; % Padding with zero values to make sue H_0 and H_1 have the same sizes
    H_0_c = Rtp_previous.X.c;
    H_1 = Xhom.G;
    H_1_c = Xhom.c;

    % Deduce the part that was added by the constant-like parts
    u_const = sys.B * utotal_hom.G;
    u_const_c = sys.B * utotal.c + sys.chi;

    % Now, it's time for the time interval solution. We first compute the
    % curvature error:
    combined_interval = [reachability.F_interval reachability.G_interval];

    combined_matrix = [H_0; u_const];
    combined_center = [H_0_c; u_const_c];
    [C_G, C_c, additional_constraints, additional_variables_multIM] = multIM_zono(combined_interval, combined_matrix, combined_center);

    additional_variables.multIM = additional_variables_multIM;

    % We can now easily construct the time interval solution
    [comb_G, comb_c] = combZonos(H_0, H_1, H_0_c, H_1_c);
    Rti.X.G = [comb_G C_G optsInternal.reachability.P_dist];
    Rti.X.c = comb_c + C_c;

    % We now do a little trick to compute the time interval solution
    % of Yexact
    H_0_Yexact = sys.C*H_0 + sys.D*[utotal_hom.G sparse(dims.nu, size(H_0,2)-size(utotal_hom.G,2))];
    H_0_Yexact_c = sys.C*H_0_c + sys.D*utotal.c;
    H_1_Yexact = sys.C*H_1 + sys.D*utotal_hom.G;
    H_1_Yexact_c = sys.C*H_1_c + sys.D*utotal.c;
    [comb_G_Yexact, comb_c_Yexact] = combZonos(H_0_Yexact, H_1_Yexact, H_0_Yexact_c, H_1_Yexact_c);

    Rti.Yexact.G = [comb_G_Yexact sys.C*C_G sys.C*optsInternal.reachability.P_dist];
    Rti.Yexact.c = comb_c_Yexact + sys.C*C_c + sys.gamma;

    
    Rti.Y.G = [Rti.Yexact.G sys.F];
    Rti.Y.c = Rti.Yexact.c;

    % Merge constraints
    additional_constraints = [pre_constraints additional_constraints];
    
end

function [comb_G, comb_c] = combZonos(zonoGenerators1, zonoGenerators2, center1, center2)
    m1 = size(zonoGenerators1,2);
    m2 = size(zonoGenerators2,2);

    comb_c = 0.5 * (center1 + center2);

    if m1 == m2
        comb_G = [0.5 * (center1 - center2) 0.5 * (zonoGenerators1 + zonoGenerators2) 0.5 * (zonoGenerators1 - zonoGenerators2)];
    else
        error('Mismatch of number of generators. This should not happen, so something deeply wrong happened within the algorithm...')
    end
end

function [G, c, additional_constraints,additional_variables] = multIM_zono(intervalMatrix, zonoGenerators, zonoCenter)
    % center and radius of interval matrix
    M_c = center(intervalMatrix);
    M_r = rad(intervalMatrix);

    n = size(zonoGenerators,1);
    m = size(zonoGenerators,2);

    zonoGenerators_abs = sdpvar(n,m,'full');
    zonoCenter_abs = sdpvar(n,1,'full');

    additional_constraints = [zonoGenerators(:) <= zonoGenerators_abs(:) -zonoGenerators(:) <= zonoGenerators_abs(:)];
    

    additional_constraints = [additional_constraints zonoCenter <= zonoCenter_abs -zonoCenter <= zonoCenter_abs];
    
    % absolute value of zonotope center and generators
    Zabssum = sum([zonoCenter_abs zonoGenerators_abs],2);
    
    G = [M_c*zonoGenerators diag(M_r*Zabssum)];
    additional_variables = {zonoGenerators_abs zonoCenter_abs};

    c = M_c * zonoCenter;
end



