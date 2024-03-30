function [Rtp, Rti, ctrl_variable, R, RCI, RCI_zono, RCIState, RCIState_zono] = cLS_inverseStdDynSys(Rtp, Rti, ctrl_variable, R, RCI, RCI_zono, RCIState, RCIState_zono, optsInternal)
% CLS_INVERSESTDDYNSYS - reverses certain transformations performed in
%                        cLS_setOptsInternal.m
%
% Syntax:
%       [Rtp, Rti, ctrl_variable, R, RCI, RCI_zono, RCIState, RCIState_zono] = cLS_inverseStdDynSys(Rtp, Rti, ctrl_variable, R, RCI, RCI_zono, RCIState, RCIState_zono, optsInternal)
%
% Description:
%       Reverses certain transformations performed in cLS_setOptsInternal.m
%
% References:
%       * *[1] "From LQR to Static Output Feedback: a New LMI Approach",
%              Luis Rodrigues, CDC 2022
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



    N = optsInternal.Opts.N;

    % We begin by cutting out the part of the state that is not necessary
    % anymore
    n_additional = optsInternal.stdDynDims.nx - optsInternal.dynDims.nx;

    P = [speye(optsInternal.dynDims.nx) sparse(optsInternal.dynDims.nx, n_additional)];

    for i=1:(N+1)
        Rtp{i}.X.G = Rtp{i}.X.G(1:end-n_additional,:);
        Rtp{i}.X.c = Rtp{i}.X.c(1:end-n_additional);
    end
    for i=1:N
        Rti{i}.X.G = Rti{i}.X.G(1:end-n_additional,:);
        Rti{i}.X.c = Rti{i}.X.c(1:end-n_additional);
    end

    RCIState = P * RCIState;
    RCIState_zono = P * RCIState_zono;

    % We continue with inverting the variable transformations we did
    % because of D
    svd_U = optsInternal.stdSys.svd_U;
    svd_V = optsInternal.stdSys.svd_V;
    for i=1:(N+1)
        Rtp{i}.Y.G = svd_U*Rtp{i}.Y.G;
        Rtp{i}.Y.c = svd_U*Rtp{i}.Y.c;
        Rtp{i}.Yexact.G = svd_U*Rtp{i}.Yexact.G;
        Rtp{i}.Yexact.c = svd_U*Rtp{i}.Yexact.c;
        Rtp{i}.YnoInput.G = svd_U*Rtp{i}.YnoInput.G;
        Rtp{i}.YnoInput.c = svd_U*Rtp{i}.YnoInput.c;
        Rtp{i}.YnoInputNoGamma.G = svd_U*Rtp{i}.YnoInputNoGamma.G;
        Rtp{i}.YnoInputNoGamma.c = svd_U*Rtp{i}.YnoInputNoGamma.c;

        Rtp{i}.X.G = svd_V*Rtp{i}.X.G;
        Rtp{i}.X.c = svd_V*Rtp{i}.X.c;
    end
    for i=1:N
        Rti{i}.Y.G = svd_U*Rti{i}.Y.G;
        Rti{i}.Y.c = svd_U*Rti{i}.Y.c;
        Rti{i}.Yexact.G = svd_U*Rti{i}.Yexact.G;
        Rti{i}.Yexact.c = svd_U*Rti{i}.Yexact.c;

        Rti{i}.X.G = svd_V*Rti{i}.X.G;
        Rti{i}.X.c = svd_V*Rti{i}.X.c;
    end
    R = svd_U * R;
    RCI = svd_U * RCI;
    RCI_zono = svd_U * RCI_zono;

    RCIState = svd_V * RCIState;
    RCIState_zono = svd_V * RCIState_zono;
end