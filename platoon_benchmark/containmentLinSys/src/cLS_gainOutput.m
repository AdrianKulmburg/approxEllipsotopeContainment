function K = cLS_gainOutput(A, B, C, D, Q, R)
% CLS_GAINOUTPUT - computes the gain matrix using the method from [1]
%
% Syntax:
%       K = cLS_gainOutput(A, B, C, D, Q, R)
%
% Description:
%       For a linear system
%       \dot{x} = Ax + Bu
%       y       = Cx + Du
%       this function computes a suitable gain matrix K to control the
%       system by means of u = -Ky.
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

    nu = size(B,2);
    ny = size(C,1);

    % We first compute P via LQR
    [~,P,~] = lqr(full(A), full(B), full(Q), full(R));

    P = sparse(P);

    Rinv = inv(R);

    K = sdpvar(nu, ny);

    N = P*B*(K*D*Rinv+Rinv*D'*K')*B'*P;
    M = Q - P*B*K*C - C'*K'*B'*P+N;

    B_bar = B*(speye(nu) + K*D);
    %R_bar = (speye(nu) - K*D)*Rinv*(eye(nu) - D'*K');
    R_bar = Rinv - K*D*Rinv  - Rinv*D'*K';

    cost = [];
    constraint = [M P*B_bar C'*K'; B_bar'*P R sparse(nu, nu); K*C sparse(nu, nu) R_bar] >= 0;

    options = sdpsettings('allownonconvex',0,'solver','mosek','verbose',0);
    diagnostics = optimize(constraint,cost,options);

    if diagnostics.problem ~= 0
        warning("The solver had problems computing the gain matrix. It is thus very likely that the result will be inadequate. Here is some more detailed information:"+sprintf('\n')+ diagnostics.info)
    end

    K = value(K);
end