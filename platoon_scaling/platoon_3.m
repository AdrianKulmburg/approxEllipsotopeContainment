function f = platoon_3(x,u,w)
% PLATOON_3 - dynamic equations for the car benchmark model, with 3
%    vehicles
%
% Syntax:
%       f = PLATOON_3(x,u,w)
%
% Description:
%
%       The benchmarks represents a platoon with 3 vehicles, where the
%       states are the relative positions and velocities betwenn the
%       vehicles, and the system inputs are the vehicle accelerations (see
%       Sec. IV in [1]).
%
%
% Input Arguments:
%
%       -x:     system states x (dimension: [nx,1])
%       -u:     control inputs u (dimension: [nu,1]) 
%       -w:     disturbances w (dimension: [nw,1])
%
% Output Arguments:
%
%       -f:     value of the dynamic function = time derivative of x
%
% See Also:
%       ---
%
% References:
%       * *[1] Schuermann et al. (2017)*, Optimal Control of Sets of 
%              Solutions to Formally Guarantee Constraints of Disturbed 
%              Linear Systems
%

N_vehicles = 3;

f(1,1) = x(2);
f(2,1) = u(1) + w(1);

for i=2:N_vehicles
    f(2*(i-1)+1,1) = x(2*(i-1)+2);
    f(2*(i-1)+2,1) = u(i-1)-u(i)+w(i-1)-w(i);
end