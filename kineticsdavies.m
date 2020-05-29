%Function to carry out the kinematic analysis of the 7 DoF upper-limb. ]
%This function calculates V(1) that results in the minimum squared sum of 
%the 7 velocities V

%Input parameters:
%V = array with the design variables (joints angular velocities)
%Mn = Network Unit Motion matrix 
%v = linear velocities desired at the hand: [vx, vy, vz]
%w = angular velocities desired at the hand: [wx, wy, wz]
%N = number of the discretized point of the trajectory considered
%velocities = object handle to export the velocities at the joints

%Output parameters:
%CostF = the squared sum of the angular velocities

function [CostF] = kineticsdavies(V, Mn, v, w, N, velocities)

%Obtains the equations for the optimization problem by moving the primary 
%variables to the right side of the system of equations
Meq = Mn(:,[3:8]);
beq = -w(1).*Mn(:,9)-w(2).*Mn(:,10)-w(3).*Mn(:,11)-v(1).*Mn(:,1)+...
      -v(2).*Mn(:,12)-v(3).*Mn(:,13)-V(1)*Mn(:,2);
vel = [V(1); inv(Meq)*beq];
      
%Determines the cost function to be minimized
CostF = sumsqr(vel);
velocities.o(:,N) = vel;

end

