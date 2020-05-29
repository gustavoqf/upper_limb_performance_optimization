%Function to vary the final Bezier point and determine the optimal 
%tracjetory (minimum squared sum of torques). 

%Input parameters:
%V = array with the design variables (joint angles)
%l = array with limb lengths: [arm,forearm, hand]
%rcg = array with proximal center of mass distances: [arm, forearm, hand]
%lmass = array with limb masses: [arm, forearm, hand]
%Q = Cutset matrix requied for the static analysis
%k = stiffness for the torsional springs considered for each DoF
%Neutral = joint`s position for the torsional springs torques to be zero
%P = coordinates describing the position of the basket in the space
%Theta = angle of arrival of the ball into the basket
%m = mass of the weight being carried out by the hand. It should be zero 
%if no weight is being hold by the hand
%r = distance of the weight's center of mass to the hand. It should be
%zero if no weight is being hold by the hand
%LRM = joint's lower range of motion
%URM = joint's upper range of motion
%s0d0 = desired position for hand/weight at the beginning of the motion
%Ini = angles describing the posture at the beginning of the motion
%velocities = object handle to export the velocities a the joints
%torques =  object handle to export the torques at the joints
%angles = object handle to export the angles at the joints
%N = number of points in which the trajectory will be discretized

%Output parameters:
%CostF = the sum of the squared sum of the torques through all trajectory

function [CostF] = optimumfinalP(V, l, rcg, lmass, Q, k, Neutral, ...
                   P, Theta, m, r, LRM, URM, s0d0, Ini, minfinalP, ...
                   maxfinalP, velocities, torques, angles, ...
                   trajectory, vel, accel, N)

%Considers the actual point for shooting
s0d1 = V';

%Defines the penalty varibables for the cost function
penalty = 0;
rhog = 100000;

%Updates the required velocities for the ball to reach the target
%For any other task being optimized, the values of v and w shall be updated
%accordingly and the input P and Theta changed as necessary
[v,w] = parabola(s0d1,P,Theta,m);

%Defines the call for the optimization of the trajectory
fun = @(T) trajplanning(T, l, rcg, lmass, Q, k, Neutral, m, r, LRM, ...
      URM, s0d0, s0d1, v, w, Ini, velocities, torques, angles, vel, ...
      accel, N);

%Defines optimization starting point
I0 = [0.05, -0.25, -0.4, 1, 0.05, -0.75];
T0 = [0, -0.2, -0.45, 1, 0.05, -0.75]; 

%Configures the numer of iterations and calls the optimization
options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000, ...
    'DiffMinChange', 0.1);
[X3, fval3, flag3, out3] = fmincon(fun,I0,[],[],[],[],[],[],[],options);

%Computes the penalization
g = [minfinalP - s0d1; s0d1 - maxfinalP];
for i=1:length(g);
    if g(i) > 0
        penalty = penalty + rhog*(g(i)^2);
    end
end

%Makes CostF equal to the value returned the trajectory planning algorithm
CostF = fval3 + penalty;

%Updates the trajectory parameters returned by this function
trajectory.o = X3;

end

