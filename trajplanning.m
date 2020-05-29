%Function to determine upper-limb's optimum trajectory as a Bezier curve, 
%added of a 5th order polynomia to represent hand/weight velocity profile 

%Input parameters:
%V = array with the design variables (joint angles)
%l = array with limb lengths: [arm,forearm, hand]
%rcg = array with proximal center of mass distances: [arm, forearm, hand]
%lmass = array with limb masses: [arm, forearm, hand]
%Q = Cutset matrix requied for the static analysis
%k = stiffness for the torsional springs considered for each DoF
%Neutral = joint`s position for the torsional springs torques to be zero
%m = mass of the weight being carried out by the hand. It should be zero 
%if no weight is being hold by the hand
%r = distance of the weight's center of mass to the hand. It should be
%zero if no weight is being hold by the hand
%LRM = joint's lower range of motion
%URM = joint's upper range of motion
%P0 = trajectory's starting point
%P1 = trajectory'd final point
%v1 = desired linear velocity at the end of the trajectory
%w1 = desired angular velocity at the end of the trajectory
%torques =  object handle to export the torques at the joints
%angles = object handle to export the angles at the joints
%X4 = angles describing the posture at the beginning of the motion
%velocities = object handle to export the velocities a the joints
%torques =  object handle to export the torques at the joints
%angles = object handle to export the angles at the joints
%N = number of points in which the trajectory will be discretized

%Output parameters:
%CostF = the sum of the squared sum of the torques through all trajectory

function [CostF] = trajplanning(V, l, rcg, lmass, Q, k, Neutral, m, r, ...
                   LRM, URM, P0, P1, v1, w1, X4, velocities, torques, ...
                   angles, vel, accel, N)

%Creates a matrix with the initial and end trajectory points 
IniEndP = [P0'; P1'];

%The derivative of the trajectory at P0 is given as the first 3 design 
%variables V(1), V(2) and V(3)
deriv0 = [V(1), V(2), V(3)];

%The derivatinve of the trajectory at P1 is given by the normalized v1,
%multiplied of a scalar V(4)
deriv1 = V(4)*(v1/(norm(v1)));

%Calculates the bezier curve for each one the three axis (x, y and z)
[xbez, xlbez, xllbez] = bezier(IniEndP(:,1), deriv0(1), deriv1(1), N);
[ybez, ylbez, yllbez] = bezier(IniEndP(:,2), deriv0(2), deriv1(2), N);
[zbez, zlbez, zllbez] = bezier(IniEndP(:,3), deriv0(3), deriv1(3), N);

%Determines the profile velocity parameters for a 5th order curve
A = [N^5, N^4, N^3; 5*(N^4), 4*(N^3), 3*(N^2); 20*(N^3), 12*(N^2), 6*N];
b = [norm(v1), V(5), V(6)];
param = A\b';

%Determines the penalty function regarding the direction of deriv1
%V(4) shall always be positive
rhog = 10000;
penalty = 0;
if V(4) <= 0
    penalty = penalty + rhog*(V(4)^2);
end

%Starts the variables that will receive the velocities and accelerations
%for each one the N discretized points of the trajectory
CostF = 0;
v = zeros(N,1);
vel.o = zeros(3, N+1);
a = zeros(N,1);
accel.o = zeros(3, N+1);

%Iterates through all N points of the discretized trajectory
%The iterator starts at 2 as the initial posture is already determined
for i=2:(N+1)
    
    %Calculates the normalized derivative for the point i
    T = [xlbez(i); ylbez(i); zlbez(i)]/norm([xlbez(i); ylbez(i); ...
         zlbez(i)]);
     
    %Calculates the normalized derivative of dr for the point i
    dnormdr = (xlbez(i)*xllbez(i) + ylbez(i)*yllbez(i) + ...
         zlbez(i)*zllbez(i))/norm([xlbez(i); ylbez(i); zlbez(i)]);
     
    %Calculates the dertivative of T for the same point i
    dT = (norm([xlbez(i); ylbez(i); zlbez(i)])*[xllbez(i); yllbez(i); ...
         zllbez(i)] - dnormdr*[xlbez(i); ylbez(i); zlbez(i)])/(...
         norm([xlbez(i); ylbez(i); zlbez(i)])^2);
    
    %Obtains the velocity magnitude (v) and the complete vector (vel)
    v(i-1) = param(1)*((i-1)^5) + param(2)*((i-1)^4) + param(3)*((i-1)^3);
    vel.o(:,i) = v(i-1)*T;
    
    %Obtains the acceleration magnitude (a) and the complete vector (accel)
    a(i-1) = 5*param(1)*((i-1)^4) + 4*param(2)*((i-1)^3) + ...
           +3*param(3)*((i-1)^2);
    accel.o(:,i) = a(i-1)*T + abs(v(i-1))*dT;  
    
    %Defines the call for the optimization of the torques for position i
    fun2 = @(A) staticsdavies(A, l, rcg, ...
                     lmass, Q, k, Neutral, m, r, LRM, URM, torques, ...
                     angles, 2, i, [xbez(i); ybez(i); zbez(i)], ...
                     vel.o(:,i), accel.o(:,i)); 
    
    %Configures the number of iterations and calls the optimization
    %The optimization initial position as given as the angles calculated
    %for the previous point i
    options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000);
    [X4,fval4,flag4,out4] = fmincon(fun2,X4,[],[],[],[],[],[],[],options);
    
    %The cost function is given as the sum of the squared sum of torques
    %for each point i
    CostF = CostF + fval4;
    
end

end