%Function to carry out the static analysis of the 7DoF upper-limb. It
%considers either situations in which the arm is at rest or cases when the
%arm is moving, holding or not some weight at the hand

%Input parameters:
%V = array with the design variables (joint angles)
%l = array with limb lengths: [arm,forearm, hand]
%rcg = array with proximal center of mass distances: [arm, forearm, hand]
%lmass = array with limb masses: [arm, forearm, hand]
%Q = Cutset matrix required for the static analysis
%k = stiffness for the torsional springs considered for each DoF
%Neutral = joint`s position for the torsional springs torques to be zero
%m = mass of the weight being carried out by the hand. It should be zero 
%if no weight is being hold by the hand
%r = distance of the weight's center of mass to the hand. It should be
%zero if no weight is being hold by the hand
%LRM = joint's lower range of motion
%URM = joint's upper range of motion
%torques =  object handle to export the torques at the joints
%angles = object handle to export the angles at the joints
%type = type of the analysis. It is 1 for the arm at rest and 2 otherwise

%Optional input parameters, when the arm is moving(type equal 2):
%N = number of the discretized point of the trajectory considered
%s0d = desired position for hand or object being moved if this is the case
%v = velocity at hand or object being moved if this is the case
%acel = acceleration at hand or object being moved if this is the case

%Output parameters:
%CostF = the sum of the optimum torques 

function [CostF] = staticsdavies(V, l, rcg, lmass, Q, k, Neutral, ...
                   m, r, LRM, URM, torques, angles, type, N, s0d, v, acel)

%Transforms the design variables into a column array
Alpha = V';

%Initiates the penalty variables for the cost function
rhogangle = 1000;
rhohtraj = 1000000;
rhogvar = 100;
penalty = 0;

%Defines the inequality constraints for the range of motion
gstat=[LRM(1)-Alpha(1); LRM(2)-Alpha(2); LRM(3)-Alpha(3); ...
       LRM(4)-Alpha(4); LRM(5)-Alpha(5); LRM(6)-Alpha(6); ...
       LRM(7)-Alpha(7); Alpha(1)-URM(1); Alpha(2)-URM(2); ...
       Alpha(3)-URM(3); Alpha(4)-URM(4); Alpha(5)-URM(5); ...
       Alpha(6)-URM(6); Alpha(7)-URM(7)];

%Calculates the penalty the inequality constraints above
for i=1:length(gstat);
    if gstat(i) > 0
        penalty = penalty + rhogangle*(gstat(i)^2);
    end
end

%Computes the direction and position vectors
[s, s0] = screw(l,rcg, Alpha,r);

%If the arm is at rest
if(type==1)

    %Defines the forces at the hand to hold any weight given
    F = [0,-m*9.81,0];
    
    %Defines the moments at the hand/weight as zero
    %If any other moment is desired, this should be changed
    M = [0,0,0];

%If the arm is moving
elseif(type==2)

    %Calculates the forces required at the hand/weight as the sum of the
    %force to hold any weight and the force due its acceleration
    %If any other force is desired, this should be changed
    F = -m*[acel(1); 9.81+acel(2); acel(3)];
    
    %Defines the moments at the hand/weight as zero
    %If any other moment is desired, this should be changed
    M = [0,0,0];
    
    %Obtains the penalty due the equality constraint of the trajectory
    %Makes sure the hand/weight obeys the desired position s0d
    penalty = penalty + rhohtraj*(sumsqr(s0d - s0(:,4)));
    
    %Obtains the penalty due the equality constraint of the plane
    %Makes sure the hand is always perpendicular to v
    %penalty = penalty + 100*((dot(s(:,9),v))^2);
    
    %Obtains the penalty due the inequality constraints of joints
    %displacements between two consecutive iterations
    %Makes sure the joints does not exceed a maximum acceleration
    gvar = abs(Alpha - angles.o(:,N-1));
    for i=1:length(gvar)
        if gvar(i) >= 0.125
            penalty = penalty + rhogvar*(gvar(i)^2);
        end
    end

end

%Determines the wrenches due limb masses
Ad = zeros(6,34);
sw = [0;1;0];
Ad(:,1) = [cross(s0(:,8),sw);sw];  %Mass of the arm: m1
Ad(:,2) = [cross(s0(:,9),sw);sw];  %Mass of the forearm: m2
Ad(:,3) = [cross(s0(:,10),sw);sw]; %Mass of the hand: m3

%Determines the wrenches at joint d
Ad(:,4) = [cross(s0(:,4),[1; 0; 0]);[1; 0; 0]];  %Produced force: Fx
Ad(:,5) = [cross(s0(:,4),[0; 1; 0]);[0; 1; 0]];  %Produced force: Fy
Ad(:,6) = [cross(s0(:,4),[0; 0; 1]);[0; 0; 1]];  %Produced force: Fz
Ad(:,7) = [[1; 0; 0];0;0;0];  %Produced moment: Mx
Ad(:,8) = [[0; 1; 0];0;0;0];  %Produced moment: My
Ad(:,9) = [[0; 0; 1];0;0;0];  %Produced moment: Mz

%Determines the wrenches at joint a
Ad(:,10) = [cross(s0(:,1),s(:,1));s(:,1)];  %Passive constraint: F1
Ad(:,11) = [cross(s0(:,1),s(:,2));s(:,2)];  %Passive constraint: F2
Ad(:,12) = [cross(s0(:,1),s(:,3));s(:,3)];  %Passive constraint: F3
Ad(:,13) = [s(:,1);0;0;0];  %Active constraint: T1
Ad(:,14) = [s(:,2);0;0;0];  %Active constraint: T2
Ad(:,15) = [s(:,3);0;0;0];  %Active constraint: T3
Ad(:,16) = [s(:,1);0;0;0];  %Active constraint due stiffness: PT1
Ad(:,17) = [s(:,2);0;0;0];  %Active constraint due stiffness: PT2
Ad(:,18) = [s(:,3);0;0;0];  %Active constraint due stiffness: PT3

%Determines the wrenches at joint b
Ad(:,19) = [cross(s0(:,2),s(:,4));s(:,4)];  %Passive constraint: F4
Ad(:,20) = [cross(s0(:,2),s(:,5));s(:,5)];  %Passive constraint: F5
Ad(:,21) = [cross(s0(:,2),s(:,6));s(:,6)];  %Passive constraint: F5x4
Ad(:,22) = [s(:,6);0;0;0];  %Passive constraint: T5x4
Ad(:,23) = [s(:,4);0;0;0];  %Active constraint: T4
Ad(:,24) = [s(:,5);0;0;0];  %Active constraint: T5
Ad(:,25) = [s(:,4);0;0;0];  %Active constraint due stiffness: PT4
Ad(:,26) = [s(:,5);0;0;0];  %Active constraint due stiffness: PT5

%Determines the wrenches at joint c
Ad(:,27) = [cross(s0(:,3),s(:,7));s(:,7)];  %Passive constraint: F6
Ad(:,28) = [cross(s0(:,3),s(:,8));s(:,8)];  %Passive constraint: F7
Ad(:,29) = [cross(s0(:,3),s(:,9));s(:,9)];  %Passive constraint: F6x7
Ad(:,30) = [s(:,9);0;0;0];  %Passive constraint: T6x7
Ad(:,31) = [s(:,7);0;0;0];  %Active constraint: T6
Ad(:,32) = [s(:,8);0;0;0];  %Active constraint: T7
Ad(:,33) = [s(:,7);0;0;0];  %Active constraint due stiffness: PT6
Ad(:,34) = [s(:,8);0;0;0];  %Active constraint due stiffness: PT7

%Determines the Network Unit Action matrix (An)
[nrows,ncolumns] = size(Q);
An = zeros(nrows*6, ncolumns);
An = Ad*diag(Q(1,:));
for j=2:nrows
    An = [An;Ad*diag(Q(j,:))];
end

%Rearragnes An so that it has the following order: Passive constraints,
%active constraints, stiffness active constraints, external forces/moments
%and limb masses
An = An(:,[10:12,19:22,27:30,13:15,23,24,31,32,16:18,25,26,33,34,4:9,1:3]);

%Determines the magnitufes of the passive torques due muscles elasticity
PT = -k.*(Alpha-Neutral);

%Moves the known variables (primary) to the right side of the equation
b = F(1).*An(:,26) + F(2).*An(:,27) + F(3).*An(:,28)...
    + M(1).*An(:,29) + M(2).*An(:,30) + M(3).*An(:,31)...
    - PT(1).*An(:,19) - PT(2).*An(:,20) - PT(3).*An(:,21)...
    - PT(4).*An(:,22) - PT(5).*An(:,23) - PT(6).*An(:,24)...
    - PT(7).*An(:,25) - lmass(1).*An(:,32) - lmass(2).*An(:,33)...
    - lmass(3).*An(:,34);

%Solves the static problem
T = (An(:,[1:18]))\b;

%Isolates the 7 active torques of the shoulder, elbow and wrist
torques.o(:,N) = T(12:end);

%Exports the angles in its object handle
angles.o(:,N) = Alpha;

%Checks the lower and upper torque limits for each one of the joints
TL = torquelimits(torques.o(:,N), Alpha);

%If the arm is at rest
if (type==1) 

    %Computes the sum of torques penalized by max. isometric limits
    CostF = sumsqr(torques.o(:,N)./TL) + penalty; 
    
%If the arm is moving
elseif (type==2)
       
    %Computes the sum of torques
    CostF = (sumsqr(torques.o(:,N))) + penalty; 
end

end

