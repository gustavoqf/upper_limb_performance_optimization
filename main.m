%Development of a 7 DoF upper-limb model for performance optmization
%Author: Gustavo Queiroz Fernandes
%ID number: 201705682
%Advisor: Professor Daniel Martins
%Coadvisor: Professor Leonardo Mejia
%Federal University of Santa Catarina

%cleans Matlab's memory and workspace
clc
clear all
close all

%Defines trajectory's final coordinates [m]
%Inform an empty array for optimization of final position and the region
%of the space in which the the final of the trajectory shall be located
finalP = [0.30; 0.7; 0.0];
%finalP = [];
minfinalP = [0.2; 0.6; -0.1];
maxfinalP = [0.4; 0.8; 0.1];

%Defines trajectoryes initial posture [rad]
%Inform an empty arra for optimization of the initial posture: X1 = []
X1 = [];

%Input parameters regarding athletes physiological properties 
bmass = 70; %Total body mass [kg]
l = [0.32; 0.25; 0.19]; %limb lengths [m]: [arm, forearm, hand]
rcg = [0.436; 0.430; 0.506]; %proximal center of mass distances [m]
lmass = bmass*9.81.*[0.028; 0.016; 0.006]; %limb masses [N]
LRM = [25; -130; -15; 20; 70; -70; -20]; %lower range of motion (degrees)
URM = [95; 0; 45; 145; 170; 0; 30]; %upper range of motion (degrees)
k = (10E-5).*[0.7; 0.7;0.7;0.7;0.1;0.1;0.1]; %torsional spring stiffiness
Neutral = [0.890; -0.524; 0.943; 0.7856; 1.606; 0.855;0]; %neutral posture 

%Input parameters regarding information about the sport
m = 0.625; %mass of the weight begin hold by the hand [kg]
r = 22.86/200; %distance from weight's center of mass to the hand [m]
g= 9.81; %gravity [m/s^2]

%Defines information regarding the basketball specifically
Theta = 44.48*pi/180; %angle the ball arrives into the basket [rad]
P = [5,(3.048-1.5),0]; %coordinates of the basket in the xyz space [m,m,m]

%Converts the lower and upper ranges of motion to radian
LRM = LRM.*pi/180; 
URM = URM.*pi/180;
 
%Defines the number of segments in which the trajectory is discretized
%and creates the object handles to store the results
N = 40;
velocities = hObj([]);
torques = hObj([]);
angles = hObj([]);
trajectory = hObj([]);
vel = hObj([]);
accel = hObj([]);
velocities.o = zeros(7,N+1);
torques.o = zeros(7,N+1);
angles.o = zeros(7,N+1);
trajectory.o = zeros(1, 6);
vel.o = zeros(3, N+1);
accel.o = zeros(3, N+1);

%Defines the fundamental Circuit matrix
B0 = [1 1 1 1 1 1 1];

%Defines the degrees of freedom for the joints: [g, a, b, c, d, e, f]
DF = [1 3 2 2 3 1 1];

%Expands each column of B0 into the Circuit Matrix, B, according to the 
%DoFs of each one of the joints 
[nrows,ncolumns]=size(B0);
column = 1;
for i=1:ncolumns
    for n=1:DF(i)
        B(:,column)=B0(:,i);
        column = column + 1;
    end    
end

%Defines the fundamental Cutset matrix
Q0 = [-1  -1  -1  -1  1  0  0;...
       0  -1  -1  -1  0  1  0;...
       0   0  -1  -1  0  0  1];

%Defines the degrees of constraints: [m1, m2, m3, d, a, b, c]  
DC = [1 1 1 6 9 8 8];

%Expands each column of the cutset matrix according to the degree of 
%constraint of each joint
[nrows,ncolumns]=size(Q0);
column = 1;
for i=1:ncolumns
    for n=1:DC(i)
        Q(:,column)=Q0(:,i);
        column = column + 1;
    end    
end

%If trajectory's initial posture was not given
if (isempty(X1) == 1)

    %Defines the call for the optimization of the initial position
    %It assumes the upper-limp starts at rest: type = 1
    fun = @(A) staticsdavies(A, l, rcg, lmass, Q, k, Neutral, ...
          m, r, LRM, URM, torques, angles, 1, 1);
    
    %Defines optimization's starting point
    I0 = [55, -45, -5, 140, 110, -45, 10]*pi/180;
    
    %Configures the number of iterations and calls the optimization
    options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000);
    [X1,fval1,flag1,out1] = fmincon(fun,I0,[],[],[],[],[],[],[],options);
   
end

%Stores the initial posture in the angles object handle
angles.o(:,1) = X1;

%Calculates screws location and directions based in the angles X1
[s, s0] = screw(l,rcg,X1,r);

%Defines the location of the hand/weight for the initial position
s0d0 = s0(:,4);

%If trajectory's final position was not given
if (isempty(finalP) == 1)
    
    %Defines the call for the optimization of the final position
    fun = @(A) optimumfinalP(A, l, rcg, lmass, Q, k, Neutral, P, ...
          Theta, m, r, LRM, URM, s0d0, X1, minfinalP, maxfinalP, ...
          velocities, torques, angles, trajectory, vel, accel, N);
    
    %Defines optimization starting point
    F0 = (minfinalP + maxfinalP)'/2;
    %F0 = [0.30, 0.7, 0.0];
    
    %Configures the number of iterations and calls the optimization
    options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000);
    [X2,fval2,flag2,out2] = fmincon(fun,F0,[],[],[],[],[],[],[],options);

%If trajectory's final position was given
else
    
    %Checks the required torques and velocities for the weight
    [v, w] = parabola(finalP, P, Theta, m);
    
    %Defines the call for the optimization of the trajectory
    fun = @(T) trajplanning(T, l, rcg, lmass, Q, k, Neutral, m, ...
          r, LRM, URM, s0d0, finalP, v, w, X1, velocities, torques, ...
          angles, vel, accel, N);

    %Defines optimization starting point
    T0 = [0, 0, 0, 1, 0.05, -0.75];  
    T0 = [0, 0, 0.05, 1, 0.05, -0.75]; 
    I0 = [0.05, -0.25, -0.4, 1, 0.05, -0.75];
    %T0 = [-0.2, 0.05, 0.1, 1, 0.05, -0.75];   
    
    %Configures the number of iterations and calls the optimization
    options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000);
    [X2,fval2,flag2,out2] = fmincon(fun,T0,[],[],[],[],[],[],[],options);

    %Stores the trajectory results into its respective object handle
    trajectory.o = X2;

end

%Iterates through all columns of the angles object handle
for i=2:N+1

    %Calculates screw directions and positions for the angles of column i
    [s, s0] = screw(l, rcg, angles.o(:,i), r);
    
    %Assembles the twists into the matrix sm
    sm = zeros(6,13);
    sm(:,1) = [s(:,1);cross(s0(:,1),s(:,1))];       %sm1
    sm(:,2) = [s(:,2);cross(s0(:,1),s(:,2))];       %sm2
    sm(:,3) = [s(:,3);cross(s0(:,1),s(:,3))];       %sm3
    sm(:,4) = [s(:,4);cross(s0(:,2),s(:,4))];       %sm4
    sm(:,5) = [s(:,5);cross(s0(:,2),s(:,5))];       %sm5
    sm(:,6) = [s(:,7);cross(s0(:,3),s(:,7))];       %sm6
    sm(:,7) = [s(:,8);cross(s0(:,3),s(:,8))];       %sm7
    sm(:,8) = [[1;0;0];cross(s0(:,4),[1;0;0])];     %sdx
    sm(:,9) = [[0;1;0];cross(s0(:,4),[0;1;0])];     %sdy
    sm(:,10) = [[0;0;1];cross(s0(:,4),[0;0;1])];    %sdz
    sm(:,11) = [0;0;0;[0;1;0]];                     %se
    sm(:,12) = [0;0;0;[0;0;1]];                     %sf
    sm(:,13) = [0;0;0;[1;0;0]];                     %sg

    %Rearranges the columns of sm to create the Motion matrix Md 
    Md = [sm(:,13),sm(:,[1:12])];

    %Determines the Network Unit Motion matrix (Mn)
    [nrows,ncolumns] = size(B);
    Mn = zeros(nrows*6, ncolumns);
    Mn = Md*diag(B(1,:));
    for j=2:nrows
        Mn = [Mn;Md*diag(B(j,:))];
    end

    %Defines the optimization function for the kinematics problem
    fun2 = @(V) kineticsdavies(V, Mn, vel.o(:,i), [0;0;0], i, velocities);
    
    %Configures the number of iterations and calls the optimization
    options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000);
    [X5,fval5,flag5,out5] = fmincon(fun2,0,[],[],[],[],[],[],[],options);

end

    

