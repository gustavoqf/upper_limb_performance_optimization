function [TL] = torquelimits(T, Alpha)
%TORQUELIMITS Summary of this function goes here
%   Detailed explanation goes here

%Coefficients for positive torque equations
%Flexion, adduction and pronation
P = [-16.891,   26.353,     -25.035,    82.591;
     -16.543,   -47.486,    -54.837,    63.89;
     24.294,    -56.006,    26.478,     63.387;
     -32.867,   71.19,      10.304,     32.146;
     1.6004,    -9.2107,    10.504,     9.186;
     0,         -0.6262,    2.4011,     10.116;
     10.639,    0,          -6.3734,    5.46];

%Coefficients for negative torque equations
%Extension, abduction and supination
N = [19.907,    -37.001,    -31.03,     -57.454;
     6.7209,    5.8784,     23.225,     -71.628;
     -11.999,   54.763,     -71.894,    -21.532;
     -9.1633,   37.357,     -41.288,    -26.13;
     3.5805,    -16.672,    15.116,     -8.1468;
     0,         0.3036,     1.4903,     -6.5318;
     17.736,    2.4744,     4.4288,     -10.943];
 
for i=1:length(T)
    if T(i) > 0
        TL(i) = P(i,1)*(Alpha(i)^3) + P(i,2)*(Alpha(i)^2) + ...
                P(i,3)*Alpha(i) + P(i,4);
    else
        TL(i) = N(i,1)*(Alpha(i)^3) + N(i,2)*(Alpha(i)^2) + ...
                N(i,3)*Alpha(i) + N(i,4);
    end
end

TL = TL';        
        
end

