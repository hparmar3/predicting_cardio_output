% Question 6
clc;
close all;

% Run question 5 from which the parlikar algorithm helper function is called
% Change file name in question 5 to produce results for patient s00020,
% s00138 and s00214
Group19_COProject_Q5;

% Convert Compliance to mL/mmHg
Cn = calibrationFactor * 1000 ./ 60; 

% Calculate Rn = Tau_n/Cn
Rn = tau./Cn;

% Convert time axis of 12 h into second
t_s = t * 60;

% Plot Rn vs time
hold off;
figure;
plot(t_s,Rn);
xlabel("time(s)");
ylabel("TPR(mmHg/(mL/s))");
title('Total peripheral resistance of patient 138 over the first 12 hours');
