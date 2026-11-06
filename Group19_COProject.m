% QUESTION 1
clear;
clc;
close all;

% Load the ABP signal
% Add the folder with the MATLAB functions
addpath('PhysioToolkitCardiacOutput_MatlabCode/2analyze');

% load the ABP text file
ABP = load(fullfile('s00020', ...
    's00020-2567-03-30-17-47_ABP.txt'));

% Make ABP is a column vector
ABP = ABP(:);

% Sampling frequency
Fs = 125;

% Display information
fprintf('Number of ABP samples: %d\n', length(ABP)); %[output:94efade5]
fprintf('Recording duration %.2f hours\n', ... %[output:group:42840193] %[output:6e4d1d14]
    length(ABP)/Fs/3600); %[output:group:42840193] %[output:6e4d1d14]

% Detect heartbeats
% Use the WABP algorithm
OnsetTimes = wabp(ABP);
fprintf('Number of detected pulse onsets: %d\n', ... %[output:group:2024614b] %[output:676cb0f8]
    length(OnsetTimes)); %[output:group:2024614b] %[output:676cb0f8]

% Extract ABP waveform features
% Use the ABPFEATURE algorithm




%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:94efade5]
%   data: {"dataType":"text","outputData":{"text":"Number of ABP samples: 14265000\n","truncated":false}}
%---
%[output:6e4d1d14]
%   data: {"dataType":"text","outputData":{"text":"Recording duration 31.70 hours\n","truncated":false}}
%---
%[output:676cb0f8]
%   data: {"dataType":"text","outputData":{"text":"Number of detected pulse onsets: 74977\n","truncated":false}}
%---
