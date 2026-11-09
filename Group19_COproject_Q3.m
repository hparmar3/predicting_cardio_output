% Task 3
% Extract data for patient 20 by 2analyze
% Could be the same as in Task 1
close all;
clc;

ABP_timeseries = load('s00020-2567-03-30-17-47_ABP.txt');
TCO_data = readmatrix('s00020-2567-03-30-17-47n.txt','NumHeaderLines',2);

fs = 125; % Sampling frequency
ts_12 = 12 * 60 * 60; % 12h in second
t_start_abp = ABP_timeseries(1,1); % Extract initial time value
t_start_TCO = TCO_data(1,1); 
ABP_window_12 = ABP_timeseries(:,1) < (ts_12 + t_start_abp); % 12 hour window in second
TCO_window_12 = TCO_data(:,1) < (ts_12 + t_start_TCO); % 12 hour window in minute

% Extract the data from the first 12h
ABP_timeseries_12h = ABP_timeseries(ABP_window_12,:);
TCO_data_12h = TCO_data(TCO_window_12,:);

% Extract the ABP pressure column
abp_pressure_12 = ABP_timeseries_12h(:,2);

%% Use helper functions from 2analyze to extract feature from the first 12h
ABP_onset = wabp(abp_pressure_12);  % Detect ABP onset
ABP_feature = abpfeature(abp_pressure_12,ABP_onset); % Extract the 12 ABP features
quality = jSQI(ABP_feature, ABP_onset, abp_pressure_12); % calculate signal quality index

%% Use helper functions from 3estimate to get continuous CO
[CO, t_CO, t_old, feat_filtered] = estimateCO_v3(ABP_onset,ABP_feature,quality,5,2);

%% C2 Calibration
TCO = TCO_data_12h(:,16);

% Extract the index and value of the first non-zero TCO
r_1_index = find(TCO > 0,1,'first');  
r_1 = TCO(r_1_index);

% Convert index into minute
t_first_valid = TCO_data_12h(r_1_index,1)/60;

% Find the estimate at the 1 min window of timestamp
window = (t_CO > t_first_valid-1) & (t_CO <= t_first_valid);
x_1 = median(CO(window)); % Take the median of the window

% Calculate the calibrating index 
k = r_1/x_1;
CO_cal = k * CO;

% Extract PP, MAP, HR
PP = feat_filtered(:,5);
MAP = feat_filtered(:,6);
BP = feat_filtered(:,7); % extract beat period first
HR = 60*125./BP; % calculate HR and convert into bpm

%% Extract episodic thermodilution CO measurements
TCO_valid = TCO(TCO > 0);
index_valid = find(TCO > 0);  
time_TCO_valid = TCO_data_12h(index_valid,1)/60; % timestamp of TCO measurement in minute
HR_TCO = TCO_data_12h(index_valid,2);
PP_TCO = TCO_data_12h(index_valid,3)-TCO_data_12h(index_valid,4);
MAP_TCO = TCO_data_12h(index_valid,5);

%% Plot continuous CO, PP, MAP, and HR across time
figure(1);
subplot(4,1,1);
plot(t_CO,CO_cal); % continuous cardiac output (curve line)
hold on;
stem(time_TCO_valid,TCO_valid,"filled",LineWidth = 1.0) % discrete TCO measurement (stem plots)
xlabel('time (min)');
ylabel('CO');

subplot(4,1,2);
plot(t_CO,PP);
hold on;
stem(time_TCO_valid,PP_TCO,"filled", LineWidth = 1.0);
xlabel('time (min)');
ylabel('PP');

subplot(4,1,3);
plot(t_CO,MAP);
hold on;
stem(time_TCO_valid,MAP_TCO,"filled",LineWidth = 1.0);
xlabel('time (min)');
ylabel('MAP');

subplot(4,1,4);
plot(t_CO,HR);
hold on;
stem(time_TCO_valid,HR_TCO,"filled",LineWidth = 1.0);
xlabel('time (min)');
ylabel('HR');

% Resize the image to be taller
fig = gcf;
fig.Position(4) = fig.Position(4) + 120;

% Add overall legend
Lgnd = legend('continuous estimate curve', 'discrete TCO measurement');
Lgnd.Position(1) = 0.7;
Lgnd.Position(2) = 0.94;

hold off;