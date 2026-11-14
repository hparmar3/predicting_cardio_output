
%% Task 4: Compare three CO estimation algorithms

close all;
clc;

addpath('PhysioToolkitCardiacOutput_MatlabCode/2analyze');
addpath('PhysioToolkitCardiacOutput_MatlabCode/3estimate');

Patients = {'s00020','s00138','s00214'};
StartHours = [0,1,1];
DurationHours = [12,12,12];

EstimatorIDs = [5,2,6];

%% Analyze each patient

for p = 1:3

    % Load ABP and TCO data
    ABP_files = dir(fullfile(Patients{p},'*_ABP.txt'));
    TCO_files = dir(fullfile(Patients{p},'*n.txt'));

    ABP_timeseries = load(fullfile(Patients{p},ABP_files(1).name));
    TCO_data = readmatrix(fullfile(Patients{p},TCO_files(1).name),'NumHeaderLines',2);

    fs = 125;

    % Select the analysis period
    t_start_abp = ABP_timeseries(1,1);
    t_start = t_start_abp + StartHours(p)*3600;
    t_end = t_start + DurationHours(p)*3600;

    ABP_window = ABP_timeseries(:,1) >= t_start & ABP_timeseries(:,1) < t_end;

    ABP_timeseries = ABP_timeseries(ABP_window,:);
    abp_pressure = ABP_timeseries(:,2);

    %% Extract ABP features

    ABP_onset = wabp(abp_pressure);
    ABP_feature = abpfeature(abp_pressure,ABP_onset);
    quality = jSQI(ABP_feature,ABP_onset,abp_pressure);

    %% Extract TCO measurements

    TCO = TCO_data(:,16);
    time_TCO = TCO_data(:,1)/60;

    index_valid = find(TCO > 0 & isfinite(TCO) & time_TCO >= StartHours(p)*60 & time_TCO < (StartHours(p)+DurationHours(p))*60);

    TCO_valid = TCO(index_valid);
    time_TCO_valid = time_TCO(index_valid);

    %% C2 calibration reference

    r_1_index = index_valid(1);
    r_1 = TCO(r_1_index);
    t_first_valid = time_TCO(r_1_index);

    
    

    %% Estimate and calibrate CO

    Colors = [0.0000 0.4470 0.7410;   % Blue: Liljestrand
              0.1800 0.5500 0.3400;   % Green: Windkessel
              0.9290 0.6940 0.1250];  % Yellow: Herd

    EstimatorNames = {'Liljestrand','Windkessel','Herd'};

    figure;
    hold on;

    for a = 1:3

        % Estimate continuous CO
        [CO,t_CO,t_old,feat_filtered] = estimateCO_v3(ABP_onset,ABP_feature,quality,EstimatorIDs(a),2);

        % Time in minutes
        t_CO = t_CO + StartHours(p)*60;

        % C2 calibration
        window = (t_CO > t_first_valid-1) & (t_CO <= t_first_valid);

        x_1 = median(CO(window));
        k = r_1/x_1;
        CO_cal = k*CO;

        % Plot CO estimate
        plot(t_CO,CO_cal,'Color',Colors(a,:), 'LineWidth',0.5);

    end

    % Plot discrete TCO measurements
    stem(time_TCO_valid,TCO_valid,'filled', 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1.0);

    xlabel('time (min)');
    ylabel('CO (L/min)');
    title(['Patient ',Patients{p},': CO Estimator Comparison']);

    xlim([0 700]);
    xticks(0:100:700);

    legend('Liljestrand','Windkessel','Herd', 'discrete TCO measurement','Location','best');

    grid on;

    fig = gcf;
    fig.Position(4) = fig.Position(4) + 120;

    hold off;

    saveas(gcf,['Q4_CO_',Patients{p},'.png']);

end
