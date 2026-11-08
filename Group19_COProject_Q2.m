
% QUESTION 2
clear;
clc;
close all;

%% Load the ABP signals
addpath('PhysioToolkitCardiacOutput_MatlabCode/2analyze');

Patients = {'s00138', 's00214'};

Files = {'s00138-2817-11-01-13-39_ABP.txt', ...
         's00214-3084-11-28-16-23_ABP.txt'};

Fs = 125;
StartHours = [1, 1];
NumberOfPulses = 20;

for k = 1:length(Patients)

    %% Load patient data
    RawData = load(fullfile(Patients{k}, Files{k}));

    Time = RawData(:,1);
    ABP = RawData(:,2);

    fprintf('\nPatient %s\n', Patients{k});
    fprintf('Recording duration: %.2f hours\n', ...
        (Time(end)-Time(1))/3600);

    %% Select the analysis period
    % Select a short interval around the requested hour
    StartSample = find(Time >= Time(1) + ...
        StartHours(k)*3600, 1, 'first');

    SegmentLength = 60*Fs;
    EndSample = min(StartSample + SegmentLength - 1, length(ABP));

    ABPSegment = ABP(StartSample:EndSample);

    %% Detect heartbeat onsets
    OnsetTimes = wabp(ABPSegment);

    %% Extract ABP waveform features
    Features = abpfeature(ABPSegment, OnsetTimes);

    %% Evaluate signal quality
    [BeatQ, QualityRatio] = jSQI(Features, OnsetTimes, ABPSegment);

    fprintf('Fraction of good beats: %.3f\n', QualityRatio);

    %% Select the first 20 complete pulses
    FirstBeat = 1;

    BeatIndices = FirstBeat:(FirstBeat + NumberOfPulses - 1);

    SelectedOnsets = OnsetTimes( ...
        FirstBeat:FirstBeat + NumberOfPulses);

    SelectedFeatures = Features(BeatIndices,:);

    SelectedQuality = BeatQ(BeatIndices,1);

    SegmentStart = SelectedOnsets(1);
    SegmentEnd = SelectedOnsets(end);

    ABPPlot = ABPSegment(SegmentStart:SegmentEnd);

    SegmentTime = (0:length(ABPPlot)-1)/Fs;

    %% Extract feature locations
    SystolicTimes = SelectedFeatures(:,1) - SegmentStart + 1;
    DiastolicTimes = SelectedFeatures(:,3) - SegmentStart + 1;

    EndSystole1 = SelectedFeatures(:,9) - SegmentStart + 1;
    EndSystole2 = SelectedFeatures(:,11) - SegmentStart + 1;

    OnsetRelative = SelectedOnsets(1:end-1) ...
        - SegmentStart + 1;

    ValidSystolic = SystolicTimes >= 1 & ...
        SystolicTimes <= length(ABPPlot);

    ValidDiastolic = DiastolicTimes >= 1 & ...
        DiastolicTimes <= length(ABPPlot);

    ValidEnd1 = EndSystole1 >= 1 & ...
        EndSystole1 <= length(ABPPlot);

    ValidEnd2 = EndSystole2 >= 1 & ...
        EndSystole2 <= length(ABPPlot);

    %% Plot ABP waveform and features
    figure;

    plot(SegmentTime, ABPPlot, 'b-', 'LineWidth', 1.2);
    hold on;

    % Pulse onsets
    plot(SegmentTime(OnsetRelative), ABPPlot(OnsetRelative), ...
        'k*', 'MarkerSize', 7);

    % Systolic peaks
    plot(SegmentTime(SystolicTimes(ValidSystolic)), ...
        SelectedFeatures(ValidSystolic,2), ...
        'r^', 'MarkerSize', 7, 'MarkerFaceColor', 'r');

    % Diastolic points
    plot(SegmentTime(DiastolicTimes(ValidDiastolic)), ...
        SelectedFeatures(ValidDiastolic,4), ...
        'bv', 'MarkerSize', 7, 'MarkerFaceColor', 'b');

    % End systole using 0.3*sqrt(RR)
    plot(SegmentTime(EndSystole1(ValidEnd1)), ...
        ABPPlot(EndSystole1(ValidEnd1)), ...
        'mx', 'MarkerSize', 9, 'LineWidth', 1.5);

    % End systole using minimum slope
    plot(SegmentTime(EndSystole2(ValidEnd2)), ...
        ABPPlot(EndSystole2(ValidEnd2)), ...
        'go', 'MarkerSize', 7, 'LineWidth', 1.5);

    %% Format plot
    xlabel('Time (seconds)');
    ylabel('Arterial Blood Pressure (mmHg)');

    title(sprintf( ...
        'Patient %s: First 20 ABP Pulses Starting at %d Hour', ...
        Patients{k}, StartHours(k)));

    legend('ABP waveform', ...
        'Pulse onset (*)', ...
        'Systolic peak (^)', ...
        'Diastolic point (v)', ...
        'End systole: 0.3 sqrt(RR) (X)', ...
        'End systole: minimum slope (o)', ...
        'Location', 'northeastoutside');

    grid on;
    box on;
    hold off;

    %% Save figure
    OutputName = sprintf('ABP_patient%s_%dhours.png', ...
        Patients{k}(2:end), StartHours(k));

    saveas(gcf, OutputName);

    fprintf('Number of pulses plotted: %d\n', NumberOfPulses);
    fprintf('Good pulses: %d\n', sum(~SelectedQuality));
    fprintf('Bad pulses: %d\n', sum(SelectedQuality));
    fprintf('Figure saved as: %s\n', OutputName);

end

%% Q2 Results and Observations
% Two additional patients, s00138 and s00214, were analyzed.
% For each patient, 20 consecutive ABP pulses were plotted starting at 1 hour after the beginning of the recording.
%
% Both patients showed relatively regular ABP waveforms.
% Patient s00138 had systolic pressures of approximately 85-105 mmHg and diastolic pressures of 55-65 mmHg.
% Patient s00214 had higher systolic pressures of approximately 115-130 mmHg and diastolic pressures of 53-58 mmHg.
