function [X, t, tau] = parlikarCO(ABP, onset, feat, beatq, Fs, windowSize)
% Uncalibrated Parlikar cardiac output estimator.
%
% ABP - arterial blood pressure waveform (mmHg)
% onset - pulse onset indices from wabp
% feat - features from abpfeature
% beatq - beat quality output from jSQI
% Fs - sampling frequency (Hz)
% windowSize - odd moving window size, e.g. 9 beats
%
% X - uncalibrated CO-related estimate
% t - time in minutes from start of ABP segment
% tau - estimated Windkessel time constant (seconds)

if nargin < 5 || isempty(Fs)
    Fs = 125;
end

if nargin < 6 || isempty(windowSize)
    windowSize = 9;
end

if mod(windowSize, 2) == 0
    error('windowSize must be odd.');
end

ABP = ABP(:);
onset = onset(:);

% Each row corresponds to an interval between two onsets
n = min([length(onset)-1, size(feat,1), size(beatq,1)]);

if n < windowSize
    error('Not enough beats for the selected moving window.');
end

onset = onset(1:n+1);
feat = feat(1:n,:);
beatq = beatq(1:n,:);

if any(onset < 1) || any(onset > length(ABP))
    error('Onset indices are outside the ABP signal.');
end

% Beat-level pressure features
DAP = feat(:,4);
Pmean = feat(:,6);

% Beat interval in sec
T = diff(onset)/Fs;

% Pressure change between onset samples
P_onset = ABP(onset);
DeltaP = diff(P_onset);

% Parlikar equation, triangular pulse approximation
PP = 2*(Pmean-DAP);

% Identify good-qualtiy beats
goodBeat = isfinite(beatq(:, 1)) & beatq(:, 1) == 0;

% Exclude beat containing invalid ABP samples
for i = 1:n
    beatSignal = ABP(onset(i):onset(i+1) - 1);

    if any(~isfinite(beatSignal)) || any(beatSignal <= 0) || any(beatSignal >= 300)
        goodBeat(i) = false;
    end
end

% Estimate 1/tau using moving-window least squares
invTau = nan(n,1);
halfWindow = floor(windowSize/2);

for i = 1:n
    idx = max(1,i-halfWindow):min(n,i+halfWindow);

    A = Pmean(idx);
    b = (PP(idx)-DeltaP(idx))./T(idx);

    valid = isfinite(A) & isfinite(b) & isfinite(T(idx)) & T(idx) > 0 & A > 0;

    if sum(valid) >= 3
        estimate = A(valid)\b(valid);

        if isfinite(estimate) && estimate > 0
            invTau(i) = estimate;
        end
    end
end

tau = 1./invTau;

% Parlikar equation
X = DeltaP./T + Pmean./tau;

% Time at the beginning of each beat, in min
t = (onset(1:n)-onset(1))/(Fs*60);

% Remove beats marked as poor quality by jSQI
bad = ~goodBeat;
X(bad) = NaN;
tau(bad) = NaN;

% Reject invalid or non-positive estimates
X(~isfinite(X) | X <= 0) = NaN;
tau(~isfinite(tau) | tau <= 0) = NaN;

end

%[appendix]{"version":"1.0"}
%---
