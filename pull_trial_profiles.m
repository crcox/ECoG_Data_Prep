function [ E ] = pull_trial_profiles(EEG, TrialInfo, window, baseline, varargin)
% PULL_TRIAL_PROFILES Returns a structure with a field for each electrode.
%
% Input:
%  EEG : An EEG formated structure, as specified in  the Clinical
%  Neurophysiology Data Analysis Tools on MATLAB (MMscripts20141207/eeg).
%  Credit to Masao Matsuhashi, HBRC, Kyoto Univ., HMCS, NINDS, Kyoto Inst.
%  of Technology.
%
%  TrialInfo : A table object, with fields Session, Trial, ItemIndex, and
%  OnsetIndex. OnsetIndex is taken from the "tag" variables saved along with
%  the EEG object on disk. ItemIndex refers to one of the 100 images in the
%  experiment, which are repeated in a different order in each session.
%
%  window : a 2 element vector which has the window onset and window
%  duration in milliseconds.
%
%  baseline : a scalar value that indicates how many milliseconds prior to
%  stimulus onset should be treated as the baseline window. The mean of the
%  baseline window will be subtracted from the trial activation.
%
%  boxcar : a scalar value that indicates a width (in milliseconds) of a
%  boxcar defined to downsample without convolution. The timeseries for a
%  trial will be broken into N contiguous pieces of size "boxcar", and the
%  data within each boxcar is then averaged over time.
%
%  electrodes : If a list of electrode labels are provided, then only these
%  electrodes will be returned in the output.
%
% Output :
%  Each field in the output structure contains a session x stimulus x time
%  array. Trials are ordered by stimulus index, not order of presentation
%  (so trial 1 in each session can be interpretted as the same stimulus.)
%
% Chris Cox 25 March 2018

    p = inputParser();
    addRequired(p, 'EEG', @isstruct);
    addRequired(p, 'TrialInfo', @istable);
    addRequired(p, 'window', @isvector);
    addRequired(p, 'baseline', @isscalar);
    addOptional(p, 'boxcar', 0, @isscalar);
    addOptional(p, 'electrodes', cellstr(EEG.DIM(2).label), @iscellstr);
    addParameter(p, 'ReturnBaseline', false, @islogical);
    parse(p, EEG, TrialInfo, window, baseline, varargin{:});
    
    Hz = 1 / EEG.DIM(1).interval; % ticks per second
    boxcar_size   = (p.Results.boxcar / 1000) * Hz;
    window_start  = (window(1) / 1000) * Hz;
    window_size   = (window(2) / 1000) * Hz; % in ticks (where a tick is a single time-step).
    baseline_size = (baseline  / 1000) * Hz; % in ticks (where a tick is a single time-step).
    
    allelectrodes = cellstr(EEG.DIM(2).label);
    sessions = sort(unique(TrialInfo.Session));
    E = struct('label',p.Results.electrodes,'data',[]);
    for k = 1:numel(p.Results.electrodes);
        currentElectrodeFilter = strcmp(p.Results.electrodes{k}, allelectrodes);
        if ~any(currentElectrodeFilter)
            continue;
        end
        X = zeros(4, max(TrialInfo.Trial), window_size);
        B = zeros(4, max(TrialInfo.Trial), baseline);
        B_raw = zeros(4, max(TrialInfo.Trial), baseline);
        for i = 1:numel(sessions)
            session = sessions(i);
            z = TrialInfo.Session == session;
            SessionInfo = TrialInfo(z,:);
            for j = 1:size(SessionInfo,1);
                z = SessionInfo.Trial == j;
                onset = SessionInfo.OnsetIndex(z);
                stimid = SessionInfo.ItemIndex(z);
                if baseline > 0;
                    a = onset - baseline_size;
                    b = onset - 1;
                    baseline_mean = mean(EEG.DATA(a:b,currentElectrodeFilter));
                    B(session,stimid,:) = EEG.DATA(a:b,currentElectrodeFilter) - baseline_mean;
                    B_raw(session,stimid,:) = EEG.DATA(a:b,currentElectrodeFilter);
                else
                    baseline_mean = 0;
                end
                a = onset + window_start;
                b = a + window_size - 1;
                X(session,stimid,:) = EEG.DATA(a:b,currentElectrodeFilter) - baseline_mean;
            end
        end
        E(k).data = boxcar_average_timepoints(X, boxcar_size);
        if p.Results.ReturnBaseline
            E(k).baseline = boxcar_average_timepoints(B, boxcar_size);
            E(k).baseline_raw = boxcar_average_timepoints(B_raw, boxcar_size);
        end
    end
end

function Xbc = boxcar_average_timepoints(X, boxcar_size)
% BOXCAR_AVERAGE_TIMEPOINTS Downsample without convolution
% 
% Algorithm :
% 1. Trim the window so that it is evenly divisible by the boxcar.
% 2. Transpose the data to make it time by item.
% 3. Reshape to be boxcar_size by adj_window_size / boxcar_size by item.
% 4. Take mean over first dimension (default of mean()).
% 5. Squeeze out the now singleton first dimension
% 6. Transpose back to item by time.
    if boxcar_size == 0
        Xbc = X;
    else
        [nsessions, nitems, window_size] = size(X);
        a = 1;
        b = window_size - rem(window_size, boxcar_size);
        c = b / boxcar_size;
        Xbc = nan(nsessions, nitems, c);
        if c < window_size
            for i = 1:nsessions
                x = squeeze(X(i,:,a:b));
                r = size(x, 1);
                Xbc(i,:,:) = squeeze(mean(reshape(x',boxcar_size,c,r),1))';
            end
        else
            Xbc = X;
        end
    end
end