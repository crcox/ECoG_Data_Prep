function [ E ] = pull_trial_profiles_derivatives(EEG, TrialInfo, window, baseline, varargin)
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
    addOptional(p, 'target_resolution', 0, @isscalar);
    addOptional(p, 'slope_interval', 0, @isscalar);
    addOptional(p, 'electrodes', cellstr(EEG.DIM(2).label), @iscellstr);
    parse(p, EEG, TrialInfo, window, baseline, varargin{:});
    
    % All times are supplied in ms and converted to ticks.
    Hz = 1 / EEG.DIM(1).interval; % ticks per second
    target_resolution   = (p.Results.target_resolution / 1000) * Hz;
    slope_interval   = (p.Results.slope_interval / 1000) * Hz;
    window_start  = (window(1) / 1000) * Hz;
    window_size   = (window(2) / 1000) * Hz; % in ticks (where a tick is a single time-step).
    
    Electrodes = p.Results.electrodes;
    
    allelectrodes = cellstr(EEG.DIM(2).label);
    sessions = sort(unique(TrialInfo.Session));
    E = struct('label',Electrodes,'data',[]);
    for k = 1:numel(Electrodes)
        currentElectrodeFilter = strcmp(Electrodes{k}, allelectrodes);
        if ~any(currentElectrodeFilter)
            continue;
        end
        p = 0;
        q = 0;
        if slope_interval > 0
            p = slope_interval - target_resolution;
            p = p + mod(p,2);
            q = p / 2;
        end
        X = zeros(4, max(TrialInfo.Trial), window_size + p);
        for i = 1:numel(sessions)
            session = sessions(i);
            z = TrialInfo.Session == session;
            SessionInfo = TrialInfo(z,:);
            for j = 1:size(SessionInfo,1)
                z = SessionInfo.Trial == j;
                onset = SessionInfo.OnsetIndex(z);
                stimid = SessionInfo.ItemIndex(z);
                a = onset + window_start;
                b = (a + window_size) - 1;
                a = a - q;
                b = b + q;
                X(session,stimid,:) = EEG.DATA(a:b,currentElectrodeFilter);
            end
        end
        E(k).data = derivative_of_timepoints(X, target_resolution, slope_interval, p);
    end
end

function Xdv = derivative_of_timepoints(X, target_resolution, slope_interval, p)
% BOXCAR_AVERAGE_TIMEPOINTS Downsample without convolution

    if slope_interval == 0
        Xdv = X;
    else
        q = p / 2;
        [nsessions, nitems, window_size_p] = size(X);
        window_size = window_size_p - p;
        a = 1;
        b = window_size - rem(window_size, target_resolution);
        c = b / target_resolution;
        Xdv = nan(nsessions, nitems, c);
        % Columns of M correspond to windows of time that the slope will be
        % computed over. The temporal resolution of the results will be
        % equal to the number of columns in M, which will be c.
        offset = target_resolution .* 0:(c-1); %#ok<BDSCA>
        offset = offset + floor(target_resolution / 2);
        tmp = (1:slope_interval) - floor(slope_interval/2);
        M = bsxfun(@plus, tmp', offset) + q;
        if c < window_size
            for i = 1:nsessions
                % Time by Item after transpose
                x = squeeze(X(i,:,a:(b+q)))';
                % Number of items
                r = size(x, 2);
                for j = 1:r
                    yy = reshape(x(M,j), size(M));
                    xx = 1:slope_interval;
                    Xdv(i,j,:) = best_fit_slope(yy,xx);
                end
            end
        else
            Xdv = X;
        end
    end
end

function [m] = best_fit_slope(y, x)
% BEST_FIT_slope Least-squares estimate of best fit slope.
    if nargin < 2
        x = linspace(-0.5,0.5,size(y,1));
    end
    x = (x - mean(x)) / 100;
    ym = mean(y);
    yc = bsxfun(@minus, y, ym);
    m = sum(bsxfun(@times, yc, x(:))) ./ sum(x .^ 2);
    % The intercept, which we intentionally do not return, is:
%     b = ym;
end
    