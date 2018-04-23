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
%  duration.
%
%  baseline : a scalar value that indicated how many time point prior to
%  stimulus onset should be treated as the baseline window. The mean of the
%  baseline window will be subtracted from the trial activation.
%
%  electrodes : If a list of electrode labels are provided, then only these
%  electrodes will be returned in the output.
%
% Output :
%  Each field in the output structure contains a session x stimulus x time
%  array. Trials are ordered by stimulus index, not order of presentation
%  (so trial 1 in each session can be interpretted as the same stimulus.)
%
% Chris Cox 01 March 2018

    p = inputParser();
    addRequired(p, 'EEG', @isstruct);
    addRequired(p, 'TrialInfo', @istable);
    addRequired(p, 'window', @isvector);
    addRequired(p, 'baseline', @isscalar);
    addOptional(p, 'electrodes', cellstr(EEG.DIM(2).label), @iscellstr);
    addParameter(p, 'ReturnBaseline', true, @islogical);
    parse(p, EEG, TrialInfo, window, baseline, varargin{:});
    
    allelectrodes = cellstr(EEG.DIM(2).label);
    sessions = sort(unique(TrialInfo.Session));
    E = struct('label',p.Results.electrodes,'data',[]);
    for k = 1:numel(p.Results.electrodes);
        currentElectrodeFilter = strcmp(p.Results.electrodes{k}, allelectrodes);
        X = zeros(4, max(TrialInfo.Trial), window(2));
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
                    a = onset - baseline;
                    b = onset - 1;
                    baseline_mean = mean(EEG.DATA(a:b,currentElectrodeFilter));
                    B(session,stimid,:) = EEG.DATA(a:b,currentElectrodeFilter);
                    B_raw(session,stimid,:) = EEG.DATA(a:b,currentElectrodeFilter);
                else
                    baseline_mean = 0;
                end
                a = onset + window(1);
                b = a + window(2) - 1;
%                 a = onset - 200;
%                 b = onset + 1000;
                X(session,stimid,:) = EEG.DATA(a:b,currentElectrodeFilter) - baseline_mean;
            end
        end
        E(k).data = X;
        if p.Results.ReturnBaseline
            E(k).baseline = B;
            E(k).baseline_raw = B_raw;
        end
    end
end

