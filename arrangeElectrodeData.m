function [ M ] = arrangeElectrodeData( X, onsetIndex, window )
%ARRANGEELECTRODEDATA Process and subset time-by-electrode matrix
%          X : t by f matrix, t time points by f features.
% onsetIndex : s by 1 cell array, where each cell is a session. A session
%              contains a n by 1 matrix of indexes into the rows of X. Each
%              index corresponds to a stimulus onset. Indexes are expected
%              to be absolute (i.e., from beginning of experiment, not
%              beginning of session)
%     window : [start, size], where start is "number of ticks post stimulus
%              onset" and size is "length of the window in ticks", where each
%              row in X corresponds to a tick.
%
% The purpose of this function is to break up a continuous time series into
% trials within experimental sessions. The window defines a block of time
% relative to each stimulus onset to treat as a trial.  Trials will likely
% have been completed over multiple sessions, and so trial onset indexes
% should be grouped into cells by session.
%
% Example:
% % 100 time points, at 4 electrodes
% % 10 trials, over 2 sessions (i.e., 5 trials per session).
% % A trial is in the arranged output will contain time points 2--8 (start
% % at time point 2, include 6)
% X = randn(100,4); 
% onsetIndex = {[1,11,21,31,41], [51,61,71,81,91]};
% window = [2,6]; 
% arrangeElectrodeData(X, onsetIndex, window)
%
% ans = 
% 
%     [5x6 double]    [5x6 double]    [5x6 double]    [5x6 double]
%     [5x6 double]    [5x6 double]    [5x6 double]    [5x6 double]

    if nargin < 2
        window_start = 1;
        window_end = inf;
    else
        window_start = window(1);
        window_size = window(2);
        window_end = window_start + window_size;
    end

    nSessions = numel(onsetIndex);
    session = cell(1,nSessions);

    iti = cell(1,nSessions);
    iti_session_max = zeros(1,nSessions);
    for iSession = 1:nSessions;
        session{iSession} = onsetIndex{iSession}(:);
        d = diff(session{iSession});
        iti{iSession} = [d;0];
        iti_session_max(iSession) = max(d);
    end

    iti_max = max(iti_session_max);
    if window_end < iti_max
        iti_max = window_end;
    end

    for iSession = 1:nSessions
        iti{iSession}(end) = iti_max;
    end
    nTicks = iti_max - window_start;

    nElectrodes = size(X, 2);

    M = cell(nSessions,nElectrodes);
    for iElectrode = 1:nElectrodes
        for iSession = 1:nSessions
            nStimuliInSession = numel(session{iSession});
            nTrials = numel(iti{iSession});
            m = nan(nTrials,nTicks);
            for iStimulus = 1:nStimuliInSession
                t = session{iSession}(iStimulus) + window_start;
                if window_end < iti{iSession}(iStimulus);
                    o = window_end;
                else
                    o = iti{iSession}(iStimulus);
                end
                if window_start > 0
                    p = window_size;
                else
                    p = o;
                end
                m(iStimulus,1:o-window_start) = X(t:t+(p-1),iElectrode);
            end
            M{iSession,iElectrode} = m;
        end
    end
end
