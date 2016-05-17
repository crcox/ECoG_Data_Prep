function [ M ] = arrangeElectrodeData( X, onsetIndex, window )
%ARRANGEELECTRODEDATA Process and subset time-by-electrode matrix
%          X : t by f matrix, t time points by f features.
% onsetIndex : s by 1 cell array, where each cell is a session. A session
%              contains a n by 1 matrix of indexes into the rows of X. Each
%              index corresponds to a stimulus onset.
%     window : [start, size], where start is "number of ticks post stimulus
%              onset" and size is "length of the window in ticks", where each
%              row in X corresponds to a tick.

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
        for iSession = 1:4
            nTrials = numel(iti{iSession});
            m = nan(nTrials,nTicks);
            for iStimulus = 1:100
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
