function [Pt] = load_and_process(subjects, datacode, datadir)
% LOAD_AND_PROCESS Load data from the Kyoto Naming ECoG studies
%
    key = initialize_datakey('raw');
    subjects = subject_id_crossref(subjects,'ptcode');
    NSUBJ = numel(subjects);

    for iSubject=1:NSUBJ
        sdir = subjects{iSubject}; % PtCode is also directory name for subject
        F = selectbyfield(key, 'dirname', sdir);
        if any(cellfun('isempty', struct2cell(F)))
            fprintf('Skipping subject %d, %s because of missing data.\n',subjects(iSubject),datacode);
            continue
        else
            fprintf('Beginning subject %s, %s.\n',subjects{iSubject},datacode);
        end
        spath = fullfile(datadir,datacode,F.dirname,F.filename);
        fprintf('Loading %s...\n', spath);
        Pt = load(spath);

        DataIsSplit = numel(F.variables)>1;
        if DataIsSplit
            Pt = reassemble_split_data(Pt, F);
        else
            cvar = F.variables{1};
            Pt.LFP = Pt.(cvar);
            Pt = rmfield(Pt,cvar);
        end

        Pt.tags = cell(numel(F.sessions), 1);
        for i = 1:numel(cell2mat(F.sessions))
            tagfield = sprintf(F.sessiontag, i);
            Pt.tags{i} = Pt.(tagfield)(:);
            Pt = rmfield(Pt,tagfield);
        end
    end
end

function key = initialize_datakey(datacode)
% INITIALIZE_DATAKEY A structure that encodes variable names for common data
% elements subjects.
%
%   In the source data, naming conventions are not consistent. The
%   following structures explicitly represent all the file and field names that
%   need to be referenced.
%
    subjects = 1:10;
    ptcodes = subject_id_crossref(subjects,'ptcode');

    key = struct('subject',num2cell(subjects), 'dirname', ptcodes, 'filename',[], 'variables', [], 'sessions', [], 'sessiontag', []);
    switch datacode
        % This bit is necessary because the raw data were not stored with a
        % consistent naming convention.
        case 'raw'
            % Subject 1
            key(1).filename = 'namingERP_Pt01.mat';
            key(1).variables = {'namingERP_data_PtYK_Pt01'};
            key(1).sessions = {1:4};
            key(1).sessiontag = 'tag_ss%02d_all';
            % Subject 2
            key(2).filename = 'namingERP_Pt02.mat';
            key(2).variables = {'namingERP_data_Pt02'};
            key(2).sessions = {1:4};
            key(2).sessiontag = 'Tag_ss%02d_all';
            % Subject 3
            key(3).filename = 'namingERP_Pt03.mat';
            key(3).variables = {'namingERP_data_Pt03'};
            key(3).sessions = {1:4};
            key(3).sessiontag = 'Tag_ss%02d_all';
            % Subject 4
            key(4).filename = [];
            key(4).variables = [];
            key(4).sessions = [];
            key(4).sessiontag = [];
            % Subject 5
            key(5).filename = 'namingERP_Pt05.mat';
            key(5).variables = {'namingERP_data_Pt05'};
            key(5).sessions = {1:4};
            key(5).sessiontag = 'tag_ss%02d';
            % Subject 6
            key(6).filename = 'namingERP_Pt06.mat';
            key(6).variables = {'namingERP_data_Pt06'};
            key(6).sessions = {1:4};
            key(6).sessiontag = [];
            % Subject 7
            key(7).filename = 'namingERP_Pt07.mat';
            key(7).variables = {'namingERPdataPt07','namingERPdataPt07_ss0304'};
            key(7).sessions = {1:2,3:4};
            key(7).sessiontag = 'Tag_ss%02d';
            % Subject 8
            key(8).filename = 'namingERP_Pt08.mat';
            key(8).variables = {'namingERPdataPt08'};
            key(8).sessions = {1:4};
            key(8).sessiontag = 'tag_%02d';
            % Subject 9
            key(9).filename = 'namingERP_Pt09.mat';
            key(9).variables = {'namingERPdataPt09'};
            key(9).sessions = {1:4};
            key(9).sessiontag = 'tag%02d';
            % Subject 10
            key(10).filename = 'namingERP_Pt10.mat';
            key(10).variables = {'namingERPdataPt10'};
            key(10).sessions = {1:4};
            key(10).sessiontag = 'tagall%02d';
        case 'ref'
            error('Nothing is implemented for referenced data yet!');
        %     key_ref = {...
        %       'namingERP_Pt01_refD14.mat',{'namingERP_data_PtYK_Pt01_refD14'},{[1,2,3,4]},'tag_ss%02d_all';...
        %       [],[],[],[];...
        %       [],[],[],[];...
        %       'namingERP_PtMA_REF4.mat',{'namingERP_data_ss01ss02_REF4','namingERP_data_ss03ss04_REF4'},{[1,2],[3,4]},'tag_ss%02d';...
        %       'namingERP_Pt05_ref02.mat',{'namingERP_data_Pt05_ref02'},{[1,2,3,4]},'tag_ss%02d';...
        %       'namingERP_Pt06_ref01.mat',{'namingERP_data_Pt06_ref01'},{[1,2,3,4]},[];...
        %       'namingERPdata_Pt07_ref02.mat',{'namingERPdataPt07_ss0102_ref02','namingERPdataPt07_ss0304_ref02'},{[1,2],[3,4]},'Tag_ss%02d';...
        %       'namingERPdata_Pt08_ref03.mat',{'namingERPdataPt08_ref03'},{[1,2,3,4]},'tag_%02d';...
        %       [],[],[],[];...
        %       'namingERPdata_Pt10_ref02.mat',{'namingERPdataPt10_ref02'},{[1,2,3,4]},'tagall%02d';...
        %     };
    end

end

function s_out = subject_id_crossref(s_in, outtype)
% SUBJECT_ID_CROSSREF Translate subject number to PtCode and vice verse
%
%  By default, the function toggles between subject number (numeric type) and a
%  PtCode string ('Pt%02d'). If the second (optional) argument is provided, it
%  constrains the output to be the specified type (either 'number' or 'ptcode').
%  This can be used in functions to ensure that subject identifiers
%  are in a specific format when receiving user input.
%
%   Usage:
%   s_in = 1:3;
%   s_out = subject_id_crossref(s_in);
%   % s_out = {'Pt01','Pt02','Pt03'}
%
%   s_out = subject_id_crossref(s_out);
%   % s_out = [1  2  3]
%
%   s_out = subject_id_crossref(s_out,'number');
%   % s_out = [1  2  3]
%
%   s_out = subject_id_crossref(s_out,'toggle');
%   % s_out = {'Pt01','Pt02','Pt03'}
%
%   s_out = subject_id_crossref(s_out,'idcode');
%   % s_out = {'Pt01','Pt02','Pt03'}
%
%   s_out = subject_id_crossref(s_out,'number');
%   % s_out = [1  2  3]
%
%   s_out = subject_id_crossref(s_out,'idcode');
%   % s_out = {'Pt01','Pt02','Pt03'}
%
    if nargin < 2
        outtype = 'toggle';
    end
    if ~ischar(outtype) || ~any(strcmpi(outtype, {'numeric','ptcode','toggle'}))
        error('Invalid outtype argument. Must be one of ''numeric'', ''ptcode'', or ''toggle''.');
    end
    if isnumeric(s_in)
        if strcmpi(outtype,'numeric');
            s_out = s_in;
        else
            s_out = arrayfun(@(x) sprintf('Pt%02d', x), s_in, 'UniformOutput', false);
        end
    else
        if strcmpi(outtype,'ptcode');
            s_out = s_in;
        else
            s_out = sscanf(s_in, 'Pt%02d');
        end
    end
end

function Pt = reassemble_split_data(Pt, F)
% REASSEMBLE_SPLIT_DATA Enforce more consistent data representation across
% subjects by combining data split over multiple sub-structures.
%
    nChunks = numel(F.variables);
    nTicks = 0;
    for iChunk = 1:nChunks
        cvar = F.variables{iChunk};
        nTicks = nTicks + size(Pt.(cvar).DATA,1);
    end
    interval = Pt.(cvar).DIM(1).interval;
    electrodeLabels = Pt.(cvar).DIM(2).label;
    Pt.LFP(1) = init_source_struct(nTicks,electrodeLabels,interval);
    psize = 0;
    pscale = 0;
    for iChunk = 1:nChunks
        cvar = F.variables{iChunk};
        sessions = F.sessions{iChunk};
        tagfmt = F.sessiontag;
        a = psize + 1;
        b = psize + size(Pt.(cvar).DATA,1);
        Pt.LFP.DATA(a:b,:) = Pt.(cvar).DATA;
        Pt.LFP.DIM(1).scale(a:b) = Pt.(cvar).DIM(1).scale + pscale;

        for iSession = sessions
            tag = sprintf(tagfmt,iSession);
            Pt.(tag) = Pt.(tag) + psize;
        end

        psize = b;
        pscale = max(Pt.(cvar).DIM(1).scale);
        Pt = rmfield(Pt,cvar);
    end
end
