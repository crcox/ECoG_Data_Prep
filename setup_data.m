function setup_data(varargin)
% DEPENDENCIES
% 1. glmnet
% 2. sqrt_truncate_r.m (from WholeBrain_RSA/src)
    if nargin == 1
        % If there is only one argument, treat it as a path to a parameter
        % file, or a path to a directory that contains a file named
        % 'params.json'.
        [p,f,e] = fileparts(varargin{1});
        if strcmp('.json', e);
            jsonpath = varargin{1};
        elseif isempty(e)
            jsonpath = fullfile(p,f,'params.json');
        else
            warning('Expected a json file. %s does not have a conventional file extension for a json file (.json).');
            jsonpath = varargin{1};
        end
    elseif nargin == 0
        % If there are no arguments, look for a parameter file in the
        % current directory.
        jsonpath = './params.json';
    else
        % If more than 1 parameter is provided, then assume that all
        % arguments are provided at the command line and parse as such.
        jsonpath = [];
    end

    % If nargin < 2, jsonpath will be set to something. It may not be a
    % valid file, so 'jsonload()' may result in a fatal error.
    if isempty(jsonpath);
        isjsoncfg = false;
    else
        jdat = loadjson(jsonpath);
        isjsoncfg = true;
    end

    % Input is handled differently depending on how it is input and whether
    % the program is being run from a Unix/Windows terminal or within
    % Matlab propper. isdeployed() checks whether the program is being run
    % after compilation with mcc.
    p = inputParser();
    if isdeployed && ~isjsoncfg
        % In this case, all arguments must be handled as strings, since
        % they will be passed directly from the Unix/Windows terminal.
        addParameter(p, 'onset',[]);
        addParameter(p, 'duration',[]);
        addParameter(p, 'baseline','0');
        addParameter(p, 'subjects', '1 2 3 5 7 8 9 10', @ischar);
        addParameter(p, 'boxcar', '0', @ischar);
        addParameter(p, 'average', '1', @ischar);
        addParameter(p, 'datacode', 'raw', @ischar);
        addParameter(p, 'dataroot', '/mnt/sw01-home01/mbmhscc4/scratch/data/Naming_ECoG');
        addParameter(p, 'metaroot', '/mnt/sw01-home01/mbmhscc4/scratch/data/Naming_ECoG/meta');
        addParameter(p, 'datarootout', []);
        addParameter(p, 'cvpath', []);
        addParameter(p, 'overwrite', '0');
        addParameter(p, 'WriteIndividualMetadata', '0');
        parse(p, varargin{:});
        if isempty(p.Results.onset)
            error('A window onset time must be provided.');
        end
        if isempty(p.Results.duration)
            error('A window duration must be provided.');
        end
        fprintf('        Window onset (ms): %s\n', p.Results.onset);
        fprintf('     Window duration (ms): %s\n', p.Results.duration);
        fprintf('                 subjects: %s\n', p.Results.subjects);
        fprintf('              Boxcar Size: %s\n', p.Results.boxcar);
        fprintf('    Average over sessions: %s\n', p.Results.average);
        fprintf('    Data root (for input): %s\n', p.Results.dataroot);
        fprintf('    Meta root (for input): %s\n', p.Results.metaroot);
        fprintf('  Overwrite existing data: %s\n', p.Results.overwrite);
        fprintf('Write individual metadata: %s\n', p.Results.WriteIndividualMetadata);

        %% Setup
        WindowStartInMilliseconds = str2double(p.Results.onset);
        WindowSizeInMilliseconds = str2double(p.Results.duration);
        BaselineSizeInMilliseconds = str2double(p.Results.baseline);
        BoxCarSize = str2double(p.Results.boxcar);
        AverageOverSessions = str2double(p.Results.average);
        OVERWRITE = str2double(p.Results.overwrite);
        WriteIndividualMetadata = str2double(p.Results.WriteIndividualMetadata);
        clean_string = regexprep(p.Results.subjects, ',? *', ' ');
        SUBJECTS = str2double(strsplit(clean_string));
    else
        % In this case, arguments can be handled in their expected types
        % (numbers as numbers, etc). The input may be coming directly from
        % the Matlab command line or script, or parsed from a json file.
        addParameter(p, 'onset',[]);
        addParameter(p, 'duration',[]);
        addParameter(p, 'baseline',0);
        addParameter(p, 'subjects', [1:3,5,7:10], @isnumeric);
        addParameter(p, 'boxcar', 0);
        addParameter(p, 'average', 1);
        addParameter(p, 'datacode', 'raw', @ischar);
        addParameter(p, 'dataroot', 'D:\ECoG\KyotoNaming\data');
        addParameter(p, 'metaroot', 'C:\Users\mbmhscc4\MATLAB\ECOG\naming\data');
        addParameter(p, 'datarootout', []);
        addParameter(p, 'cvpath', []);
        addParameter(p, 'overwrite', 0);
        addParameter(p, 'WriteIndividualMetadata', 0);

        if isjsoncfg
            % Input read from a json file will be represented as a
            % structure. The inputParser cannot handle this structure, and
            % so we need to get the structure into an amenable format.
            varargin = [fieldnames(jdat); struct2cell(jdat)];
        end
        parse(p, varargin{:});
        fprintf('        Window onset (ms): %d\n', p.Results.onset);
        fprintf('     Window duration (ms): %d\n', p.Results.duration);
        fprintf('                 subjects: %s\n', strjoin(strsplit(num2str(p.Results.subjects)), ', '));
        fprintf('              Boxcar Size: %d\n', p.Results.boxcar);
        fprintf('    Average over sessions: %d\n', p.Results.average);
        fprintf('    Data root (for input): %s\n', p.Results.dataroot);
        fprintf('    Meta root (for input): %s\n', p.Results.metaroot);
        fprintf('  Overwrite existing data: %d\n', p.Results.overwrite);
        fprintf('Write individual metadata: %d\n', p.Results.WriteIndividualMetadata);

        %% Setup
        WindowStartInMilliseconds = p.Results.onset;
        WindowSizeInMilliseconds = p.Results.duration;
        BaselineSizeInMilliseconds = p.Results.baseline;
        BoxCarSize = p.Results.boxcar;
        AverageOverSessions = p.Results.average;
        OVERWRITE = p.Results.overwrite;
        WriteIndividualMetadata = p.Results.WriteIndividualMetadata;
        SUBJECTS = p.Results.subjects;
    end
    DATA_DIR = p.Results.dataroot;
    if isempty(p.Results.datarootout)
        DATA_ROOT_OUT = DATA_DIR;
    else
        DATA_ROOT_OUT = p.Results.datarootout;
    end
    META_DIR = p.Results.metaroot;
    STIM_DIR = fullfile(META_DIR,'stimuli');
    FILTER_DIR = fullfile(META_DIR,'filters');
    COORD_DIR = fullfile(META_DIR,'coords');
    SIM_DIR = fullfile(META_DIR,'similarity');
    TARGET_DIR = fullfile(META_DIR,'targets');
    CV_DIR = fullfile(META_DIR,'cv');
    cvpath = p.Results.cvpath;
    datacode = p.Results.datacode;


    %% Define Output directory
    if AverageOverSessions == 1
        base_dir = 'avg';
    else
        base_dir = 'full';
    end
    DATA_DIR_OUT = fullfile(...
        DATA_ROOT_OUT,...
        base_dir,...
        'BoxCar',sprintf('%03d',BoxCarSize),...
        'WindowStart',sprintf('%04d',WindowStartInMilliseconds),...
        'WindowSize',sprintf('%04d',WindowSizeInMilliseconds)...
    );

    if ~exist(DATA_DIR_OUT,'dir')
        mkdir(DATA_DIR_OUT)
    end

    fprintf('Looking for source ECoG data in:\n\t%s\n',fullfile(DATA_DIR,'raw'));
    fprintf('Looking for stimulus labels and trial orders in:\n\t%s\n',STIM_DIR);
    fprintf('Looking for looking for (basal) label-to-coordinate mapping in:\n\t%s\n',COORD_DIR);
    fprintf('Looking for similarity structures in:\n\t%s\n\n',SIM_DIR);

    fprintf('Average over sessions: %d\n', AverageOverSessions);
    fprintf('Box-car size for time-series averaging: %d\n', BoxCarSize);
    fprintf('Beginning of time selection window: %d ms post stim onset\n', WindowStartInMilliseconds);
    fprintf('Length of time window selection: %d ms\n', WindowSizeInMilliseconds);
    fprintf('Processed data will be written to:\n\t%s\n', DATA_DIR_OUT);
    fprintf('Subjects:\n');
    disp(SUBJECTS)

    %% Read presentation order and stim labels from file
    % All subjects have the same order
    file_stim_order = fullfile(STIM_DIR, 'picture_namingERP_order.csv');
    file_stim_key = fullfile(STIM_DIR, 'picture_namingERP_key.csv');
    
    stim_order = readtable(file_stim_order);
    stim_key = readtable(file_stim_key);
% 
%     fid = fopen(file_stim_order);
%     [~] = textscan(fid,'%s %s %s',1,'Delimiter',','); % stim_order_head
%     stim_order = textscan(fid,'%d %d %d','Delimiter',',');
%     fclose(fid);
% 
%     fid = fopen(file_stim_key);
%     [~] = textscan(fid,'%s %s',1,'Delimiter',','); % stim_key_head
%     stim_key = textscan(fid,'%d %s','Delimiter',',');
%     fclose(fid);

    nitems = numel(stim_key.Stimulus);
    nsessions = numel(unique(stim_order.Session));

    %% Ensure stimulus labels are sorted by their index
    [stim_key.ItemIndex,ix] = sort(stim_key.ItemIndex);
    stim_key.Stimulus = stim_key.Stimulus(ix);

    %% Generate sort indexes
    x = tabulate(stim_order.Session);
    stim_order_ix = mat2cell(stim_order.ItemIndex, x(:,2), 1);
    stim_sort_ix = cell(size(x,1),1);
    for i = 1:size(x,1)
        [~,stim_sort_ix{i}] = sort(stim_order_ix{i});
    end

    %% Load (basal) coordinates
    [subject,electrode,x,y,z] = importcoordinates(fullfile(COORD_DIR,'MNI_basal_electrodes_Pt01_10_w_label.csv'));
    xyz = [x,y,z]; clear x y z
    x = tabulate(subject);
    XYZ = mat2cell(xyz,cell2mat(x(:,2)),3);
    ELECTRODE = mat2cell(electrode,cell2mat(x(:,2)),1);

    %% Prep metadata structure
    metadata = struct(...
      'subject',num2cell(SUBJECTS),...
      'targets',struct(),...
      'stimuli',{stim_key.Stimulus},...
      'filters',struct(),...
      'coords',struct(),...
      'cvind',[],...
      'nrow',[],...
      'ncol',[],...
      'sessions',[],...
      'samplingrate',[],...
      'AverageOverSessions',[],...
      'BoxCarSize',[],...
      'WindowStartInMilliseconds',[],...
      'WindowSizeInMilliseconds',[]...
    );

    % TARGETS
    t_type_fmt = {'category','embedding','similarity'};
    n_tf = numel(t_type_fmt);
    for i_tf = 1:n_tf
        type_fmt = t_type_fmt{i_tf};
        switch type_fmt
            case 'category'
                t_type_categories = list_files(fullfile(TARGET_DIR,type_fmt));
                n_c = numel(t_type_categories);
                for i_c = 1:n_c
                    category_file = fullfile(TARGET_DIR,type_fmt,t_type_categories{i_c});
                    category_label = strip_extension(t_type_categories{i_c});
                    fid = fopen(category_file);
                    tmp = textscan(fid, '%s %u8','Delimiter',',');
                    category_labels = tmp{1};
                    category_targets = tmp{2};
                    fclose(fid);
                    metadata = installCategoryStructure(metadata, category_targets, category_labels, category_label, []);
                end
            case {'embedding','similarity'}
                t_type_sims = list_dirs(fullfile(TARGET_DIR,type_fmt));
                n_ts = numel(t_type_sims);
                for i_ts = 1:n_ts
                    type_sim = t_type_sims{i_ts};
                    t_sim_sources = list_dirs(fullfile(TARGET_DIR,type_fmt,type_sim));
                    n_s = numel(t_sim_sources);
                    for i_s = 1:n_s
                        source = t_sim_sources{i_s};
                        t_source_metrics = list_files(fullfile(TARGET_DIR,type_fmt,type_sim,source));
                        t_source_metrics = t_source_metrics(~strcmp('labels.txt',t_source_metrics));
                        n_m = numel(t_source_metrics);
                        for i_m = 1:n_m
                            metric_file = t_source_metrics{i_m};
                            metric_label = strip_extension(metric_file);
                            structure_file = fullfile(TARGET_DIR,type_fmt,type_sim,source,metric_file);
                            structure_label_file = fullfile(TARGET_DIR,type_fmt,type_sim,source,'labels.txt');
                            structure_matrix = csvread(structure_file);
                            fid = fopen(structure_label_file);
                            tmp = textscan(fid, '%s');
                            structure_labels = tmp{1};
                            fclose(fid);
                            metadata = installSimilarityStructure(metadata, structure_matrix, structure_labels, type_sim, source, metric_label, type_fmt);
                        end
                    end
                end
        end
    end

	animate = [ones(50,1);zeros(50,1)];
%    TARGETS = struct(...
%        'label', {'animate','semantic','semantic'},...
%        'type', {'category','similarity','embedding'},...
%        'sim_source',{[],'NEXT','NEXT'},...
%        'sim_metric',{[],'correlation','euclidean'},...
%        'target',{animate,S,embedding}...
%    );
%     if AverageOverSessions == 0;
%         for iTarget = 1:numel(TARGETS)
%             switch lower(TARGETS(iTarget).type)
%             case {'category','embedding'}
%                 TARGETS(iTarget).target = repmat(TARGETS(iTarget).target,nsessions,1);
%             case 'similarity'
%                 TARGETS(iTarget).target = repmat(TARGETS(iTarget).target,nsessions,nsessions);
%             end
%             SCHEMES = repmat(SCHEMES,nsessions,1);
%         end
%     end

    % CV
    if isempty(cvpath)
        nschemes = 10;
        nfolds = 10;
        SCHEMES = zeros(nitems, nschemes);
        for iScheme = 1:nschemes
            c = cvpartition(animate,'KFold', nfolds);
            for iFold = 1:nfolds
                SCHEMES(:,iScheme) = SCHEMES(:,iScheme) + (test(c, iFold) * iFold);
            end
        end
    else
        if exist(cvpath, 'file')
            load(cvpath, 'CV');
        else
            cvfile = cvpath;
            cvpath = fullfile(CV_DIR, cvfile);
            load(cvpath, 'CV');
        end
        SCHEMES = CV;
    end

    %% Define metadata
    NSUBJ = numel(SUBJECTS);
    for iSubject = 1:NSUBJ
        M = selectbyfield(metadata, 'subject', SUBJECTS(iSubject));
        % FILTERS
        % This is kind of a place holder. When we remove outliers, we'll want
        % to create filters for that.
        FILTERS = struct('label',[],'dimension',[],'filter',[]);
        % ---------
        if AverageOverSessions == 1;
            M.sessions = [];
            M.nrow = 100;
        else
            M.sessions = stim_order.Session;
            M.nrow = 400;
        end
        M.AverageOverSessions = AverageOverSessions;
        M.BoxCarSize = BoxCarSize;
        M.WindowStartInMilliseconds = WindowStartInMilliseconds;
        M.WindowSizeInMilliseconds = WindowSizeInMilliseconds;
        if numel(fieldnames(M.filters)) == 0
            M.filters = FILTERS;
        end
        M.cvind = SCHEMES;
        M.ncol = 0; % will be set later
        metadata = replacebyfield(metadata, M, 'subject', SUBJECTS(iSubject));
    end

    %% Load And Process Data
    % NOTE: In the source data, naming conventions are not consistent. The
    % following structures explicitly represent all the file and field names that
    % need to be referenced.
    filelist = struct('subject',num2cell(1:10), 'filename',[], 'variables', [], 'sessions', [], 'sessiontag', []);
    switch datacode
        % This bit is necessary because the raw data were not stored with a
        % consistent naming convention.
        case 'raw'
            % Subject 1
            filelist(1).filename = 'namingERP_Pt01.mat';
            filelist(1).variables = {'namingERP_data_PtYK_Pt01'};
            filelist(1).sessions = {1:4};
            filelist(1).sessiontag = 'tag_ss%02d_all';
            % Subject 2
            filelist(2).filename = 'namingERP_Pt02.mat';
            filelist(2).variables = {'namingERP_data_Pt02'};
            filelist(2).sessions = {1:4};
            filelist(2).sessiontag = 'Tag_ss%02d_all';
            % Subject 3
            filelist(3).filename = 'namingERP_Pt03.mat';
            filelist(3).variables = {'namingERP_data_Pt03'};
            filelist(3).sessions = {1:4};
            filelist(3).sessiontag = 'Tag_ss%02d_all';
            % Subject 4
            filelist(4).filename = [];
            filelist(4).variables = [];
            filelist(4).sessions = [];
            filelist(4).sessiontag = [];
            % Subject 5
            filelist(5).filename = 'namingERP_Pt05.mat';
            filelist(5).variables = {'namingERP_data_Pt05'};
            filelist(5).sessions = {1:4};
            filelist(5).sessiontag = 'tag_ss%02d';
            % Subject 6
            filelist(6).filename = 'namingERP_Pt06.mat';
            filelist(6).variables = {'namingERP_data_Pt06'};
            filelist(6).sessions = {1:4};
            filelist(6).sessiontag = [];
            % Subject 7
            filelist(7).filename = 'namingERP_Pt07.mat';
            filelist(7).variables = {'namingERPdataPt07','namingERPdataPt07_ss0304'};
            filelist(7).sessions = {1:2,3:4};
            filelist(7).sessiontag = 'Tag_ss%02d';
            % Subject 8
            filelist(8).filename = 'namingERP_Pt08.mat';
            filelist(8).variables = {'namingERPdataPt08'};
            filelist(8).sessions = {1:4};
            filelist(8).sessiontag = 'tag_%02d';
            % Subject 9
            filelist(9).filename = 'namingERP_Pt09.mat';
            filelist(9).variables = {'namingERPdataPt09'};
            filelist(9).sessions = {1:4};
            filelist(9).sessiontag = 'tag%02d';
            % Subject 10
            filelist(10).filename = 'namingERP_Pt10.mat';
            filelist(10).variables = {'namingERPdataPt10'};
            filelist(10).sessions = {1:4};
            filelist(10).sessiontag = 'tagall%02d';
        case 'ref'
        %     filelist_ref = {...
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

    for iSubject=1:NSUBJ
        sdir = sprintf('Pt%02d',SUBJECTS(iSubject));
        F = selectbyfield(filelist, 'subject', SUBJECTS(iSubject));
        if any(cellfun('isempty', struct2cell(F)))
            fprintf('Skipping subject %d, %s because of missing data.\n',SUBJECTS(iSubject),datacode);
            continue
        else
            fprintf('Beginning subject %d, %s.\n',SUBJECTS(iSubject),datacode);
        end
        dpath_out = fullfile(DATA_DIR_OUT, sprintf('s%02d_%s.mat',SUBJECTS(iSubject),datacode));
        spath = fullfile(DATA_DIR,datacode,sdir,F.filename);
        fprintf('Loading %s...\n', spath);
        Pt = load(spath);

        nChunks = numel(F.variables);
        if nChunks > 1
            nTicks = 0;
            for iChunk = 1:nChunks
                cvar = F.variables{iChunk};
                nTicks = nTicks + size(Pt.(cvar).DATA,1);
            end
            interval = Pt.(cvar).DIM(1).interval;
            electrodeLabels = Pt.(cvar).DIM(2).label;
            Pt.LFP(1) = init_source_struct(nTicks,electrodeLabels,interval);
            for iChunk = 1:nChunks
                cvar = F.variables{iChunk};
                sessions = F.sessions{iChunk};
                tagfmt = F.sessiontag;
                if iChunk == 1
                    a = 1;
                    b = size(Pt.(cvar).DATA,1);
                    Pt.LFP.DATA(a:b,:) = Pt.(cvar).DATA;
                    Pt.LFP.DIM(1).scale(a:b) = Pt.(cvar).DIM(1).scale;
                    psize = b;
                    pscale = max(Pt.(cvar).DIM(1).scale);
                    Pt = rmfield(Pt,cvar);
                else
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
        else
            cvar = F.variables{1};
            Pt.LFP = Pt.(cvar);
            Pt = rmfield(Pt,cvar);
        end

        ecoord = ELECTRODE{SUBJECTS(iSubject)};
        edata = cellstr(Pt.LFP.DIM(2).label);
        zc = ismember(ecoord, edata);
        ecoord = ecoord(zc);
        
        ntp = (WindowSizeInMilliseconds / BoxCarSize);
        xyz = XYZ{SUBJECTS(iSubject)}(zc,:);
        xyzc = mat2cell(xyz, ones(size(xyz,1),1), size(xyz,2));
        xyzcm = repmat(xyzc(:)', ntp, 1);
        xyzcv = xyzcm(:);
        xyz_repeated = cell2mat(xyzcv);
        
        trode = ELECTRODE{SUBJECTS(iSubject)}(zc);
        trodecm = repmat(trode(:)', ntp, 1);
        trodecv = trodecm(:);
        trode_repeated = cell2mat(trodecv);
        
        COORDS = struct('orientation','mni','labels',{trode_repeated},'ijk',[],'ind',[],'xyz',xyz_repeated);

        tagfmt = F.sessiontag;
        stim_order.OnsetIndex = zeros(size(stim_order,1),1);
        for iSession = 1:nsessions
            tagname = sprintf(tagfmt,iSession);
            z = stim_order.Session == iSession;
            stim_order.OnsetIndex(z) = Pt.(tagname);
        end

        Hz = 1 / Pt.LFP.DIM(1).interval; % ticks per second
        boxcar_size  = (BoxCarSize / 1000) * Hz;
        window_start = (WindowStartInMilliseconds / 1000) * Hz;
        window_size  = ( WindowSizeInMilliseconds / 1000) * Hz; % in ticks (where a tick is a single time-step).
        baseline_size  = ( BaselineSizeInMilliseconds / 1000) * Hz; % in ticks (where a tick is a single time-step).

        % Will return a session -by- electrode cell array, each containing a
        % trial -by- time matrix.
        ECA = pull_trial_profiles(Pt.LFP, stim_order, [window_start, window_size], baseline_size, ecoord);
%         ECA = arrangeElectrodeData(Pt.LFP.DATA, onsetIndex, [window_start, window_size]);
        % Average over sessions?
        if AverageOverSessions
            % data array is sessions x items x time, and mean operates over
            % the first dimension by default.
            for iElectrode = 1:numel(ECA)
                ECA(iElectrode).data = squeeze(mean(ECA(iElectrode).data));
            end
        else
            % permute will make items the first dimension, and the reshape
            % will then block rows by session
            for iElectrode = 1:numel(ECA)
                ECA(iElectrode).data = reshape((permute(ECA(iElectrode).data,[2,1,3])), [], window_size);
            end
        end

        % Boxcar Average time-points
        % ==========================
        % Algorithm :
        % 1. Trim the window so that it is evenly divisible by the boxcar.
        % 2. Transpose the data to make it time by item.
        % 3. Reshape to be boxcar_size by adj_window_size / boxcar_size by
        % item.
        % 4. Take mean over first dimension (default of mean()).
        % 5. Squeeze out the now singleton first dimension
        % 6. Transpose back to item by time.
        a = 1;
        b = window_size - rem(window_size, boxcar_size);
        c = b / boxcar_size;
        for iElectrode = 1:numel(ECA)
            x = ECA(iElectrode).data(:,a:b);
            r = size(x, 1);
            ECA(iElectrode).data = squeeze(mean(reshape(x',boxcar_size,c,r)))';
        end
        
        % Compose a single matrix, where rows are items and columns are
        % time points nested within electrodes. This is done by simply
        % concatenating the matrices associated with each electrode
        % columnwise.
        X = cat(2, ECA.data);
        
        % This is the old approach... will be deleted soon.
%         nElectrodes = size(ECA,2);
%         for iElectrode = 1:nElectrodes
%             for iSession = 1:4
%                 ECA{iSession,iElectrode} = ECA{iSession,iElectrode}(stim_sort_ix{iSession},:);
%                 if boxcar_size > 1
%                     ECA{iSession,iElectrode} = boxcarmean(ECA{iSession,iElectrode},boxcar_size,'KeepPartial',0);
%                 end
%             end
%             % Average Sessions
%             if AverageOverSessions
%                 tmp = cat(3,ECA{:,iElectrode});
%                 ECA{1,iElectrode} = mean(tmp,3);
%             end
%         end
%         if AverageOverSessions
%             ECA(2:end,:) = [];
%         end
%         X = cell2mat(ECA);

        [~,reduxFilter] = removeOutliers(X);
        y = COORDS.xyz(:,2);
        m = median(y);
        
        M = selectbyfield(metadata, 'subject', SUBJECTS(iSubject));
        M = registerFilter(M, 'rowfilter', 1, reduxFilter.words);
        M = registerFilter(M, 'colfilter', 2, reduxFilter.voxels);
        M = registerFilter(M, 'anterior',  2, y > m);
        M = registerFilter(M, 'posterior', 2, y <= m);
        tmp = load(fullfile(FILTER_DIR,'DilkineaSplit','filter.mat'),'filter');
        if AverageOverSessions
            M = registerFilter(M, 'animate', 1, [true(50,1);false(50,1)]);
            M = registerFilter(M, 'inanimate', 1, [false(50,1);true(50,1)]);
            M = registerFilter(M, 'DilkineaSplit1', 1, tmp.filter == 0);
            M = registerFilter(M, 'DilkineaSplit2', 1, tmp.filter == 1);
        else
            M = registerFilter(M, 'animate', 1, repmat([true(50,1);false(50,1)],nsessions,1));
            M = registerFilter(M, 'inanimate', 1, repmat([false(50,1);true(50,1)],nsessions,1));
            M = registerFilter(M, 'DilkineaSplit1', 1, repmat(tmp.filter,nsessions,1) == 0);
            M = registerFilter(M, 'DilkineaSplit2', 1, repmat(tmp.filter,nsessions,1) == 1);
        end

        M.coords = COORDS;
        M.ncol = size(X,2);
        M.samplingrate = Hz;
        metadata = replacebyfield(metadata, M, 'subject', SUBJECTS(iSubject));

        if exist(dpath_out,'file') && ~OVERWRITE;
            fprintf('Subject %d not written to disk, %s because output already exists.\n',SUBJECTS(iSubject),datacode)
        else
            save(dpath_out, 'X');
            fprintf('Subject written to %s\n', dpath_out);
            if WriteIndividualMetadata
                OneSubject.metadata = M;
                metapathout = fullfile(DATA_DIR_OUT,sprintf('metadata_%s_%02d.mat',datacode,SUBJECTS(iSubject)));
                save(metapathout, '-struct', OneSubject);
                fprintf('Metadata written to %s\n', metapathout);
            end
        end
    end
    %% Save metadata
    if ~WriteIndividualMetadata
        metapathout = fullfile(DATA_DIR_OUT,sprintf('metadata_%s.mat',datacode));
        fprintf('Metadata written to %s\n', metapathout);
        save(metapathout,'metadata');
    end
end
