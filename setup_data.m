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
        addParameter(p, 'subjects', '1 2 3 5 7 8 9 10', @ischar);
        addParameter(p, 'boxcar', '0', @ischar);
        addParameter(p, 'average', '1', @ischar);
        addParameter(p, 'dataroot', '/mnt/sw01-home01/mbmhscc4/scratch/data/Naming_ECoG/avg');
        addParameter(p, 'metaroot', '/mnt/sw01-home01/mbmhscc4/scratch/data/Naming_ECoG/meta');
        addParameter(p, 'cvpath', []);
        addParameter(p, 'overwrite', '0');
        parse(p, varargin{:});
        if isempty(p.Results.onset)
            error('A window onset time must be provided.');
        end
        if isempty(p.Results.duration)
            error('A window duration must be provided.');
        end
        fprintf('      Window onset (ms): %s\n', p.Results.onset);
        fprintf('   Window duration (ms): %s\n', p.Results.duration);
        fprintf('               subjects: %s\n', p.Results.subjects);
        fprintf('            Boxcar Size: %s\n', p.Results.boxcar);
        fprintf('  Average over sessions: %s\n', p.Results.average);
        fprintf(' Data root (for output): %s\n', p.Results.dataroot);
        fprintf('  Meta root (for input): %s\n', p.Results.metaroot);
        fprintf('Overwrite existing data: %s\n', p.Results.overwrite);
        
        %% Setup
        WindowStartInMilliseconds = str2double(p.Results.onset);
        WindowSizeInMilliseconds = str2double(p.Results.duration);
        BoxCarSize = str2double(p.Results.boxcar);
        AverageOverSessions = str2double(p.Results.average);
        OVERWRITE = str2double(p.Results.overwrite);
        clean_string = regexprep(p.Results.subjects, ',? *', ' ');
        SUBJECTS = str2double(strsplit(clean_string));
    else
        % In this case, arguments can be handled in their expected types
        % (numbers as numbers, etc). The input may be coming directly from
        % the Matlab command line or script, or parsed from a json file.
        addParameter(p, 'onset',[]);
        addParameter(p, 'duration',[]);
        addParameter(p, 'subjects', [1:3,5,7:10], @isnumeric);
        addParameter(p, 'boxcar', 0);
        addParameter(p, 'average', 1);
        addParameter(p, 'dataroot', 'D:\ECoG\KyotoNaming\data');
        addParameter(p, 'metaroot', 'C:\Users\mbmhscc4\MATLAB\ECOG\naming\data');
        addParameter(p, 'cvpath', []);
        addParameter(p, 'overwrite', 0);
        if isjsoncfg
            % Input read from a json file will be represented as a
            % structure. The inputParser cannot handle this structure, and
            % so we need to get the structure into an amenable format.
            varargin = [fieldnames(jdat); struct2cell(jdat)];
        end
        parse(p, varargin{:});
        fprintf('      Window onset (ms): %d\n', p.Results.onset);
        fprintf('   Window duration (ms): %d\n', p.Results.duration);
        fprintf('               subjects: %s\n', strjoin(strsplit(num2str(p.Results.subjects)), ', '));
        fprintf('            Boxcar Size: %d\n', p.Results.boxcar);
        fprintf('  Average over sessions: %d\n', p.Results.average);
        fprintf(' Data root (for output): %s\n', p.Results.dataroot);
        fprintf('  Meta root (for input): %s\n', p.Results.metaroot);
        fprintf('Overwrite existing data: %d\n', p.Results.overwrite);
        
        %% Setup
        WindowStartInMilliseconds = p.Results.onset;
        WindowSizeInMilliseconds = p.Results.duration;
        BoxCarSize = p.Results.boxcar;
        AverageOverSessions = p.Results.average;
        OVERWRITE = p.Results.overwrite;
        SUBJECTS = p.Results.subjects;
    end
    DATA_DIR = p.Results.dataroot;
    META_DIR = p.Results.metaroot;
    STIM_DIR = fullfile(META_DIR,'stimuli');
    COORD_DIR = fullfile(META_DIR,'coords');
    SIM_DIR = fullfile(META_DIR,'similarity');
    TARGET_DIR = fullfile(META_DIR,'targets');
    CV_DIR = fullfile(META_DIR,'cv');
    cvpath = p.Results.cvpath;

        
    %% Define Output directory
    if AverageOverSessions == 1
        bdir = 'avg';
    else
        bdir = 'full';
    end
    DATA_DIR_OUT = fullfile(...
        DATA_DIR,...
        bdir,...
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

    fid = fopen(file_stim_order);
    [~] = textscan(fid,'%s %s %s',1,'Delimiter',','); % stim_order_head
    stim_order = textscan(fid,'%d %d %d','Delimiter',',');
    fclose(fid);

    fid = fopen(file_stim_key);
    [~] = textscan(fid,'%s %s',1,'Delimiter',','); % stim_key_head
    stim_key = textscan(fid,'%d %s','Delimiter',',');
    fclose(fid);

    nitems = numel(stim_key{2});
    nsessions = numel(unique(stim_order{1}));

    %% ensure stimulus labels are sorted by their index
    [stim_key{1},ix] = sort(stim_key{1});
    stim_key{2} = stim_key{2}(ix);

    %% Generate sort indexes
    x = tabulate(stim_order{1});
    stim_order_ix = mat2cell(stim_order{3}, x(:,2), 1);
    stim_sort_ix = cell(size(x,1),1);
    for i = 1:size(x,1)
        [~,stim_sort_ix{i}] = sort(stim_order_ix{i});
    end

    %% load similarity structure
    % NEXT (judge similarity in kind based on word)
    % The rows in the embedding are already in an order that matches the order
    % in the stim_key.
    % S = csvread(fullfile(SIM_DIR,'LeuvenNorm_cosine_similarity.csv'));
%     if usejava('jvm')
%         plot_similarity_decompositions(embedding);
%     end
%     S = corr(embedding','type','Pearson');

    %% Check dimensionality of decomposition
%     addpath('~/src/WholeBrain_RSA/src');
%     tau = 0.2;
%     [~,r] = sqrt_truncate_r(S, tau);
%     fprintf('-----\n');
%     fprintf('To approximate S with error tolerance %.2f, %d dimensions are required.\n', tau, r);
%     fprintf('-----\n');

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
      'stimuli',stim_key(2),...
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
              metadata = installSimilarityStructure(metadata, structure_matrix, structure_labels, type_sim, source, metric_label);
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
    NCOND = 2;
    for iSubject = 1:NSUBJ
        iSubj = SUBJECTS(iSubject);
        % FILTERS
        % This is kind of a place holder. When we remove outliers, we'll want
        % to create filters for that.
        FILTERS = struct('label',[],'dimension',[],'filter',[]);

        % ---------
        if AverageOverSessions == 1;
            metadata(iSubj).sessions = [];
            metadata(iSubj).nrow = 100;
        else
            metadata(iSubj).sessions = stim_order{1};
            metadata(iSubj).nrow = 400;
        end
        metadata(iSubj).AverageOverSessions = AverageOverSessions;
        metadata(iSubj).BoxCarSize = BoxCarSize;
        metadata(iSubj).WindowStartInMilliseconds = WindowStartInMilliseconds;
        metadata(iSubj).WindowSizeInMilliseconds = WindowSizeInMilliseconds;
        if numel(fieldnames(metadata(iSubject).filters)) == 0
          metadata(iSubj).filters = FILTERS;
        end
%         metadata(iSubj).targets = TARGETS; % handled above
        metadata(iSubj).cvind = SCHEMES;
        metadata(iSubj).ncol = 0; % will be set later
    end

    %% Load And Process Data
    % NOTE: In the source data, naming conventions are not consistent. The
    % following structures explicitly represent all the file and field names that
    % need to be referenced.
    % COLUMN KEY:
    %  1. filename
    %  2. variable names
    %  3. sessions per variable
    %  4. format of tag variables
    filelist_hdr = {'filename','varnames','sessions','tagfmt'};
    filelist_raw = {...
      'namingERP_Pt01.mat',{'namingERP_data_PtYK_Pt01'},{[1,2,3,4]},'tag_ss%02d_all';...
      'namingERP_Pt02.mat',{'namingERP_data_Pt02'},{[1,2,3,4]},'Tag_ss%02d_all';...
      'namingERP_Pt03.mat',{'namingERP_data_Pt03'},{[1,2,3,4]},'Tag_ss%02d_all';...
      [],[],[],[];...
      'namingERP_Pt05.mat',{'namingERP_data_Pt05'},{[1,2,3,4]},'tag_ss%02d';...
      'namingERP_Pt06.mat',{'namingERP_data_Pt06'},{[1,2,3,4]},[];...
      'namingERP_Pt07.mat',{'namingERPdataPt07','namingERPdataPt07_ss0304'},{[1,2],[3,4]},'Tag_ss%02d';...
      'namingERP_Pt08.mat',{'namingERPdataPt08'},{[1,2,3,4]},'tag_%02d';...
      'namingERP_Pt09.mat',{'namingERPdataPt09'},{[1,2,3,4]},'tag%02d';...
      'namingERP_Pt10.mat',{'namingERPdataPt10'},{[1,2,3,4]},'tagall%02d';...
    };
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

    for iRef = 1%:2
        if iRef == 1
            filelist = filelist_raw;
            fmt = 's%02d_raw.mat';
            mode = 'raw';
        else
            filelist = filelist_ref;
            fmt = 's%02d_ref.mat';
            mode = 'ref';
        end
        for iSubject=1:NSUBJ
            iSubj = SUBJECTS(iSubject);
            sdir = sprintf('Pt%02d',iSubj);
            sfile = filelist{iSubj,1};
            if isempty(sfile)||isempty(filelist{iSubj,4});
                fprintf('Skipping subject %d, %s because of missing data.\n',iSubj,mode);
                continue;
            else
                fprintf('Beginning subject %d, %s.\n',iSubj,mode);
            end
            dpath_out = fullfile(DATA_DIR_OUT, sprintf(fmt,iSubj));
            if exist(dpath_out,'file') && ~OVERWRITE;
                fprintf('Skipping subject %d, %s because output already exists.\n',iSubj,mode)
                continue
            end
            spath = fullfile(DATA_DIR,'raw',sdir,sfile);
            fprintf('Loading %s...\n', spath);
            Pt = load(spath);

            nChunks = numel(filelist{iSubj,2});
            if nChunks > 1
                nTicks = 0;
                for iChunk = 1:nChunks
                    cvar = filelist{iSubj,2}{iChunk};
                    nTicks = nTicks + size(Pt.(cvar).DATA,1);
                end
                interval = Pt.(cvar).DIM(1).interval;
                electrodeLabels = Pt.(cvar).DIM(2).label;
                Pt.LFP(1) = init_source_struct(nTicks,electrodeLabels,interval);
                for iChunk = 1:nChunks
                    cvar = filelist{iSubj,2}{iChunk};
                    sessions = filelist{iSubj,3}{iChunk};
                    tagfmt = filelist{iSubj,4};
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
                cvar = filelist{iSubj,2}{1};
                Pt.LFP = Pt.(cvar);
                Pt = rmfield(Pt,cvar);
            end

            ecoord = ELECTRODE{iSubj};
            edata = cellstr(Pt.LFP.DIM(2).label);

            [zc,zd] = ismember(ecoord, edata);
            zd = zd(zc);

            COORDS = struct('orientation','mni','labels',{ELECTRODE{iSubj}(zc)},'ijk',[],'ind',[],'xyz',XYZ{iSubj}(zc,:));

            Pt.LFP.DATA = Pt.LFP.DATA(:,zd);

            onsetIndex = cell(1,4);
            tagfmt = filelist{iSubj,4};
            for iSession = 1:nsessions
                tagname = sprintf(tagfmt,iSession);
                onsetIndex{iSession} = Pt.(tagname);
            end

            Hz = 1 / Pt.LFP.DIM(1).interval; % ticks per second
            boxcar_size = (BoxCarSize / 1000) * Hz;
            window_start = (WindowStartInMilliseconds / 1000) * Hz;
            window_size = (WindowSizeInMilliseconds / 1000) * Hz; % in ticks (where a tick is a single time-step).

            % Will return a session -by- electrode cell array, each containing a
            % trial -by- time matrix.
            M = arrangeElectrodeData(Pt.LFP.DATA, onsetIndex, [window_start, window_size]);

            % Sort and average time-points
            nElectrodes = size(M,2);
            for iElectrode = 1:nElectrodes
                for iSession = 1:4
                    M{iSession,iElectrode} = M{iSession,iElectrode}(stim_sort_ix{iSession},:);
                    if boxcar_size > 1
                        M{iSession,iElectrode} = boxcarmean(M{iSession,iElectrode},boxcar_size,'KeepPartial',0);
                    end
                end
                % Average Sessions
                if AverageOverSessions
                    tmp = cat(3,M{:,iElectrode});
                    M{1,iElectrode} = mean(tmp,3);
                end
            end
            if AverageOverSessions
                M(2:end,:) = [];
            end
            X = cell2mat(M);
            [~,reduxFilter] = removeOutliers(X);
            y = COORDS.xyz(:,2);
            m = median(y);
            metadata(iSubject) = registerFilter(metadata(iSubject), 'rowfilter', 1, reduxFilter.words);
            metadata(iSubject) = registerFilter(metadata(iSubject), 'colfilter', 2, reduxFilter.voxels);
            metadata(iSubject) = registerFilter(metadata(iSubject), 'anterior',  2, y > m);
            metadata(iSubject) = registerFilter(metadata(iSubject), 'posterior', 2, y <= m);
            metadata(iSubject).coords = COORDS;
            metadata(iSubject).ncol = size(X,2);
            metadata(iSubject).samplingrate = Hz;

            save(dpath_out, 'X');
        end
      %% Save metadata
      save(fullfile(DATA_DIR_OUT,sprintf('metadata_%s.mat',mode)),'metadata');
    end
end
