%% Setup
DATA_DIR = '/media/chris/Data/ECOG';
META_DIR = '/media/chris/Data/ECOG/everything_else';
SUBJECTS = 1:10; % subjects 11 and 12 do not have labeled coordinates.

% IMPORTANT: The following 4 variables will effect how the data are
% processed.
AverageOverSessions = true;
BoxCarSize = 20; % If this is greater than 1, this is the number of 
                 % subsequent time points that will be averaged together.
WindowStartInMilliseconds = 0;  % zero means "window starts at stimulus onset",
                  % one means "window starts one tick post stimulus onset".
WindowSizeInMilliseconds = 1000;

%% Read presentation order and stim labels from file
% All subjects have the same order
file_stim_order = fullfile(META_DIR, 'picture_namingERP_order.csv');
file_stim_key = fullfile(META_DIR, 'picture_namingERP_key.csv');

fid = fopen(file_stim_order);
stim_order_head = textscan(fid,'%s %s %s',1,'Delimiter',',');
stim_order = textscan(fid,'%d %d %d','Delimiter',',');
fclose(fid);

fid = fopen(file_stim_key);
stim_key_head = textscan(fid,'%s %s',1,'Delimiter',',');
stim_key = textscan(fid,'%d %s','Delimiter',',');
fclose(fid);

nitems = numel(stim_key{2});

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
embedding = csvread('Similarity/NEXT_CK_KIND_5D.csv');
plot_similarity_decompositions(embedding);
S = corr(embedding','type','Pearson');

%% Check dimensionality of decomposition
addpath('~/src/WholeBrain_RSA/src');
tau = 0.2;
[~,r] = sqrt_truncate_r(S, tau);
fprintf('-----\n');
fprintf('To approximate S with error tolerance %.2f, %d dimensions are required.\n', tau, r);
fprintf('-----\n');

%% Load (basal) coordinates
[subject,electrode,x,y,z] = importcoordinates('MNI_basal_electrodes_Pt01_10_w_label.csv');
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

%% Define metadata
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

NSUBJ = numel(SUBJECTS);
NCOND = 2;
animate = [zeros(50,1);ones(50,1)];
for iSubj = 1:NSUBJ
    % TARGETS
    TARGETS = struct(...
        'label', {'animate','semantic'},...
        'type', {'category','similarity'},...
        'sim_source',{[],'NEXT'},...
        'sim_metric',{[],'correlation'},...
        'target',{animate,S}...
    );

    % FILTERS
    % This is kind of a place holder. When we remove outliers, we'll want
    % to create filters for that.
    FILTERS = struct('label',[],'dimension',[],'filter',[]);

    % CV
    nschemes = 10;
    nfolds = 10;
    SCHEMES = zeros(nitems, nschemes);
    for iScheme = 1:nschemes
        c = cvpartition(animate,'KFold', nfolds);
        for iFold = 1:nfolds
            SCHEMES(:,iScheme) = SCHEMES(:,iScheme) + (test(c, iFold) * iFold);
        end
    end

    COORDS = struct('orientation','mni','labels',ELECTRODE{iSubj},'ijk',[],'ind',[],'xyz',XYZ{iSubj});

    % ---------
    metadata(iSubj).AverageOverSessions = AverageOverSessions;
    metadata(iSubj).BoxCarSize = BoxCarSize;
    metadata(iSubj).WindowStartInMilliseconds = WindowStartInMilliseconds;
    metadata(iSubj).WindowSizeInMilliseconds = WindowSizeInMilliseconds;
    metadata(iSubj).targets = TARGETS;
    metadata(iSubj).filters = FILTERS;
    metadata(iSubj).coords = COORDS;
    metadata(iSubj).cvind = SCHEMES;
    if AverageOverSessions == 1;
        metadata(iSubj).sessions = [];
        metadata(iSubj).nrow = 100;
    else
        metadata(iSubj).sessions = stim_order{1};
        metadata(iSubj).nrow = 400;
    end
    metadata(iSubj).ncol = 0; % will be set later
end

%% Load And Process Data
% NOTE: The filenames and field names will not work for all subjects!!!
for iSubj=1:NSUBJ
    sdir = sprintf('Pt%02d',iSubj);
    sfile = sprintf('namingERP_Pt%02d_refD14.mat',iSubj);
    spath = fullfile(DATA_DIR,sdir,sfile);
    fprintf('Loading %s...\n', spath);
    Pt = load(spath);
    Pt.LFP = Pt.namingERP_data_PtYK_Pt01_refD14;
    Pt = rmfield(Pt,'namingERP_data_PtYK_Pt01_refD14');

    ecoord = ELECTRODE{iSubj};
    edata = cellstr(Pt.LFP.DIM(2).label);

    zd = ismember(edata, ecoord);
    zc = ismember(ecoord, edata);

    Pt.LFP.DATA = Pt.LFP.DATA(:,zd);

    onsetIndex = {Pt.tag_ss01_all,Pt.tag_ss02_all,Pt.tag_ss03_all,Pt.tag_ss04_all};

    Hz = 1 / Pt.LFP.DIM(1).interval; % ticks per second
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
            if BoxCarSize > 1
                M{iSession,iElectrode} = boxcarmean(M{iSession,iElectrode},BoxCarSize,'KeepPartial',0);
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
    metadata(iSubj).filters(end+1) = struct('label','rowfilter','dimension',1,'filter',reduxFilter.words);
    metadata(iSubj).filters(end+1) = struct('label','colfilter','dimension',2,'filter',reduxFilter.voxels);
    metadata(iSubj).ncol = size(X,2);
    metadata(iSubj).samplingrate = Hz;

    dpath_out = fullfile(DATA_DIR_OUT, sprintf('s%02d.mat',iSubj));
    save(dpath_out, 'X');
end

%% Save metadata
save(fullfile(DATA_DIR_OUT,'metadata.mat'));
