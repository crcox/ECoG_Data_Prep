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
nsessions = numel(unique(stim_order{1}));;

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
    metadata(iSubj).filters = FILTERS;
    metadata(iSubj).coords = COORDS;
    if AverageOverSessions == 1;
        for iTarget = 1:numel(TARGETS)
            switch lower(TARGETS(iTarget).type)
            case 'category'
                TARGETS(iTarget).target = repmat(TARGETS(iTarget).target,nsessions,1);
            case 'similarity'
                TARGETS(iTarget).target = repmat(TARGETS(iTarget).target,nsessions,nsession);
            end
            SCHEMES = repmat(SCHEMES,nsessions,1);
        end
        metadata(iSubj).sessions = [];
        metadata(iSubj).nrow = 100;
    else
        metadata(iSubj).sessions = stim_order{1};
        metadata(iSubj).nrow = 400;
    end
    metadata(iSubj).targets = TARGETS;
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
  'namingERPdata_Pt07.mat',{'namingERPdataPt07','namingERPdataPt07_ss0304'},{[1,2],[3,4]},'Tag_ss%02d';...
  'namingERPdata_Pt08.mat',{'namingERPdataPt08'},{[1,2,3,4]},'tag_%02d';...
  'namingERPdata_Pt09.mat',{'namingERPdataPt09'},{[1,2,3,4]},'tag%02d';...
  'namingERPdata_Pt10.mat',{'namingERPdataPt10'},{[1,2,3,4]},'tagall%02d';...
};
filelist_ref = {...
  'namingERP_Pt01_refD14.mat',{'namingERP_data_PtYK_Pt01_refD14'},{[1,2,3,4]},'tag_ss%02d_all';...
  [],[],[],[];...
  [],[],[],[];...
  'namingERP_PtMA_REF4.mat',{'namingERP_data_ss01ss02_REF4','namingERP_data_ss03ss04_REF4'},{[1,2],[3,4]},'tag_ss%02d';...
  'namingERP_Pt05_ref02.mat',{'namingERP_data_Pt05_ref02'},{[1,2,3,4]},'tag_ss%02d_all';...
  'namingERP_Pt06_ref01.mat',{'namingERP_data_Pt06_ref01'},{[1,2,3,4]},[];...
  'namingERPdata_Pt07_ref02.mat',{'namingERPdataPt07_ss0102_ref02','namingERPdataPt07_ss0304_ref02'},{[1,2],[3,4]},'Tag_ss%02d';...
  'namingERPdata_Pt08_ref03.mat',{'namingERPdataPt08_ref03'},{[1,2,3,4]},'tag_%02d';...
  [],[],[],[];...
  'namingERPdata_Pt10_ref02.mat',{'namingERPdataPt10_ref02'},{[1,2,3,4]},'tagall%02d';...
};

for iRef = 1:2
  if iRef == 1
    filelist = filelist_raw;
    fmt = 's%02d_raw.mat';
  else
    filelist = filelist_ref;
    fmt = 's%02d_ref.mat';
  end
  for iSubj=1:NSUBJ
      sdir = sprintf('Pt%02d',iSubj);
      sfile = filelist{iSubj,1};
      spath = fullfile(DATA_DIR,sdir,sfile);
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
        Pt.LFP = init_source_struct(nTicks,electrodeLabels,interval);
        for iChunk = 1:nChunks
          cvar = filelist{iSubj,2}{iChunk};
          sessions = filelist{iSubj,3}{iChunk};
          nSessions = numel(sessions);
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

            for iSession = 1:nSessions
              tag = sprintf(tagfmt,iSession);
              Pt.(tag) = Pt.(tag) + psize;
            end

            psize = b;
            pscale = max(Pt.(cvar).DIM(1).scale);
            Pt = rmfield(Pt,cvar);
          end
        end
      else
        cvar = varnames{1};
        Pt.LFP = Pt.(cvar);
        Pt = rmfield(Pt,cvar);
      end

      ecoord = ELECTRODE{iSubj};
      edata = cellstr(Pt.LFP.DIM(2).label);

      zd = ismember(edata, ecoord);
      zc = ismember(ecoord, edata);

      Pt.LFP.DATA = Pt.LFP.DATA(:,zd);

      onsetIndex = cell(1,4);
      tagfmt = filelist{iSubj,3};
      for iSession = 1:4
        tagname = sprintf(tagfmt,iSession);
        onsetIndex{iSession} = Pt.(tagname);
      end

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

      dpath_out = fullfile(DATA_DIR_OUT, sprintf(fmt,iSubj));
      save(dpath_out, 'X');
  end
end

%% Save metadata
save(fullfile(DATA_DIR_OUT,'metadata.mat'));
