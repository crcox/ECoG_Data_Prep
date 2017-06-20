function update_metadata()
  datadir = '/home/chris/MRI/ECOG/avg/BoxCar/001';
  onset = dirlist(fullfile(datadir,'WindowStart'));
  for iOnset = 1:numel(onset);
    t0 = onset{iOnset};
    n = fprintf('t0: %d\n', str2double(t0));
    fprintf('%s\n',repmat('-',1,n));

    metadatafile_raw = fullfile(datadir,'WindowStart',t0,'WindowSize','0050','metadata_raw.mat');
    metadatafile_ref = fullfile(datadir,'WindowStart',t0,'WindowSize','0050','metadata_ref.mat');
    metadatafile_raw_bck = fullfile(datadir,'WindowStart',t0,'WindowSize','0050','metadata_raw.mat.bck');
    metadatafile_ref_bck = fullfile(datadir,'WindowStart',t0,'WindowSize','0050','metadata_ref.mat.bck');

    if exist(metadatafile_raw,'file')
      copyfile(metadatafile_raw,metadatafile_raw_bck);
      load(metadatafile_raw, 'metadata');
      for iSubj = 1:numel(metadata)
        y = metadata(iSubj).coords(1).xyz(:,2);
        m = median(y);
        metadata(iSubj).filters(end+1) = struct('label','anterior','dimension',2,'filter',y > m);
        metadata(iSubj).filters(end+1) = struct('label','posterior','dimension',2,'filter',y <= m);
        metadata(iSubj).filters = squeeze(metadata(iSubj).filters);
      end
      save(metadatafile_raw, 'metadata');
    end

    if exist(metadatafile_ref,'file')
      copyfile(metadatafile_ref,metadatafile_ref_bck);
      load(metadatafile_ref, 'metadata');
      for iSubj = 1:numel(metadata)
        y = metadata(iSubj).coords(1).xyz(:,2);
        m = median(y);
        metadata(iSubj).filters(end+1) = struct('label','anterior','dimension',2,'filter',y > m);
        metadata(iSubj).filters(end+1) = struct('label','posterior','dimension',2,'filter',y <= m);
        metadata(iSubj).filters = squeeze(metadata(iSubj).filters);
      end
      save(metadatafile_ref, 'metadata');
    end
  end
end

function d = dirlist(path)
  x = dir(path);
  % Select directories ...
  z = [x.isdir];
  % ... that do not have a '.' as the first character in their name.
  z = z & ~strcmp(cellfun(@(x) x(1), {x.name},'Unif',0),'.');
  d = {x(z).name};
end
