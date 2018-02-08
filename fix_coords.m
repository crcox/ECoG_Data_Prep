% fix_coordinates
% Coordinates are not replicated to account for the number of timepoints in
% a data file.
%
% Notes from ECoG_Data_Setup/setup_data
% % Will return a session -by- electrode cell array, each containing a
% % trial -by- time matrix.
% ECA = arrangeElectrodeData(Pt.LFP.DATA, onsetIndex, [window_start, window_size]);
%  ...
% X = cell2mat(ECA);
%
% This shows that columns are arranged by *time* first, then by electrode.

boxcar = 10; %ms per time step
windowsize = 50:10:1000;
ntp = windowsize / boxcar;
metadata_repeated = cell(numel(windowsize),8);
for i = 1:numel(windowsize);
    ws = windowsize(i);
    d = sprintf('%04d', ws);
    dpath = fullfile(d, 'metadata_raw.mat');
    copyfile(dpath, strcat(dpath,'.orig'));
    disp(dpath);
    load(dpath);
    disp(cellfun(@(x) size(x.xyz, 1), {metadata.coords}) * ntp(i));
    for j = 1:numel(metadata)
        xyz = metadata(j).coords.xyz;
        xyzc = mat2cell(xyz, ones(size(xyz,1),1), size(xyz,2));
        xyzcm = repmat(xyzc(:)', ntp(i), 1);
        xyzcv = xyzcm(:);
        xyz_repeated = cell2mat(xyzcv);
        metadata(j).coords.xyz = xyz_repeated;
    end
    save(dpath, 'metadata');
end

    
    
end