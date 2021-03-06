function [files] = list_files(path)
  if nargin == 0; path = '.'; end
  if ~exist(path,'dir')
    error('ECoG_DataPrep:IO','%s: cannot access %s: No such file or directory', mfilename, path);
  end
  FileInfo = dir(path);
  FileInfo = FileInfo([FileInfo.isdir]==0);
  FileInfo = FileInfo(~cellfun('isempty', {FileInfo.date}));
  FileInfo = FileInfo(~cellfun(@(x) strncmp('.', x, 1), {FileInfo.name}));
  files = {FileInfo.name};
end
