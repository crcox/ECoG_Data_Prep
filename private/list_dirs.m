function [dirs] = list_dirs(path)
  if nargin == 0; path = '.'; end
  if ~exist(path,'dir')
    error('ECoG_DataPrep:IO','%s: cannot access %s: No such file or directory', mfilename, path);
  end
  FileInfo = dir(path);
  FileInfo = FileInfo([FileInfo.isdir]==1);
  FileInfo = FileInfo(~cellfun('isempty', {FileInfo.date}));
  % FileInfo = FileInfo(~cellfun(@(x) strcmp('.', x), {FileInfo.name}));
  % FileInfo = FileInfo(~cellfun(@(x) strcmp('..', x), {FileInfo.name}));
  FileInfo = FileInfo(~cellfun(@(x) strncmp('.', x, 1), {FileInfo.name}));
  dirs = {FileInfo.name};
end