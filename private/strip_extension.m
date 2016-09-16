function [ basepath ] = strip_extension( path )
  if iscell(path)
    basepath = cellfun(@stripext, path, 'unif', 0);
  else
    basepath = stripext(path);
  end
end
function b = stripext(p)
  [path,name,~] = fileparts(p);
  b = fullfile(path, name);
end
