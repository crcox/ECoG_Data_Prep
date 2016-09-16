function [ metadata ] = registerFilter(metadata, label, dimension, filter)
  key = {'label','dimension'};
  val = {label, dimension};
  nkey = numel(key);
  tmp = reshape([key;val], nkey*2, 1);
  for i = 1:numel(metadata)
    if numel(fieldnames(metadata(i).filters)) > 0;
      nfilters = numel(metadata(i).filters);
      z = false(nfilters, nkey);
      for j = 1:nkey
        k = key{j};
        v = val{j};
        z(:,j) = strcmp({metadata(i).filters.(k)}, v);
      end
      if any(all(z))
        metadata(i).filters(all(z)) = struct(tmp{:},'filter', filter);
      else
        metadata(i).filters(end+1) = struct(tmp{:},'filter', filter);
      end
    else
      metadata(i).filters = struct(tmp{:},'filter', filter);
    end
  end
end


