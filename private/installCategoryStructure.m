function metadata = installCategoryStructure(metadata, catmat, stimlabels, label, source)
  % Check that simmat is a similarity matrix or (probably) an embedding.
  [m,n] = size(catmat);
  if all([m,n]>1)
      error('catmat is multidimensional! Expecting a vector of category labels.');
  end
  type = 'category';

  % Confirm that stimulus order is the same for all subjects
  z = false(1, numel(metadata));
  for i = 2:numel(metadata)
    if ~isequal(metadata(1).stimuli, metadata(i).stimuli)
      z(i) = 1;
    end
  end
  if ~any(z(2:end))
    fprintf('PASS: All stimuli orders are the same.\n')
  else
    fprintf('FAIL: Stimuli orders for the following subjects differ from subject 1:\n')
    disp(find(z));
    fprintf('Exiting...\n')
    return
  end
        
  % The raw labels from the Kyoto group have some unusual spelling, typos,
  % and low-frequency word choices that standardized for the purpose of
  % this function
  labels_presentation_order = metadata(1).stimuli;
  if strcmp('rooster', stimlabels)
    labels_presentation_order{strcmp(labels_presentation_order, 'cockerel')} = 'rooster';
  end
  if strcmp('ladybug', stimlabels)
    labels_presentation_order{strcmp(labels_presentation_order, 'ladybird')} = 'ladybug';
  end
  if strcmp('butterfly', stimlabels)
    labels_presentation_order{strcmp(labels_presentation_order, 'batterfly')} = 'butterfly';
  end
  if strcmp('airplane', stimlabels)
    labels_presentation_order{strcmp(labels_presentation_order, 'aeroplane')} = 'airplane';
  end
  N = numel(labels_presentation_order);

  a = ismember(labels_presentation_order, stimlabels);
  [b, loc] = ismember(stimlabels, labels_presentation_order);

  % This should be empty, indicating that the labels associated with rows in
  % the similarity matrix are a subset of the full set of labels.
  if all(b)
    fprintf('PASS: Similarity matrix labels are a subset of stimuli list.\n')
  else
    fprintf('FAIL: Similarity matrix labels include some items not found in the stimuli list.\n');
    disp(stimlabels(~b));
  end

  % Need to enforce that the order of rows (and columns) in the similarity
  % matrix in presentation order. 
  [~,ix] = sort(loc);
  if isequal(labels_presentation_order(a), stimlabels(ix))
    fprintf('PASS: Lists confirmed identical after sorting.\n')
  else
    fprintf('FAIL: Lists are not identical after sorting... something broke.\n')
    %disp([sort(labels_presentation_order(a)), sort(simlabels)]);
  end
  stimlabels = stimlabels(ix);
  
  % In order to filter the data to match this subset while running the
  % analyses, we will need to:
  % 1. Store a 100x100 similarity matrix (i.e., enough slots for all items)
  %    with zeros where we do not have similarity information.
  % 2. Create a "leuven" filter, so we can subset the data (and the
  %    structure) to only include the items we have similarity data for.

  % Register the category structure
  switch type
    case 'category'
      S = zeros(N,n);
      S(a,1) = catmat(ix,1);
  end

  key = {'label','type','sim_source','sim_metric'};
  val = {label, type, source, []};
  nkey = numel(key);
  tmp = reshape([key;val], nkey*2, 1);
  for i = 1:numel(metadata)
    if numel(fieldnames(metadata(i).targets)) > 0;
      ntargets = numel(metadata(i).targets);
      z = false(ntargets, nkey);
      for j = 1:nkey
        k = key{j};
        v = val{j};
        z(:,j) = strcmp({metadata(i).targets.(k)}, v);
      end
      if any(all(z))
        metadata(i).targets(all(z)) = struct(tmp{:},'target',S);
      else
        metadata(i).targets(end+1) = struct(tmp{:},'target',S);
      end
    else
      metadata(i).targets = struct(tmp{:},'target',S);
    end
  end

  % Register a filter, if needed
  if numel(stimlabels) < numel(labels_presentation_order)
    key = {'label','dimension'};
    val = {source, 1};
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
          metadata(i).filters(all(z)) = struct(tmp{:},'filter',a);
        else
          metadata(i).filters(end+1) = struct(tmp{:},'filter',a);
        end
      else
        metadata(i).filters = struct(tmp{:},'filter',a);
      end
    end
  end
end