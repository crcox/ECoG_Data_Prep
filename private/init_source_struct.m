function [ s ] = init_source_struct(nTicks,electrodeLabels,interval)
  nElectrodes = size(electrodeLabels,1);
  s = struct(...
    'SUBJECT','','NAVE', [],'DATA', zeros(nTicks,nElectrodes),...
    'DIM', [...
      struct(...
        'type','interval',...
        'name','time',...
        'interval',interval,...
        'scale',zeros(1,nTicks),...
        'label',[]...
      ),...
      struct(...
        'type','nominal',...
        'name','name',...
        'interval',[],...
        'scale',[],...
        'label',electrodeLabels...
      ),...
    ],...
    'DATE',0,'HISTORY',[],'MISC',[]...
  );
end
