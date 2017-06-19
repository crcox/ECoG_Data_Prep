function [ Y ] = boxcarmean( X, boxsize, varargin)
%BOXCARMEAN Bin columns and reduce to mean of each bin.
%      X : an n x m matrix
% boxize : number of columns per bin. All bins are equal size.
% OPTIONAL
% --------
% KeepPartial : True or [false]. If the number of columns is not divisible by
%               the bin size, then partial bin is either [dropped] or averaged
%               as all other bins.
    % Defaults
    KeepPartial = 0;

    % Check varargin
    if nargin > 2
        v = reshape(varargin, 2, []);
        for i = 1:size(v,2);
            key = v{1,i};
            val = v{2,i};
            switch lower(key)
                case 'keeppartial'
                    KeepPartial = val;
                otherwise
                    error('invalid option %s.\n', key);
            end
        end
    end

    n = size(X, 2);
    nOverflow = mod(n, boxsize);
    nBoxes = (n - nOverflow) / boxsize;
    ix = repmat(1:nBoxes,boxsize,1);
    ix = ix(:);

    if KeepPartial
        nBoxes = nBoxes + 1;
        ix = [ix; repmat(nBoxes, nOverflow, 1)];
    else
        X = X(:,1:(end-nOverflow));
    end

    Y = zeros(size(X,1), nBoxes);
    for iBox = 1:nBoxes
        z = ix == iBox;
        Y(:,iBox) = mean(X(:,z),2);
    end
end
