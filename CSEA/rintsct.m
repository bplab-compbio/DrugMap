function [nl] = rintsct(nd,sd,vp,vi)
% "Random intersect" draws random elements from a vector
% ... and intersects draws with a pre-specified vector of interest
%   nd:= number of draws
%   sd:= size of draws
%   p:= vector or "pile" of elements to be drawn from
%   i:= vector or "interesting" elements to be intersected

    % Pre-allocate memory
    nl = nan(nd,1);
    
    % Intersect random draws of elements with a vector of interest
    for j = 1:nd
        nl(j) = length(intersect(p(randi([1,length(p)],[1,sd])),i)); 
    end

end
