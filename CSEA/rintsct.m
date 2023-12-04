function [nl] = rintsct(nd,sd,vp,vi)
% Draw random elements from a vector
% Intersect random draws with vector of interest
% Tabulate size of intersection along the way
%   nd:= number of draws
%   sd:= size of draws
%   vp:= vector or "pile" of elements to be drawn from (i.e. vector_pile)
%   vi:= vector to be intersected with (i.e. vector_interest)

    % Pre-allocate memory
    nl = nan(nd,1);
    
    % Intersect random draws of elements with a vector of interest
    for j = 1:nd
        nl(j) = length(intersect(vp(randi([1,length(vp)],[1,sd])),vi)); 
    end

end
