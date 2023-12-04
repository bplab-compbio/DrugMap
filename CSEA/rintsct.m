function [nl] = rintsct(nd,sd,p,i)
% "Random intersect" draws random elements from a vector
% ... and tracks size of overlap w/ elements of interest
%   nd:= number of draws
%   sd:= size of draws
%   p:= vector or "pile" of elements to be drawn from
%   i:= vector of "interest"

    % Pre-allocate memory
    nl = nan(nd,1);
    
    % Intersect random draws of elements with vector of interest
    for j = 1:nd
        nl(j) = length(intersect(p(randi([1,length(p)],[1,sd])),i)); 
    end

end
