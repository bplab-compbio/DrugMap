function [result] = csea(n,s,np,names,sets,sidedness)
% Perform Cysteine Set Enrichment Analysis (CSEA)
%   n:= string vector of cysteines which will form the null distribution
%   s:= string vector of cysteines of interest
%   np:= number of permutations requested by user
%   names:= string vector which contains set names
%   sets:= cell vector, wherein each cell contains a string vector of cysteines
%   sidedness:= 1 or 2; if 2, performs two-tailed integration for p-value; elseif sidedness == 1, ignore density in left tail


    if nargin ~= 6
        % check for full complement of inputs
        error("incorrect number of inputs, please check!")
    else
        % initialize empty matrix
        d = nan(length(sets),4);

        % for each set under study
        for x = 1:length(names)            

            % grab cysteines in set 
            c = sets{x};

            % get intersection btw. cys of interest and cysteine set
            t = length(intersect(s,c));
        
            % get null distribution
            nl = nan(np,1);
            for j = 1:np, nl(j) = length(intersect(n(randi([1,length(n)],[1,length(s)])),c)); end
    
            % integrate null distribution and report results
            d(x,1:3) = [perm_integrate(nl,t,sidedness),t,length(c)];

            % keep track of progress
            disp(x/length(sets) + " % done!")
        end    
    
        % get Benjamini-Hochberg FDR 
        % (Benjamini, Y., and Hochberg, Y. 1995. Controlling the false discovery rate: A practical and powerful approach to multiple testing. J. Royal Stat. Soc. 57:289–300.)        
        d(:,4) = mafdr(d(:,1),"BHFDR","true");
    
        % summarize results in table
        result = [table(names),splitvars(table(d(:,[1,4,2,3])))]; 
        result.Properties.VariableNames = ["set name","p-value","FDR","overlap (#)","size of cysteine set"];
    end
end
