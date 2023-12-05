function [result] = csea(n,s,np,names,sets)
% Perform Cysteine Set Enrichment Analysis (CSEA)
%   n:= string vector of cysteines which will form the null distribution
%   s:= string vector of cysteines of interest
%   np:= number of permutations requested by user
%   names:= string vector which contains set names
%   sets:= cell vector, wherein each cell contains a string vector of cysteines

    if nargin ~= 5
        % check for full complement of inputs
        error("incorrect number of inputs, please check!")
    else
        % initialize empty matrix
        d = nan(length(sets),5);

        % for each set under study
        for x = 1:length(names)            

            % grab cysteines in set 
            c = sets{x};

            % get intersection btw. cys of interest and cysteine set
            t = length(intersect(s,c));
        
            % if intersection is not empty
            if t > 0                 

                % generate null distribution
                nl = rintsct(np,length(s),n,c);

                % integrate
                p = perm_integrate(nl,t);

                % check if result is potentially significant; if so, dial in p-value with 1000 perms
                % this step is optional and can be omitted to improve speed
                if -log10(p) > 2
                    nl = rintsct(1000,length(s),n,c);
                    p = perm_integrate(rintsct(1000,length(s),n,c),t); 
                end

                % collect results
                d(x,1:4) = [p,t,(t+1)/(median(nl,"omitnan")+1),length(c)];

            end

            % track progress
            disp(100 * x/length(sets) + "% done!")
        end    

        % delete results with 0 overlap or skipped p-value calculation
        names(any(isnan(d(:,1:2)),2)) = [];
        d(any(isnan(d(:,1:2)),2),:) = [];  

        % calculate Benjamini-Hochberg FDR 
        % (Benjamini, Y., and Hochberg, Y. 1995. Controlling the false discovery rate: A practical and powerful approach to multiple testing. J. Royal Stat. Soc. 57:289–300.)        
        d(:,5) = mafdr(d(:,1),"BHFDR","true");
      
        % summarize results in table
        result = [table(names),splitvars(table(d))]; 
        result.Properties.VariableNames = ["set name"
            "p-value"
            "overlap (#)"
            "ES (enrichment score)"
            "size of cysteine set"
            "FDR"];
    end
end
