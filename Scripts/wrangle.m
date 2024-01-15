% get directories
[in,out] = ld;


% pre-wrangle files so that we don't have to repeat in future
dirr = "DrugMap";
fullfile = dir(in + dirr + "\");
subdir = string({fullfile.name})'; subdir(1:2) = [];

for xx = 1:length(subdir)
    
    % only do this if we haven't made a .mat already    
%     if ~isfile(in + "DrugMap\" + subdir(xx) + "\" + subdir(xx) + ".mat")
        tic
        clear X

        % load data
        tt = readcell(in + dirr + "\" + subdir(xx) + "\" + subdir(xx) + ".txt");
        ms2 = readcell(in + dirr + "\" + subdir(xx) + "\" + subdir(xx) + ".MS2.txt");       

        % delete empty rows        
        del = [];
        for i = 2:height(tt), if ~isnumeric(cell2mat(tt(i,end-15:end))), del = [del;i]; end; end
        tt(del,:) = [];

        del = [];
        for i = 2:height(ms2), if ~isnumeric(cell2mat(ms2(i,end-15:end))), del = [del;i]; end; end
        ms2(del,:) = [];        

        % filter FDR
        i1 = join(string(ms2(:,[3,4,6]))," ");
        i2 = join(string(tt(:,[2,4,6]))," ");

        for i = 2:length(i1) % skip i = 1 because this is the header
            r = find(strcmp(i1,i1(i)));            
            fdr = cell2mat(ms2(r,14));

            % delete peptide groups with only 1 PSM which fails FDR threshold
            if nnz(fdr > 0.01) == 1 && length(fdr) == 1
                s = find(strcmp(i2,i1(i)));

                if ~isempty(s)
                    tt(s,:) = [];
                    i2(s) = [];
                end

            % subtract abundances of PSM
            elseif nnz(fdr > 0.01) > 0 && length(fdr) > 1

                f_i = cell2mat(ms2(i,14));

                if f_i > 0.01                
                    s = find(strcmp(i2,i1(i)));                                     
    
                    if ~isempty(s)
                        tt(s,end-15:end) = num2cell(cell2mat(tt(s,end-15:end)) - cell2mat(ms2(i,end - 15:end)));
                    end
                end

            end
        
        end                
       
        % join relevant data
        tt = tt(:,[2,3,4,5,6,15:size(tt,2)]);
        
        % create initial data struct
        fld = string(tt(1,1:5));
        for i = 1:length(fld), X.(fld(i)) = string(tt(2:end,i));end
        X.Modifications(ismissing(X.Modifications)) = "N/A";
        
        % rename variants
        idx1 = find(contains(X.Proteins,"Variant"));
        lines = [];
        for i = 1:length(idx1)
            sp = split(X.Proteins(idx1(i)),';'); sp(cellfun('isempty',sp)) = [];
            if length(sp) == 1, continue; 
            else
                del = contains(sp,"Variant");
                sp(del) = [];
                X.Proteins(idx1(i)) = join(sp,";") + ";";
            end
                   
        end
    
        % create array of peptide abundances
        X.a = cell2mat(tt(2:end,6:end));

        % label batch
        X.batch = repmat(subdir(xx),[1,size(X.a,2)]);

        % get unambiguous identifier for a particular peptide
        X.id = X.Proteins + "&" + X.Sequence + "&" + X.ModSeq + "&" + X.Modifications + "&" + X.Charge;

        % save
        save(in + "DrugMap\" + subdir(xx) + "\" + subdir(xx) + ".mat",'X')
        toc
        disp(xx)        
%     end
end

% % loop to merge
dirr = "DrugMap";
fullfile = dir(in + dirr + "\");
subdir = string({fullfile.name})'; subdir(1:2) = []; 

load(in+dirr + "\" + subdir(1) + "\" + subdir(1) + ".mat");
L = load(in+dirr + "\" + subdir(2) + "\" + subdir(2) + ".mat"); L = L.X;
X.tid = union(X.id,L.id);
X.("a2") = nan(height(X.tid),size(X.a,2) + size(L.a,2)); 

[~,i1,i2] = intersect(X.tid, X.id);
X.a2(i1,1:size(X.a,2)) = X.a(i2,:);

[~,i1,i2] = intersect(X.tid, L.id);
X.a2(i1,size(X.a,2)+ 1:size(X.a2,2)) = L.a(i2,:);
X.batch = [X.batch,L.batch];
X.id = X.tid; 
X.a = X.a2; X = rmfield(X,["a2","tid","Sequence","Proteins","ModSeq","Modifications","Charge"]);

% now loop and collect all TMTs into one array
for x = 3:length(subdir)    
    L = load(in+dirr + "\" + subdir(x) + "\" + subdir(x) + ".mat"); L = L.X;
    X.tid = union(X.id,L.id);
    X.("a2") = nan(height(X.tid),size(X.a,2) + size(L.a,2)); 
    
    [~,i1,i2] = intersect(X.tid, X.id);
    X.a2(i1,1:size(X.a,2)) = X.a(i2,:);
    
    [~,i1,i2] = intersect(X.tid, L.id);
    X.a2(i1,size(X.a,2)+ 1:size(X.a2,2)) = L.a(i2,:);
    X.batch = [X.batch,L.batch];
    X.id = X.tid; 
    X.a = X.a2; X = rmfield(X,["a2","tid"]);    
    
    disp(height(X.a));disp(x)
end
