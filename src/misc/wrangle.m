% pre-wrangle files so that we don't have to repeat in future

% get base of directory
dirr = "DrugMap";

% get list of files
fullfile = dir(dirr + "\");

% get folder names
subdir = string({fullfile.name})'; subdir(1:2) = [];

% for each folder
for xx = 1:length(subdir)
    
    % only do this if we haven't made a .mat already    
    if ~isfile("Data\" + subdir(xx) + "\" + subdir(xx) + ".mat")

        % just to check time
        tic

        % clear last struct
        clear X
        
        % load data
        tt = readcell(dirr + "\" + subdir(xx) + "\" + subdir(xx) + ".txt");
        ms2 = readcell(dirr + "\" + subdir(xx) + "\" + subdir(xx) + ".MS2.txt");       
        
        % delete empty rows for MS3 data 
        del = []; for i = 2:height(tt), if ~isnumeric(cell2mat(tt(i,end-15:end))), del = [del;i]; end; end
        tt(del,:) = [];
        
        % delete empty rows for MS2 data
        del = []; for i = 2:height(ms2), if ~isnumeric(cell2mat(ms2(i,end-15:end))), del = [del;i]; end; end
        ms2(del,:) = [];        
        
        % get FDR
        i1 = join(string(ms2(:,[3,4,6]))," ");
        i2 = join(string(tt(:,[2,4,6]))," ");
        
        for i = 2:length(i1) % skip i = 1 because this is the header

            % find peptide
            r = find(strcmp(i1,i1(i)));            

            % get FDR
            fdr = cell2mat(ms2(r,14));
        
            % delete peptide groups with only 1 PSM which fails FDR threshold
            if nnz(fdr > 0.01) == 1 && length(fdr) == 1
                s = find(strcmp(i2,i1(i)));
        
                % if present in the MS3 sheet
                if ~isempty(s)

                    % delete from both
                    tt(s,:) = [];
                    i2(s) = [];
                end
        
            % subtract abundances of PSM with FDR > 0.01
            elseif nnz(fdr > 0.01) > 0 && length(fdr) > 1
        
                % get FDR
                f_i = cell2mat(ms2(i,14));
        
                if f_i > 0.01     

                    % find peptide with bad PSM in MS3 sheet
                    s = find(strcmp(i2,i1(i)));                                     
        
                    if ~isempty(s)
                        
                        % subtract
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

        % fill missing
        X.Modifications(ismissing(X.Modifications)) = "N/A";
        
        % rename variants
        v = find(contains(X.Proteins,"Variant"));
        
        % for each variant
        for i = 1:length(v)

            % parse the text
            sp = split(X.Proteins(v(i)),';'); sp(cellfun('isempty',sp)) = [];
            if length(sp) == 1, continue; 
            else

                % rename the variant
                del = contains(sp,"Variant");
                sp(del) = [];
                X.Proteins(v(i)) = join(sp,";") + ";";
            end
                   
        end
        
        % create array of peptide abundances
        X.a = cell2mat(tt(2:end,6:end));
        
        % label batch
        X.batch = repmat(subdir(xx),[1,size(X.a,2)]);
        
        % get unambiguous identifier for a particular peptide
        X.id = X.Proteins + "&" + X.Sequence + "&" + X.ModSeq + "&" + X.Modifications + "&" + X.Charge;
        
        % save
        save("Data\" + subdir(xx) + "\" + subdir(xx) + ".mat",'X')
        toc

        % display folder #
        disp(xx)        
    end
end

% now loop to merge the .mat we just created

% get base of directory
dirr = "Data";

% get list of files
fullfile = dir(dirr + "\");

% get list of folders
subdir = string({fullfile.name})'; subdir(1:2) = []; 

% load first TMT
load(dirr + "\" + subdir(1) + "\" + subdir(1) + ".mat");

% load second TMT into struct L
L = load(dirr + "\" + subdir(2) + "\" + subdir(2) + ".mat"); L = L.X;

% merge peptide id
X.tid = union(X.id,L.id);

% initialize new array for abundances
X.("a2") = nan(height(X.tid),size(X.a,2) + size(L.a,2)); 

% backfill data from first TMT
[~,i1,i2] = intersect(X.tid, X.id);
X.a2(i1,1:size(X.a,2)) = X.a(i2,:);

% backfill data from new TMT
[~,i1,i2] = intersect(X.tid, L.id);
X.a2(i1,size(X.a,2)+ 1:size(X.a2,2)) = L.a(i2,:);

% add batch information
X.batch = [X.batch,L.batch];

% replace old peptide id field with the concatenated version
X.id = X.tid; 
X.a = X.a2; X = rmfield(X,["a2","tid","Sequence","Proteins","ModSeq","Modifications","Charge"]);

% now loop and collect all TMTs into one array
for x = 3:length(subdir)    

    % load the xth TMT
    L = load(dirr + "\" + subdir(x) + "\" + subdir(x) + ".mat"); L = L.X;

    % union peptide id
    X.tid = union(X.id,L.id);

    % make new abundance array with space for new TMT
    X.("a2") = nan(height(X.tid),size(X.a,2) + size(L.a,2)); 
    
    % backfill new abundance array with original 
    [~,i1,i2] = intersect(X.tid, X.id);
    X.a2(i1,1:size(X.a,2)) = X.a(i2,:);
    
    % backfill new abundance array with new TMT
    [~,i1,i2] = intersect(X.tid, L.id);
    X.a2(i1,size(X.a,2)+ 1:size(X.a2,2)) = L.a(i2,:);

    % append batch info
    X.batch = [X.batch,L.batch];

    % update peptide ids
    X.id = X.tid; 

    % replace old arrays with updated arrays
    X.a = X.a2; X = rmfield(X,["a2","tid"]);    
    
    % display the height of the abundance array for fun
    disp(height(X.a));
    
    % display TMT number
    disp(x)
end

% find rows with missing id
idx = find(ismissing(X.id));

% delete
X.id(idx) = [];
X.a(idx,:) = [];

% save
save("CDM.v.1.1.mat","X","-v7.3")

%% now wrangle the abundances and remove bad peptides

% replace 0 abundance with NaN
X.a(X.a == 0) = NaN;


str = strrep(X.batch,'.mat','');
str2 = str;

% load metadata
t = string(readcell("DrugMap.metadata.xlsx"));

% organize columns of final table by name of cell line, collect metadata
t(1,:) = [];
u = unique(t(:,1));
str = regexprep(str,'.txt','');
for i = 1:length(u)
    i1 = strcmp(str,u(i));
    i2 = find(strcmp(t(:,1),u(i)));
    str2(i1) = t(i2,2);
end

% add metadata
X.line.name = str2;
X.line.batch = X.batch;
X = rmfield(X,'batch');
X.pep.a = X.a; X =rmfield(X,'a');
X.pep.id = X.id; X = rmfield(X,'id');
X.line.trt = repmat(["KB05","KB03","KB02","DMSO"],[1,size(X.pep.a,2)/4]);

% remove empty TMT channels
idx = find(strcmp(X.line.name,"empty"));
X.pep.a(:,idx) = []; 
fld = fieldnames(X.line);for i =1 :length(fld), X.line.(fld{i})(idx) = []; end
X.line.batch = strrep(X.line.batch,'.txt'," ");

% 2D array into 3D
u = ["KB05","KB03","KB02","DMSO"];
X.pep.a2 = nan(height(X.pep.a),size(X.pep.a,2)/4,4);
for i = 1:length(u), X.pep.a2(:,:,i) = X.pep.a(:,strcmp(X.line.trt,u(i))); end
X.pep.a = X.pep.a2; X.pep = rmfield(X.pep,"a2");

i1 = find(strcmp(X.line.trt,"DMSO"));
X.line.name = X.line.name(i1); X.line.batch = X.line.batch(i1);
X.line.trt = u;

% delete peptides without cysteine
sp = split(X.pep.id,'&');
idx = ~contains(sp(:,2),"C");
fld = fieldnames(X.pep); for i = 1:length(fld), X.pep.(fld{i})(idx,:,:) = []; end

save("CDM.v.1.2.mat","X","-v7.3")

%% annotate peptides with metadata 
% ... like gene names, protein names, oxidation state, oncogenic variants
P = load("fasta.mat");P=P.X;
sp = split(X.pep.id,"&");
X.pep.peptide = sp(:,2);
X.pep.annotated = sp(:,3);
X.pep.modification = sp(:,4);
X.pep.charge = sp(:,5);

% delete peptides which matched to the reverse (i.e. decoy) sequences during post-hoc search
idx = contains(X.pep.id,"rev_");
fld = fieldnames(X.pep);
for i = 1:length(fld), X.pep.(fld{i})(idx,:,:) = []; end

% delete peptides which are not labelled with IADTB (the chemical probe)
idx = ~contains(X.pep.id,"IADTB");
fld = fieldnames(X.pep);
for i = 1:length(fld), X.pep.(fld{i})(idx,:,:) = []; end

% create empty fields
X.pep.accession = repmat("",[length(X.pep.charge),1]);
X.pep.gene = repmat("",[length(X.pep.charge),1]);
X.pep.protein = repmat("",[length(X.pep.charge),1]);
X.pep.seqvariant = repmat("",[length(X.pep.accession),1]);
X.pep.oncovariant = repmat("",[length(X.pep.accession),1]);
X.pep.cys = repmat("",[length(X.pep.accession),1]);
X.pep.gene_cys = repmat("",[length(X.pep.accession),1]);
X.pep.acc_cys = repmat("",[length(X.pep.accession),1]);
X.pep.metox = false(length(X.pep.cys),1);
X.pep.metox(contains(X.pep.id,"Oxidation")) = 1;


% for each peptide, this loop is collecting metadata from fasta.mat 
% (fasta.mat was assembled by looping over all uniprot accessions via curl, i.e.
% curl -H "Accept: application/xml" "https://rest.uniprot.org/uniprotkb/P12345" --) 
sp = split(X.pep.id,"&");
for i = 1:height(sp)
    sp2 = split(sp(i,1),';');    
    sp2(cellfun(@isempty,sp2)) = [];
    for j = 1:length(sp2)

        % if peptide is variant from canonical fasta
        if contains(sp2,"Variant")
            l = split(X.pep.modification(i),';'); l = l(contains(l,"IADTB"));
            s = string(regexp(l,'(?<=IADTB of C[\(]).+(?=[\)])','match'));                                    
            for k = 1:length(s)
                v = split(sp2(j),"Variant");
                t = upper(strrep(v(1),'-',''));
                X.pep.oncovariant(i) = X.pep.oncovariant(i) + ";" + regexprep(upper(t),'-','');
                t = split(v(2),'(');
                idx = find(strcmp(P.md.accession,t(1)));
    
                X.pep.gene(i) = X.pep.gene(i) + ";" + P.md.gene(idx);
                X.pep.protein(i) = X.pep.protein(i) + ";" + P.md.protein(idx);
                X.pep.accession(i) = X.pep.accession(i) + ";" + t(1);
    
                site = string(regexp(sp2(j),'(?<=[\(]).+(?=[\)])','match'));
                X.pep.cys(i) = X.pep.cys(i) + ";" + string(str2double(site) + 1 + str2double(s(k)));
                X.pep.gene_cys(i) = X.pep.gene_cys(i) + ";" + P.md.gene(idx) + " C" + string(str2double(site) + 1 + str2double(s(k)));
                X.pep.acc_cys(i) = X.pep.acc_cys(i) + ";" + P.md.accession(idx) + " C" + string(str2double(site) + 1 + str2double(s(k)));
            end

        % if not
        else            
            l = split(X.pep.modification(i),';'); l = l(contains(l,"IADTB"));
            s = string(regexp(l,'(?<=IADTB of C[\(]).+(?=[\)])','match'));            
            for k = 1:length(s)
                t = split(sp2(j),'(');
                if any(contains(t,'_'))
                    v = split(t(1),'_'); X.pep.seqvariant(i) = v(2) + "_" + v(3);                    
                    t = v(2);                    
                end
                idx = find(strcmp(P.md.accession,t(1)));
        
                X.pep.gene(i) = X.pep.gene(i) + ";" + P.md.gene(idx);
                X.pep.protein(i) = X.pep.protein(i) + ";" + P.md.protein(idx);
                X.pep.accession(i) = X.pep.accession(i) + ";" + t(1);
        
                site = string(regexp(sp2(j),'(?<=[\(]).+(?=[\)])','match'));
                X.pep.cys(i) =  X.pep.cys(i) + ";" + string(str2double(site) + 1 + str2double(s(k)));
                X.pep.gene_cys(i) = X.pep.gene_cys(i) + ";" + P.md.gene(idx) + " C" + string(str2double(site) + 1 + str2double(s(k)));
                X.pep.acc_cys(i) = X.pep.acc_cys(i) + ";" + P.md.accession(idx) + " C" + string(str2double(site) + 1 + str2double(s(k)));
                idx = find(strcmp(P.md.accession,sp2(j)));            
            end
        end
    end

    % keep track of progress
    disp(i / height(sp))
end

% convert cell vectors to string vectors
fld = ["gene","protein","accession","cys","acc_cys","gene_cys"];
for i = 1:length(fld), X.pep.(fld(i)) = X.pep.(fld(i)) + ";"; end

% add a field which tells you how many TMT runs a peptide was detected in
X.pep.det = sum(~isnan(X.pep.a(:,:,4)),2);

save("CDM.v.1.3.mat","X","-v7.3")

%% where possible, add DepMap annotations 
% and corresponding metadata for cell lines

% load .mat which contains CCLE metadata
M=load("CCLE.metadata.mat");M=M.X;

% list out field names
fld = {'DepMap_ID','cell_line_name','stripped_cell_line_name','CCLE_Name','alias','COSMICID','sex','source','RRID','WTSI_Master_Cell_ID','sample_collection_site','primary_or_metastasis','primary_disease','Subtype','age','Sanger_Model_ID','depmap_public_comments','lineage','lineage_subtype','lineage_sub_subtype','lineage_molecular_subtype','default_growth_pattern','model_manipulation','model_manipulation_details','patient_id','parent_depmap_id','Cellosaurus_NCIt_disease','Cellosaurus_NCIt_id','Cellosaurus_issues'};

% this is where we add metadata
for j = 1:length(X.line.name)
    [~, i1, i2] = intersect(X.line.name(j),M.md.stripped_cell_line_name);
    if ~isempty(i1)
        for i = 1:length(fld), X.line.(fld{i})(j,1) = M.md.(fld{i})(i2); end
    end
end

% manually fix lineages of cell lines which are not in DepMap

% bile duct
X.line.lineage(contains(X.line.name,"ICC")) = "bile_duct";
X.line.lineage(strcmp(X.line.name,"OZ")) = "bile_duct";
X.line.lineage(strcmp(X.line.name,"MG984")) = "bile_duct";

% lung
X.line.lineage(contains(X.line.name,"MGH")) = "lung";

% breast
X.line.lineage(contains(X.line.name,"BRX")) = "breast";

% skin
X.line.lineage(strcmp(X.line.name,"COLO853")) = "skin";
X.line.lineage(strcmp(X.line.name,"WM793B")) = "skin";

% ovary
X.line.lineage(strcmp(X.line.name,"OVCA429")) = "ovary";

% skin
X.line.lineage(strcmp(X.line.name,"A375S2")) = "skin";
X.line.lineage(strcmp(X.line.name,"MEL182")) = "skin";
X.line.lineage(strcmp(X.line.name,"MEL167")) = "skin";

% brain/glioma neurosphere
X.line.lineage(strcmp(X.line.name,"MGG75")) = "brain";
X.line.lineage(strcmp(X.line.name,"MGG23")) = "brain";
X.line.lineage(strcmp(X.line.name,"MGG123")) = "brain";

% make sure all fields are tall
X.line.name = X.line.name';
X.line.batch = X.line.batch';

% save
save("CDM.v.1.4.mat","X","-v7.3")

%% calculate engagement

% d is a "damp," like a pseudocount used in RNA-seq analysis
% ... is meant to quench extremely low abundance fluctuations
d = 1000;

X.pep.e = (X.pep.a(:,:,4) + d)./(X.pep.a(:,:,1:3) + X.pep.a(:,:,4) + 2 * d);

%% eliminate peptides which ...
% were only rarely detected AND only present with methionine oxidations

% get a list of peptides
u = unique(X.pep.peptide);


% grow vector which contains indices of potentially bad peptides
del = [];

% for each peptide
for i = 1:length(u)
    row = find(strcmp(X.pep.peptide,u(i)));

    % if peptide is only detected w/ oxidized methionine, add to del
    if nnz(row) == 1, if contains(X.pep.annotated(row),"[15.995(M)]"), del = [del;row]; end; end
end

% this tells us whether the peptide under question was rarely detected
kp = [];
for i = 1:length(del), if X.pep.det(del(i)) >= 50, kp = [kp;i]; end; end

% find rarely detected, oxidized peptides
% we assume that (most) methionine oxidations represent technical artifact
del = setdiff(del,del(kp));

% delete
fld = string(fieldnames(X.pep));
for i = 1:length(fld), X.pep.(fld(i))(del,:,:) = []; end

% calculate a unique peptide abundance for all peptides
% e.g., median the abundances of different oxoforms of same peptide

% get unique peptides
[u,row] = unique(X.pep.peptide);

% initialize empty arrays
X.pep.a2 = nan(length(u),size(X.pep.a, 2),4);
X.pep.e2 = nan(length(u),size(X.pep.a, 2),3);

% for each peptide
for i = 1:length(u)

    % find all entries with a certain peptide
    j = strcmp(X.pep.peptide,u(i));

    % median
    X.pep.a2(i,:,:) = median(X.pep.a(j,:,:),1,'omitnan');    
    X.pep.e2(i,:,:) = median(X.pep.e(j,:,:),1,'omitnan');    

    % track progress
    disp(i / length(u))
end

% now trim the .pep struct
fld = setdiff(string(fieldnames(X.pep)),["a","a2","e","e2"]);
for i = 1:length(fld), X.pep.(fld(i)) = X.pep.(fld(i))(row,:,:); end
X.pep.a = X.pep.a2; X.pep.e = X.pep.e2; X.pep = rmfield(X.pep,["a2","e2"]); 

% recalculate peptide detection
X.pep.det = sum(~isnan(X.pep.a(:,:,4)),2);

% add a field which tells us how long a peptide is
X.pep.seqlen = cellfun(@length,X.pep.peptide);

% et voila!
save("CDM.v.1.5.mat","X","-v7.3")

%% add domain, class, pathway information
load("cysteine.ontology.mat")

[X.pep.domain, X.pep.class, X.pep.pathway] = deal(cell(height(X.pep.a),1));

fld = ["domain","class","pathway"];
for i = 1:height(X.pep.a)
    ac = split(X.pep.acc_cys(i),';');
    ac = ac(~cellfun(@isempty,ac));
    ac = unique(ac);
    for l = 1:length(ac)
        for j = 1:length(fld)
            for k = 1:length(S.(fld(j)).name)
                if any(strcmp(S.(fld(j)).cys{k},ac(l)))
                    X.pep.(fld(j)){i} = [X.pep.(fld(j)){i};S.(fld(j)).name(k)];
                end
            end
        end
    end
    disp(i/height(X.pep.a))
end

save("CDM.v.1.6.mat","X","-v7.3")

%% normalize engagement

% quantile normalize
for i = 1:3, X.pep.eq(:,:,i) = quantilenorm(X.pep.e(:,:,i)); end

% at this point, the engagement values are un-scaled and run from 0 to 1
% we will now collapse replicates and re-scale the engagement values 

[u,col] = unique(X.line.name);
lin = X.line.lineage(col);
b = X.line.batch(col);


% pre-allocate memory
qnt = nan(height(X.pep.e),length(u),3);
cc = nan(height(X.pep.a),length(u),3);


% calculate cv across replicates for all cell lines
for i = 1:length(u)
    cols = strcmp(X.line.name,u(i));    
    dat = X.pep.eq(:,cols,:);
    cc(:,i,:) = std(dat,0,2,'omitnan')./mean(dat,2,'omitnan');    
end
    
% use cv to filter out noisy measurements, calculate median engagement
% ... for each cell line
for i = 1:length(u)        

    % get replicates
    dat = X.pep.eq(:,strcmp(X.line.name,u(i)),:);

    % take median
    m = median(dat,2,'omitnan');

    % don't retain noisy values or peptides which were only detected once
    m(sum(~isnan(dat),2) < 2 | cc(:,i,:) > 0.1) = NaN;
    qnt(:,i,:) = m;
end

% re-scale
qnt = 200*(qnt - 0.5);

% this array was used for analysis

%% Generate table for DrugMap.net

% quantile-normalize
for i = 1:3, X.pep.eq(:,:,i) = quantilenorm(X.pep.e(:,:,i)); end
X.pep.eq = 200 * (X.pep.eq - 0.5);

% collect unique column information for each cell line
[u,col] = unique(X.line.name);
b = X.line.batch(col);
b = split(b,'_'); b = b(1,:,1)'; ub = unique(b);
ach_id = X.line.DepMap_ID(col);
lin = X.line.lineage(col);

% pre-allocate memory for different arrays
qnt = nan(height(X.pep.e),length(u),3);
st = nan(height(X.pep.e),length(u),3);
nd = nan(height(X.pep.e),length(u),3);
for i = 1:length(u), qnt(:,i,:) = median(X.pep.eq(:,strcmp(X.line.name,u(i)),:),2,"omitnan"); end
for i = 1:length(u), st(:,i,:) = std(X.pep.eq(:,strcmp(X.line.name,u(i)),:),0,2,"omitnan"); end
for i = 1:length(u), nd(:,i,:) = sum(~isnan(X.pep.eq(:,strcmp(X.line.name,u(i)),:)),2); end

% delete cysteines with unconfident detection (i.e. detected in less than 2 runs at the spectrometer)
del = find(sum(isnan(qnt(:,:,2)),2) == size(qnt,2));
[~,i1] = unique(X.pep.acc_cys);
i1 = setdiff(i1,del);

% get unique peptides for methionine oxidations
idcs = [];
for i = 1:length(i1)
    f = find(strcmp(X.pep.acc_cys,X.pep.acc_cys(i1(i))));

    if length(f) == 1
        idcs = [idcs;f];
    else
        m = find(~X.pep.metox(f));
        if length(m) == 1
            idcs = [idcs;f(m)];
        elseif length(m) > 1
            idcs = [idcs;f(m(1))];
        elseif isempty(m)
            m = find(X.pep.metox(f));
            if length(m) == 1
                idcs = [idcs;f(m)];
            elseif length(m) > 1
                idcs = [idcs;f(m(1))];
            end
        end
    end    

end

i2 = find(X.pep.metox);
i3 = intersect(i2,idcs);

mod = X.pep.modification(i3);
cys = X.pep.cys(i3);

strs = repmat("",[length(idcs),1]);
for i = 1:length(mod)
    
    s = split(cys(i),";");
    s = s(~cellfun(@isempty,s));
    c = str2double(s);
    c = min(c);

    s = split(mod(i),";");
    s = s(~cellfun(@isempty,s));
    nums = str2double(string(regexp(s,'(?<=[\(]).+(?=[\)])','match')));


    f = find(contains(s,"Oxidation of "));
    g = find(contains(s,"IADTB of "));

    dtb = min(nums(g));

    oxpos = nums(f);

    pos = [];
    if length(f) == 1
        pos = c + oxpos - dtb;    
        
    elseif nnz(f) > 1
        for j = 1:length(f)            
            pos = [pos;c + oxpos(j) - dtb];
        end
        
    end

    pos = eraseBetween(strrep(join(";M"+ pos)," ",""),1,1);

    strs(idcs == i3(i)) = pos;
end

% assemble final tables

q2 = nan(length(idcs),size(qnt,2),3);
for i = 1:length(idcs)
    q2(i,:,:) = median(qnt(strcmp(X.pep.acc_cys,X.pep.acc_cys(idcs(i))),:,:),1,"omitnan");
end

s2 = nan(length(idcs),size(st,2),3);
for i = 1:length(idcs)
    s2(i,:,:) = median(st(strcmp(X.pep.acc_cys,X.pep.acc_cys(idcs(i))),:,:),1,"omitnan");
end

d2 = nan(length(idcs),size(nd,2),3);
for i = 1:length(idcs)
    d2(i,:,:) = round(median(nd(strcmp(X.pep.acc_cys,X.pep.acc_cys(idcs(i))),:,:),1,"omitnan"));
end

d = [eraseBetween(X.pep.accession(idcs),1,1), ...
    eraseBetween(X.pep.gene(idcs),1,1), ...
    eraseBetween(X.pep.protein(idcs),1,1), ...
    X.pep.peptide(idcs), ...
    eraseBetween(X.pep.cys(idcs),1,1), ...
    strs];

for i = 1:height(d), for j = 1:size(d,2), s = split(d(i,j),";"); if length(s) > 1, d(i,j) = strjoin(s(1:end - 1),";"); end; end; end

f = splitvars(table(d));
f.Properties.VariableNames = ["Accession","Gene","Protein","Peptide","Cysteine","Methionine Oxidation"];

g2 = repmat("",[length(f.Gene),1]);
p2 = repmat("",[length(f.Gene),1]);
for i = 1:length(f.Gene)
    s1 = split(f.Gene(i),";");
    ss2 = split(f.Protein(i),";");
    if length(s1) > 1
        l = unique(s1);
        if length(s1) > length(l)
            for j = 1:length(l)
                if nnz(strcmp(s1,l(j))) > 1
                    ss2(strcmp(s1,l(j))) = ss2(strcmp(s1,l(j))) + "-" + string(1:nnz(strcmp(s1,l(j))))';
                    s1(strcmp(s1,l(j))) = s1(strcmp(s1,l(j))) + "-" + string(1:nnz(strcmp(s1,l(j))))';                    
                end
            end
        end
    end

    g2(i) = join(s1,";");
    p2(i) = join(ss2,";");
end

f.Gene = g2;
f.Protein = p2;

sct = ["KB05","KB03","KB02"];
for jj = 1:3
    ft = [];
    for i = 1:size(q2,2), ft = [ft,[q2(:,i,jj),s2(:,i,jj),d2(:,i,jj)]]; end
    
    ach_id(ismissing(ach_id)) = "N/A";
    fl = [];
    for i = 1:size(q2,2), fl = [fl,[u(i) + " (" + ach_id(i) + ")" + " Median Engagement (%)", u(i) + " (" + ach_id(i) + ")" + " Standard Deviation", u(i) + " (" + ach_id(i) + ")" + " # of Replicates"]]; end
    t = table(ft);
    t = splitvars(t);
    t.Properties.VariableNames = fl;
    
    tb = [f,t];
    writetable(tb,sct(jj) + ".xlsx")
    disp(jj)
end


X.pep = rmfield(X.pep,"eq");
fld = string(fieldnames(X.pep));
for i = 1:length(fld), X.pep.(fld(i)) = X.pep.(fld(i))(idcs,:,:); end
X.pep.e = q2;
X.pep.std = s2;
X.pep.detections = d2;
save("CDM.v.1.7.mat","X","-v7.3")

%% integrate CDM with mutations
% to run this code block, first generate the array which is the output of the block entitled "normalize engagement"

load("CCLE.mutations.mat");

% get unique information for each cell line (cell line name, project Achilles ID, batch)
[u,col] = unique(X.line.name);

% get batch information
b = X.line.batch(col); b = split(b,'_'); b = b(:,1)'; ub = unique(b);

% write new vector with ...

% achilles ID
ach_id = X.line.DepMap_ID(col);

% lineage
lin = X.line.lineage(col);

% get cell lines with profiled mutations in CCLE
m = unique(M.ModelID);

% intersect DrugMap which have recorded mutations in CCLE
[~,i1,i2] = intersect(m,ach_id);

% subset
q = qnt(:,i2,:); 

% also subset column information
b = b(i2); ach_id = ach_id(i2); u = u(i2);

% pre-allocate vector which will contain gene labels
gn = repmat("",[height(q),1]);

% for each peptide
for i = 1:height(q)   
    
    % split gene labels into all identifications output by search engine
    s = split(X.pep.gene(i),';');

    % take unique
    s = unique(s(~cellfun(@isempty,s)));

    % if there's only one gene
    if length(s) == 1

        % accept
        gn(i) = s;
    end   
end

% pre-allocate vector which will contain accession labels
acn = repmat("",[height(q),1]);

% for each peptide
for i = 1:height(q)   

    % split accessions into all identifications output by search engine
    s = split(X.pep.accession(i),';');

    % take unique
    s = unique(s(~cellfun(@isempty,s)));

    % for each unique accession
    t = [];
    for j = 1:length(s)

        % if alternative isoform, just take accession
        sp = split(s(j),'-');

        % append
        t = [t;sp(1)];
    end

    % take unique
    s = unique(t);

    % if only one unique accession (allowing for alternative isoforms)
    if length(s) == 1

        % accept
        acn(i) = s;
    end   
end

% initialize logical which indicates mutational status of a gene encoding a peptide
m = false(height(X.pep.gene_cys),length(ach_id));

% get unique UNP accessions
ua = unique(acn); 

% first one is empty
ua(1) = [];

% for each accession
for i = 1:length(ua)    

    % if any mutations recorded for this accession in CCLE
    if any(strcmp(M.UniprotID,ua(i)))

        % get model names
        id = M.ModelID(strcmp(M.UniprotID,ua(i)));

        % if any models (there must be by construction)
        if ~isempty(id)

            % intersect with models profiled in DrugMap
            [~,i1,i2] = intersect(ach_id,id);

            % if overlap is not empty
            if ~isempty(i1)          

                % find rows in initial matrix which contain accession ua(i)
                row = strcmp(acn,ua(i));

                % assign mutation
                m(row,i1) = 1;
            end        
        end
    end

    % display progress for the impatient
    disp(i/length(ua))
end

% get unique batches
ub = unique(b);

% don't consider batches which are represented by less than 3 cell lines
db = []; for i = 1:length(ub),if nnz(strcmp(b,ub(i))) < 3, db = [db;ub(i)];end;end

% find these under-represented batches
del = find(contains(b,db));

% also delete from engagement data
q(:,del,:) = [];

% also delete from mutation matrix
m(:,del) = []; 

% delete from batch string, achilles_id string, stripped_cell_line_name string
b(del) = []; ach_id(del) = []; u(del) = [];

% make new struct

clear X;

% assign engagement
X.dat.qnt = q;

% assign mutations
X.dat.mutated = m;

% assign gene
X.dat.gene = gn;

% assign accession
X.dat.accession = acn;

% assign gene + cysteine labels
X.dat.gene_cys = X.pep.gene_cys;

% assign accession + cysteine labels
X.dat.accession_cys = X.pep.acc_cys;

% assign batch label
X.line.batch = b;

% assign DepMap_ID (Project Achilles ID)
X.line.DepMap_ID = ach_id;

% assign name of cell line (i.e. UACC257)
X.line.stripped_cell_line_name = u';

% save
save("mutations.x.CDM.v.1.1.mat",'X')

% now add new fields which will store the exact amino acid variants
fld = ["new","old","mutation"];

% initialize empty strings
for i = 1:length(fld), X.dat.(fld(i)) = repmat("",[height(X.dat.mutated),size(X.dat.mutated,2)]); end

% get unique accessions
ua = unique(X.dat.accession); ua(1) = [];

% for each accessions
for i = 1:length(ua)

    % find peptides which have been assigned accession ua(i)
    a = find(strcmp(X.dat.accession,ua(i)));

    % get accession
    b = a(1);

    % for each cell line
    for j = 1:size(X.dat.mutation,2)        

        % if mutated 
        if any(X.dat.mutated(a,j))

            % find mutation in CCLE
            r = find(strcmp(M.UniprotID,X.dat.accession(b))&strcmp(M.ModelID,X.line.DepMap_ID(j)));

            % if it's actually there
            if ~isempty(r)

                % assign information
                X.dat.new(row,j) = strjoin(M.new(r));
                X.dat.old(row,j) = strjoin(M.old(r));
                X.dat.mutation(row,j) = strjoin(M.ProteinChange(r));
            end        
        end
    end

    % display progress for the impatient
    disp(i/length(ua))
end

% get unique batch labels
u = unique(X.line.batch);

% initialize 
X.dat.delta = nan(height(X.dat.qnt),size(X.dat.qnt,2),3);

% for each batch
for i = 1:length(u)

    % find the cell lines in batch u(i)
    f = find(strcmp(X.line.batch,u(i)));

    % if at least 3
    if length(f) >= 3
        
        % for each cell line
        for j = 1:length(f)

            % get cell line f(j)
            i1 = f(j);

            % get batch mates
            i2 = setdiff(f,f(j));

            % calculate local difference between individual cell line and
            % its batch-mates
            X.dat.delta(:,i1,:) = median(X.dat.qnt(:,i1,:) - X.dat.qnt(:,i2,:),2,'omitnan');
        end
    end   
end

% save
save("mutations.CDM.v.1.2.mat",'X','-v7.3')

pos = str2double(M.Pos);
mp = M.UniprotID + " " + M.ProteinChange;
[~,ii] =unique(mp,'stable');
mp = mp(ii);
pos2 = pos(ii);
X.dat.pos = nan(length(X.dat.gene),size(X.dat.mutated,2));

ac = unique(X.dat.accession); ac(1) = [];
for i = 1:length(ac)
    ii = find(strcmp(X.dat.accession,ac(i)));
    ij = ii(1);
    f = find(X.dat.mutated(ij,:));
    if ~isempty(f)
        uu = X.dat.mutation(ij,f);
        for j = 1:length(f)
            ky = X.dat.accession(ij) + " " + uu(j);
            g = find(strcmp(mp,ky));
            if ~isempty(g)
                X.dat.pos(ii,f(j)) = pos2(g);
            end
        end
        disp(i/length(ac))
    end
end

save("mutations.CDM.v.1.3.mat",'X','-v7.3')

% load gene coordinates
G = load("genomic.coordinates.mat");G = G.X;

% pre-allocate array which will store chromosome name of each gene
X.dat.chromosome_name = nan(length(X.dat.gene),1);

% for each gene
for i = 1:height(X.dat.gene)

    % find gene i in G
    f = find(strcmp(G.dat.hgnc_symbol,X.dat.gene(i)));

    % if exists
    if ~isempty(f)

        % if only one exists
        if length(f) == 1            
            try
                
                % assign X --> 23
                if strcmp(G.dat.chromosome_name(f),"X")
                    X.dat.chromosome_name(i) = 23;

                % assign Y --> 24
                elseif strcmp(G.dat.chromosome_name(f),"Y")
                    X.dat.chromosome_name(i) = 24;

                % assign 1:22 --> 1:22
                else
                    X.dat.chromosome_name(i) = str2double(G.dat.chromosome_name(f));
                end                
            catch
            end
        end
    end
end

X.dat.gene_cys1 = repmat("",[length(X.dat.gene_cys),1]);
for i = 1:length(X.dat.gene_cys)
    sp = split(X.dat.gene_cys(i),';'); sp = sp(2);
    X.dat.gene_cys1(i) = sp;
end

X.dat.missense = false(height(X.dat.gene_cys),size(X.dat.qnt,2));
for i = 1:height(X.dat.missense), for j = 1:size(X.dat.missense,2), if X.dat.mutated(i,j), if ~strcmp(X.dat.old(i,j),X.dat.new(i,j)), X.dat.missense(i,j) = 1; end;end;end; disp(i/length(X.dat.gene_cys)); end

X.dat.chromosome_name = nan(length(X.dat.gene),1);
for i = 1:height(X.dat.gene)
    f = find(strcmp(G.dat.hgnc_symbol,X.dat.gene(i)));
    if ~isempty(f)
        if length(f) == 1            
            try
                if strcmp(G.dat.chromosome_name(f),"X")
                    X.dat.chromosome_name(i) = 23;
                elseif strcmp(G.dat.chromosome_name(f),"Y")
                    X.dat.chromosome_name(i) = 24;
                else
                    X.dat.chromosome_name(i) = str2double(G.dat.chromosome_name(f));
                end                
            catch
            end
        end
    end
end

X.dat.gene_cys1 = repmat("",[length(X.dat.gene_cys),1]);
for i = 1:length(X.dat.gene_cys)
    sp = split(X.dat.gene_cys(i),';'); sp = sp(2);
    X.dat.gene_cys1(i) = sp;
end

X.dat.missense = false(height(X.dat.gene_cys),size(X.dat.qnt,2));
for i = 1:height(X.dat.missense), for j = 1:size(X.dat.missense,2), if X.dat.mutated(i,j), if ~strcmp(X.dat.old(i,j),X.dat.new(i,j)), X.dat.missense(i,j) = 1; end;end;end; disp(i/length(X.dat.gene_cys)); end

save("mutations.CDM.v.1.4.mat","X","-v7.3")
%% integrate cysteines with structural data

S = load("structural.db.mat"); S = S.X;
load("CDM.v.1.6.mat")

% find peptides with only one protein assignment (isoforms allowed)
kp = [];
for i = 1:length(X.pep.accession)

    % split peptide into all possible UNP accession attributions
    sp = split(X.pep.accession(i),';');

    % delete empty entries
    sp = sp(~cellfun(@isempty,sp)); 

    % for each accession
    for j = 1:length(sp)

        % if it has a "-" character, delete it for now
        if contains(sp(j),'-'), sp(j) = regexp(sp(j),'.+(?=[-])','match'); end
    end

    % get unique accessions
    sp = unique(sp(~contains(sp,'-')));

    % if we only have one unique accession
    if length(sp) == 1, kp = [kp;i]; end
end

% get field names
fld = string(fieldnames(X.pep));

% subset our struct on these peptides
for i = 1:length(fld), X.pep.(fld(i)) = X.pep.(fld(i))(kp,:,:); end

q = qnt(kp,:,:);

% get field names
fld = string(fieldnames(X.pep));

% find noisy ligandability estimates--maybe driven by biology of individual cell lines
kp = find(X.pep.det >= 71);

% delete 
for i = 1:length(fld), X.pep.(fld(i)) = X.pep.(fld(i))(kp,:,:); end

% temporary array with filtered engagement data
q = q(kp,:,:); 

% for each probe
fld = ["KB05","KB03","KB02"];

% initialize empty engagement array
for i = 1:length(fld), S.(fld(i)) = nan(length(S.uninum),1); end

% for each PDB
for i = 1:length(S.pdb)

    % if gene name present
    if ~isempty(S.gene{i})

        % find corresponding cys in DrugMap data
        idx = find(contains(X.pep.gene_cys,";" + S.gene(i) + " C" + S.uninum(i) + ";"));    

        % if we have it
        if ~isempty(idx)        

            % find the version that is detected most frequently
            idx = argMax(idx,X.pep.det); idx = idx(1);
            
            % get median engagement
            m = median(median(q(idx,:,:),2,'omitnan'),1,'omitnan');

            % assign to our structural database
            for j = 1:length(fld), S.(fld(j))(i,1) = m(1,1,j); end
            
        end
    end
end

% save
X = S; save("structural.db.v.2.mat",'X')
