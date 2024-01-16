% get directories; these are just strings which shorten file paths and enable navigation of file structure
[in,out] = ld;


% pre-wrangle files so that we don't have to repeat in future
dirr = "DrugMap";
fullfile = dir(in + dirr + "\");
subdir = string({fullfile.name})'; subdir(1:2) = [];

for xx = 1:length(subdir)
    
    % only do this if we haven't made a .mat already    
    if ~isfile(in + "Data\" + subdir(xx) + "\" + subdir(xx) + ".mat")
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
        
            % subtract abundances of PSM with FDR > 0.01
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
        save(in + "Data\" + subdir(xx) + "\" + subdir(xx) + ".mat",'X')
        toc
        disp(xx)        
    end
end

% loop to merge
dirr = "Data";
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

% delete rows with missing id
idx = find(ismissing(X.id));
X.id(idx) = [];
X.a(idx,:) = [];

save("CDM.v.1.1.mat","X","-v7.3")

%% now wrangle the abundances and remove bad peptides

% replace 0 abundance with NaN
X.a(X.a == 0) = NaN;
str = strrep(X.batch,'.mat','');
str2 = str;

% load metadata
t = string(readcell(in+"DrugMap.metadata.xlsx"));

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
X.line.trt = repmat(["KB05","KB03","KB02","DMSO"],[1,sl(X.pep.a)/4]);

% remove empty TMT channels
idx = find(strcmp(X.line.name,"empty"));
X.pep.a(:,idx) = []; 
fld = fieldnames(X.line);for i =1 :length(fld), X.line.(fld{i})(idx) = []; end
X.line.batch = strrep(X.line.batch,'.txt'," ");

% 2D array into 3D
u = ["KB05","KB03","KB02","DMSO"];
X.pep.a2 = nan(height(X.pep.a),sl(X.pep.a)/4,4);
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
            splitty = split(X.pep.modification(i),';'); splitty = splitty(contains(splitty,"IADTB"));
            s = string(regexp(splitty,'(?<=IADTB of C[\(]).+(?=[\)])','match'));                                    
            for k = 1:length(s)
                sp3 = split(sp2(j),"Variant");
                sp4 = upper(strrep(sp3(1),'-',''));
                X.pep.oncovariant(i) = X.pep.oncovariant(i) + ";" + regexprep(upper(sp4),'-','');
                sp4 = split(sp3(2),'(');
                idx = find(strcmp(P.md.accession,sp4(1)));
    
                X.pep.gene(i) = X.pep.gene(i) + ";" + P.md.gene(idx);
                X.pep.protein(i) = X.pep.protein(i) + ";" + P.md.protein(idx);
                X.pep.accession(i) = X.pep.accession(i) + ";" + sp4(1);
    
                site = string(regexp(sp2(j),'(?<=[\(]).+(?=[\)])','match'));
                X.pep.cys(i) = X.pep.cys(i) + ";" + string(str2double(site) + 1 + str2double(s(k)));
                X.pep.gene_cys(i) = X.pep.gene_cys(i) + ";" + P.md.gene(idx) + " C" + string(str2double(site) + 1 + str2double(s(k)));
                X.pep.acc_cys(i) = X.pep.acc_cys(i) + ";" + P.md.accession(idx) + " C" + string(str2double(site) + 1 + str2double(s(k)));
            end

        % if not
        else            
            splitty = split(X.pep.modification(i),';'); splitty = splitty(contains(splitty,"IADTB"));
            s = string(regexp(splitty,'(?<=IADTB of C[\(]).+(?=[\)])','match'));            
            for k = 1:length(s)
                sp4 = split(sp2(j),'(');
                if any(contains(sp4,'_'))
                    sp5 = split(sp4(1),'_'); X.pep.seqvariant(i) = sp5(2) + "_" + sp5(3);                    
                    sp4 = sp5(2);                    
                end
                idx = find(strcmp(P.md.accession,sp4(1)));
        
                X.pep.gene(i) = X.pep.gene(i) + ";" + P.md.gene(idx);
                X.pep.protein(i) = X.pep.protein(i) + ";" + P.md.protein(idx);
                X.pep.accession(i) = X.pep.accession(i) + ";" + sp4(1);
        
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

% now delete rarely detected, oxidized peptides
% we assume that (most) methionine oxidations represent technical artifact
del = setdiff(del,del(kp));

fld = string(fieldnames(X.pep));
for i = 1:length(fld), X.pep.(fld(i))(del,:,:) = []; end

% we now calculate a unique peptide abundance for all peptides
% e.g., median the abundances of different oxoforms of same peptide

[u,row] = unique(X.pep.peptide);
X.pep.a2 = nan(length(u),sl(X.pep.a),4);
X.pep.e2 = nan(length(u),sl(X.pep.a),3);

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
load(in + "cysteine.ontology.mat")

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

%% further normalize engagement

% quantile normalize
for i = 1:3, X.pep.eq(:,:,i) = quantilenorm(X.pep.e(:,:,i)); end

% at this point, the engagement values are un-scaled and run from 0 to 1
% we will now collapse replicates and re-scale the engagement values 

[u,col] = unique(X.line.name);
lin = X.line.lineage(col);
bch = X.line.batch(col);


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
bch = X.line.batch(col);
bch = split(bch,'_'); bch = bch(1,:,1)'; ubch = unique(bch);
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
del = find(sum(isnan(qnt(:,:,2)),2) == sl(qnt));
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

    strs(find(idcs == i3(i))) = pos;
end

% assemble final tables

q2 = nan(length(idcs),sl(qnt),3);
for i = 1:length(idcs)
    q2(i,:,:) = median(qnt(strcmp(X.pep.acc_cys,X.pep.acc_cys(idcs(i))),:,:),1,"omitnan");
end

s2 = nan(length(idcs),sl(st),3);
for i = 1:length(idcs)
    s2(i,:,:) = median(st(strcmp(X.pep.acc_cys,X.pep.acc_cys(idcs(i))),:,:),1,"omitnan");
end

d2 = nan(length(idcs),sl(nd),3);
for i = 1:length(idcs)
    d2(i,:,:) = round(median(nd(strcmp(X.pep.acc_cys,X.pep.acc_cys(idcs(i))),:,:),1,"omitnan"));
end

d = [eraseBetween(X.pep.accession(idcs),1,1), ...
    eraseBetween(X.pep.gene(idcs),1,1), ...
    eraseBetween(X.pep.protein(idcs),1,1), ...
    X.pep.peptide(idcs), ...
    eraseBetween(X.pep.cys(idcs),1,1), ...
    strs];

for i = 1:height(d), for j = 1:sl(d), s = split(d(i,j),";"); if length(s) > 1, d(i,j) = strjoin(s(1:end - 1),";"); end; end; end

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
    for i = 1:sl(q2), ft = [ft,[q2(:,i,jj),s2(:,i,jj),d2(:,i,jj)]]; end
    
    ach_id(ismissing(ach_id)) = "N/A";
    fl = [];
    for i = 1:sl(q2), fl = [fl,[u(i) + " (" + ach_id(i) + ")" + " Median Engagement (%)", u(i) + " (" + ach_id(i) + ")" + " Standard Deviation", u(i) + " (" + ach_id(i) + ")" + " # of Replicates"]]; end
    t = table(ft);
    t = splitvars(t);
    t.Properties.VariableNames = fl;
    
    tb = [f,t];
    writetable(tb,out+sct(jj) + ".v.3.xlsx")
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
M = load("CCLE.mutations.mat"); M = M.X;

% get unique information for each cell line (cell line name, project Achilles ID, batch)
[u,col] = unique(X.line.name);
bch = X.line.batch(col);
bch = split(bch,'_'); bch = bch(:,1)'; ubch = unique(bch);
ach_id = X.line.DepMap_ID(col);
lin = X.line.lineage(col);
mutlines = unique(M.ModelID);

% subset DrugMap on cell lines which have recorded mutations in CCLE
[~,i1,i2] = intersect(mutlines,ach_id);
qnt2 = qnt(:,i2,:); 
bch = bch(i2); ach_id = ach_id(i2); u = u(i2);

gn = repmat("",[height(qnt2),1]);
for i = 1:height(qnt2)   
    s = split(X.pep.gene(i),';');
    s = unique(s(~cellfun(@isempty,s)));
    if length(s) == 1
        gn(i) = s;
    end   
end

acn = repmat("",[height(qnt2),1]);
for i = 1:height(qnt2)   
    s = split(X.pep.accession(i),';');
    s = unique(s(~cellfun(@isempty,s)));
    ss = [];
    for j = 1:length(s)
        sp = split(s(j),'-');
        ss = [ss;sp(1)];
    end
    s = unique(ss);
    if length(s) == 1
        acn(i) = s;
    end   
end


mut = false(height(X.pep.gene_cys),length(ach_id));
ugn = unique(gn); ugn(1) = [];
uacn = unique(acn); uacn(1) = [];
for i = 1:length(uacn)    
    if any(strcmp(M.UniprotID,uacn(i)))
        mtd = M.ModelID(strcmp(M.UniprotID,uacn(i)));
        if ~isempty(mtd)
            [~,i1,i2] = intersect(ach_id,mtd);
            if ~isempty(i1)            
                row = strcmp(acn,uacn(i));
                mut(row,i1) = 1;
            end        
        end
    end
    disp(i/length(uacn))
end

%bch = split(bch,'_'); bch = bch(1,:,1)'; ubatch = unique(bch);

ubatch = unique(bch);
db = [];for i = 1:length(ubatch),if nnz(strcmp(bch,ubatch(i))) < 3, db = [db;ubatch(i)];end;end

del = find(contains(bch,db));
mut(:,del) = []; qnt2(:,del,:) = [];
bch(del) = []; ach_id(del) = []; u(del) = [];


A.dat.qnt = qnt2;
A.dat.mutated = mut;
A.dat.gene = gn;
A.dat.accession = acn;
A.dat.gene_cys = X.pep.gene_cys;
A.dat.accession_cys = X.pep.acc_cys;

A.line.batch = bch;
A.line.DepMap_ID = ach_id;
A.line.stripped_cell_line_name = u';
X = A;save("mutations.x.CDM.v.1.1.mat",'X')


X.dat.new = repmat("",[height(X.dat.mutated),sl(X.dat.mutated)]);
X.dat.old = repmat("",[height(X.dat.mutated),sl(X.dat.mutated)]);
X.dat.mutation = repmat("",[height(X.dat.mutated),sl(X.dat.mutated)]);


ug = unique(X.dat.accession); ug(1) = [];
for i = 1:length(ug)
    row = find(strcmp(X.dat.accession,ug(i)));
    row2 = row(1);
    for j = 1:sl(X.dat.mutation)        
        if any(X.dat.mutated(row,j))
            r = find(strcmp(M.UniprotID,X.dat.accession(row2))&strcmp(M.ModelID,X.line.DepMap_ID(j)));
            if ~isempty(r)
                X.dat.new(row,j) = strjoin(M.new(r));
                X.dat.old(row,j) = strjoin(M.old(r));
                X.dat.mutation(row,j) = strjoin(M.ProteinChange(r));
            end        
        end
    end
    disp(i/length(ug))
end

u = unique(X.line.batch);
X.dat.delta = nan(height(X.dat.qnt),sl(X.dat.qnt),3);
for i = 1:length(u)
    f = find(strcmp(X.line.batch,u(i)));

    if length(f) >= 3
        
        for j = 1:length(f)
            i1 = f(j);
            i2 = setdiff(f,f(j));
            X.dat.delta(:,i1,:) = median(X.dat.qnt(:,i1,:) - X.dat.qnt(:,i2,:),2,'omitnan');
        end
    end   
end

save("mutations.CDM.v.1.2.mat",'X','-v7.3')

pos = str2double(M.Pos);
mp = M.UniprotID + " " + M.ProteinChange;
[~,ii] =unique(mp,'stable');
mp = mp(ii);
pos2 = pos(ii);
X.dat.pos = nan(length(X.dat.gene),sl(X.dat.mutated));

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

G = load(in+"genomic.coordinates.mat");G = G.X;

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

X.dat.missense = false(height(X.dat.gene_cys),sl(X.dat.qnt));
for i = 1:height(X.dat.missense), for j = 1:sl(X.dat.missense), if X.dat.mutated(i,j), if ~strcmp(X.dat.old(i,j),X.dat.new(i,j)), X.dat.missense(i,j) = 1; end;end;end; disp(i/length(X.dat.gene_cys)); end

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

X.dat.missense = false(height(X.dat.gene_cys),sl(X.dat.qnt));
for i = 1:height(X.dat.missense), for j = 1:sl(X.dat.missense), if X.dat.mutated(i,j), if ~strcmp(X.dat.old(i,j),X.dat.new(i,j)), X.dat.missense(i,j) = 1; end;end;end; disp(i/length(X.dat.gene_cys)); end

save("mutations.CDM.v.1.4.mat","X","-v7.3")
%% integrate cysteines with structural data

S = load("structural.db.mat"); S = S.X;
load("CDM.v.1.6.mat")

% find peptides with only one protein assignment (isoforms allowed)
kp = [];
for i = 1:length(X.pep.accession)
    sp = split(X.pep.accession(i),';');
    sp = sp(~cellfun(@isempty,sp)); 
    for j = 1:length(sp)
        if contains(sp(j),'-'), sp(j) = regexp(sp(j),'.+(?=[-])','match'); end
    end
    sp = unique(sp(~contains(sp,'-')));
    if length(sp) == 1, kp = [kp;i]; end
end

% subset our data table on these peptides
fld = string(fieldnames(X.pep));
for i = 1:length(fld), X.pep.(fld(i)) = X.pep.(fld(i))(kp,:,:); end
qnt2 = qnt(kp,:,:);

% eliminate noisy ligandability estimates which may driven by biology of 
% individual cell lines --> remove peptides w/ low detection
kp = find(X.pep.det >= 71);
fld = string(fieldnames(X.pep));
for i = 1:length(fld), X.pep.(fld(i)) = X.pep.(fld(i))(kp,:,:); end
qnt2 = qnt2(kp,:,:); 

fld = ["KB05","KB03","KB02"];
for i = 1:length(fld), S.(fld(i)) = nan(length(S.uninum),1); end
for i = 1:length(S.pdb)
    if ~isempty(S.gene{i})
        idx = find(contains(X.pep.gene_cys,";" + S.gene(i) + " C" + S.uninum(i) + ";"));    
        if ~isempty(idx)        
            idx = argMax(idx,X.pep.det); idx = idx(1);
            m = median(median(qnt2(idx,:,:),2,'omitnan'),1,'omitnan');
            for j = 1:length(fld), S.(fld(j))(i,1) = m(1,1,j); end
            
        end
    end
end


X = S; save("structural.db.v.2.mat",'X')
