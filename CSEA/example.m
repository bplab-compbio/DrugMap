%% First, let's calculate engagement over our pre-wrangled .mat (see https://drugmap.net) and find cysteines to analyze.
%% This code block is not necessary to run CSEA but demonstrates and example of how we analyzed our data.

% quantile normalize
for i = 1:3, X.pep.eq(:,:,i) = quantilenorm(X.pep.e(:,:,i)); end

% get unique identifiers (cell line name, project Achilles ID) for each cell line
[u,col] = unique(X.line.name);
id = X.line.DepMap_ID(col);
lin = X.line.lineage(col);

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

%% find heterogeneous cysteines

% perform "batch-aware" normalization
b = X.line.batch(col); b = split(b,'_'); b = b(:,1);

% get vector of unique batches
ub = unique(b);

% get vector of unique achilles_id
id = X.line.DepMap_ID(col);

% for each probe, initialize array which we will fill with "batch-aware" data
fld = ["KB05","KB03","KB02"];
for i = 1:3, X.het.diff.(fld(i)) = []; X.het.std.(fld(i)) = []; end

% for each probe
for j = 1:3
    
    % for each batch
    for i = 1:length(ub)

        % get engagement data
        d = qnt(:,strcmp(b,ub(i)),j);

        % if at least two cell lines
        if size(d,2) > 2       
            
            % calculate std
            X.het.std.(fld(j)) = [X.het.std.(fld(j)),std(d,0,2,'omitnan')];

            % for each cell line
            for k = 1:size(d,2)

                % subtract batch median
                d2 = d(:,k) - median(d,2,'omitnan');

                % append to new array
                X.het.diff.(fld(j)) = [X.het.diff.(fld(j)),d2];
            end
        end
    end
end

% create struct with replicate-compressed representation of DrugMap (KB03)
T.e = qnt(:,:,2); T.line = lin'; T.batch = b'; T.name = u'; T.achilles_id = id';

% initialize data structure for batch-aware data
flds = ["e","line","batch","name","achilles_id"];
for i = 1:length(flds), S.(flds(i)) = []; end

clear j;

% for each unique batch
for i = 1:length(ub)

    % make a temporary struct with batch information
    for j = 1:length(flds)
        R.(flds(j)) = T.(flds(j))(:,strcmp(b,ub(i)));
    end

    % if we have at least two cell lines in batch
    if size(R.e,2) > 2  

        % for each cell line
        for k = 1:size(R.e,2)

            % append line/batch/name/achilles_id to new struct
            for l = 1:length(flds), S.(flds(l)) = [S.(flds(l));R.(flds(l))(k)]; end
        end
    end
end

% add new information to original struct
for i = 1:length(flds), X.het.(flds(i)) = S.(flds(i)); end

% NaN-out the 0 std cysteines (means low detection, not good for this analysis)
for i = 1:3, X.het.std.(fld(i))(X.het.std.(fld(i)) == 0) = NaN; end

% count number of times a cysteine has high intra-batch std
for i = 1:3, X.het.nnz.(fld(i)) = sum(X.het.std.(fld(i)) >= 10,2,'omitnan'); end

% rank order by this count
for i = 1:3, X.het.rank.(fld(i)) = 1 - tiedrank(X.het.nnz.(fld(i)))/length(find(~isnan(X.het.nnz.(fld(i))))); end

% calculate global std for each std, across DrugMap
for i = 1:3, X.het.globalstd.(fld(i)) = std(qnt(:,:,i),0,2,'omitnan'); X.het.globalstd.(fld(i))(X.het.globalstd.(fld(i)) == 0) = NaN; end

% rank-order these std
for i = 1:3, X.het.globalrank.(fld(i)) = 1 - tiedrank(X.het.globalstd.(fld(i)))/length(find(~isnan(X.het.globalstd.(fld(i))))); end

% get detection rate for each probe
for i = 1:3, X.het.det.(fld(i)) = sum(~isnan(qnt(:,:,i)),2)./size(qnt,2); end

% get global median ligandability for each cysteine
m = median(qnt(:,:,2),2,'omitnan');

% this array is a logical which identifies the heterogeneous cysteines
row = X.het.rank.KB03 < 0.05 & X.het.globalrank.KB03 < 0.05 & X.het.det.KB03 >= 0.5;

%% Run CSEA! We will find terms which are enriched among heterogeneous cysteines
% recall that, since CSEA is a stochastic algorithm, p-value estimates will not reproduce exactly across runs
C = load("CSEA.repository.mat"); C = C.X;

% get accession + cysteine (i.e. O60701 C276) for the heterogeneous cysteines
a = find(row);
t = []; for i = 1:height(a), s = split(X.pep.acc_cys(a(i)),';'); s = s(~cellfun(@isempty,s)); t = unique([t;unique(s)]); end

% get accession + cysteine (i.e. O60701 C276) for the non-heterogeneous cysteines
b = setdiff(find(X.het.rank.KB03 > 0.25 & X.het.globalrank.KB03 > 0.25 & X.het.det.KB03 >= 0.5),a);
n = []; for i = 1:height(b), s = split(X.pep.acc_cys(i),';'); s = s(~cellfun(@isempty,s)); n = unique([n;unique(s)]); end

% get cysteine sets of interest
r = find(contains(C.set.type,"molecular")|contains(C.set.name,"zhang_et"));

% run CSEA
res = csea(n,t,100,C.set.name(r),C.set.cys(r));

% sort the output
res = sortrows(res,2,"ascend");

% plot the output
x = fliplr(1:height(res));
y = -log10(res.("FDR"));
scatter(x,y,20,'k','filled');

% find the redox sets, from Zhang, J., Simpson, C. M., Berner, J., Chong, H. B., Fang, J., Ordulu, Z., ... & Bar-Peled, L. (2023). Systematic identification of anticancer drug targets reveals a nucleus-to-mitochondria ROS-sensing pathway. Cell, 186(11), 2361-2379.
f = find(contains(res.("set name"),"zhang") & y > 10);
hold on

% plot them
scatter(x(f),y(f),20,[1,0,0],"filled","markeredgecolor",[1,1,0])

% label axes
ff('signatures',"-log_{10}(FDR)",'',12,'tex','out')

% set axis limits
xlim([-1,175]);ylim([0,180])

% set axis ticks
set(gca,"ytick",[0:20:180]);set(gca,'xtick',[0:25:175]);

% set tick length
set(gca,'ticklength',[0.025,0.025]);

% make axis tight
axis tight
