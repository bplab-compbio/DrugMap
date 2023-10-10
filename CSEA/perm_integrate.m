function [p] = perm_integrate(null_dist,observed,sidedness)
% Integrates a null distribution of permutations to calculate p value
% Need to provide observed as the value in the right tail
    
    swtch = 0;

    % make histogram of null distribution
    h = histogram(null_dist);
    b = h.BinEdges;
    bb = [];

    % get midpoints of left- and right-endpoints in histogram; this forms
    % abscissae of interpolated gaussian
    for j = 1:length(b) - 1, bb = [bb;(b(j + 1) + b(j))/2]; end    

    % if we don't have enough abscissae, add a few
    if length(bb) == 2, swtch = 1; bb = [bb;2]; end
    if length(bb)  == 1, swtch = 2; bb = [bb;[1;2]]; end

    % get ordinates of null distribution
    v = h.Values;
    v = v/sum(v);

    % if no intersections in null distribution, add a bit of density
    if swtch == 1, v = [v,0.001]; end
    if swtch == 2, v = [v,0.001,0.001]; end

    % fit gaussian
    f = fit(bb,v','gauss1');
    fun = @(x,a,b,c) a*exp(-1*((x-b)/c).^2);


    % integrate either right tail or both tails
    if sidedness == 1
        tot = integral(@(x) fun(x,f.a1,f.b1,f.c1),-1*inf,inf);
        nom = integral(@(x) fun(x,f.a1,f.b1,f.c1),observed,inf);
        p = nom/tot;
    elseif sidedness == 2
        tot = integral(@(x) fun(x,f.a1,f.b1,f.c1),-1*inf,inf);
        nom = 2 * integral(@(x) fun(x,f.a1,f.b1,f.c1),observed,inf);
        p = nom/tot;        
    end
end
