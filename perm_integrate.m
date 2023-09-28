function [p] = perm_integrate(null_dist,observed,sidedness)
% Integrates a null distribution of permutations to calculate p value
% Need to provide observed for lower limit in integral
    
    swtch = 0;
    h = histogram(null_dist);
    b = h.BinEdges;
    bb = [];
    for j = 1:length(b) - 1, bb = [bb;(b(j + 1) + b(j))/2]; end
    if length(bb) == 2, swtch = 1; bb = [bb;2]; end
    if length(bb)  == 1, swtch = 2; bb = [bb;[1;2]]; end
    v = h.Values;
    v = v/sum(v);
    if swtch == 1, v = [v,0.001]; end
    if swtch == 2, v = [v,0.001,0.001]; end
    f = fit(bb,v','gauss1');
    fun = @(x,a,b,c) a*exp(-1*((x-b)/c).^2);

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
