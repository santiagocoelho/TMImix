function validate(model,p,isSpec)
% The model schema also catches misspelled constraints.
assert(isscalar(p) && isequal(sort(fieldnames(p)),sort(fieldnames(model.defaults()))), ...
    'Parameter fields must match model.defaults() exactly.');
for k=string(fieldnames(p))'
    v=p.(k); isFrac=ismember(k,model.fractions);
    assert(isnumeric(v) && isreal(v) && isrow(v) && ...
        (isscalar(v) || (isSpec && numel(v)==2 && v(1)<v(2))) && ...
        all(isfinite(v) | (k=="tau" && isscalar(v) && v==Inf)), ...
        '%s must be a scalar or a finite increasing fit interval.',k);
    if isFrac || k=="p2"
        assert(all(v>=0 & v<=1),'%s must be in [0,1].',k);
    elseif k=="DePerp" || k=="DePar"
        assert(all(v>=0),'%s must be nonnegative.',k);
    else
        assert(all(v>0),'%s must be positive.',k);
    end
end
lo=0; hi=0;
for k=model.fractions, a=p.(k); lo=lo+a(1); hi=hi+a(end); end
assert(lo<=1+1e-12 && hi>=1-1e-12,'Fraction bounds must admit a sum of one.');
end
