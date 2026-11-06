function draw = tissueSampler(model,bounds,laws,fixed)
% Compile a sampler for one population. Each call to draw() returns one tissue.
% Laws are scalars, [lo hi] intervals, or @() draws in physical units.
% Optional fractionAlpha follows model.parts order and sets a joint Dirichlet law.
if nargin < 4, fixed = struct(); end
assert(isstruct(laws) && isscalar(laws),'Population laws must be a scalar struct.');
names = string(fieldnames(bounds))';
specified = string(fieldnames(laws))';
assert(all(ismember(specified,[names,"fractionAlpha"])), ...
    'Unknown parameter in sampling distribution.');
assert(all(ismember(string(fieldnames(fixed)),names)), 'Unknown fixed parameter.');
fractions = model.fractions;
% Explicit individual fraction laws come before unspecified residual pools.
first = specified(ismember(specified,fractions));
order = [first, fractions(~ismember(fractions,first))];
spec = bounds;
custom = struct();
for name = names
    if isfield(fixed,name)
        spec.(name) = fixed.(name);
    elseif isfield(laws,name)
        law = laws.(name);
        if isa(law,'function_handle')
            custom.(name) = law;
            if ismember(name,fractions), spec.(name) = [0 1]; end
        else
            spec.(name) = law;
        end
    end
end
plan = tmix.FitPlan(model,spec,order);
customNames = string(fieldnames(custom))';
customFractions = customNames(ismember(customNames,fractions));
alpha = [];
if isfield(laws,'fractionAlpha')
    alpha = laws.fractionAlpha;
    validateattributes(alpha,{'numeric'},{'vector','real','finite','positive','numel',numel(fractions)});
    alpha = reshape(alpha,1,[]);
end
draw = @sample;

    function p = sample()
        % Start with uniform free coordinates in the numeric bounds. Rebuild
        % the fraction plan only when individual stochastic fractions require it.
        localSpec = spec;
        localPlan = plan;
        if ~isempty(customFractions)
            for name = customFractions
                value = drawScalar(custom.(name),name);
                assert(value>=0 && value<=1,'Fraction draw for %s must lie in [0,1].',name);
                localSpec.(name) = value;
            end
            % Incompatible user-drawn fractions must be reported, never clipped.
            localPlan = tmix.FitPlan(model,localSpec,order);
        end
        p = localPlan.decode(rand(1,localPlan.dimension));
        if ~isempty(alpha)
            p = dirichletFractions(p,localSpec,fractions,alpha,localPlan);
        end
        for name = customNames(~ismember(customNames,fractions))
            p.(name) = drawScalar(custom.(name),name);
        end
        if ~isempty(customNames), tmix.validate(model,p,false); end
    end
end

function x = drawScalar(fun,name)
x = fun();
assert(isnumeric(x) && isscalar(x) && isreal(x) && ~isnan(x), ...
    'Sampling function for %s must return one real numeric value.',name);
end

function p = dirichletFractions(p,spec,fractions,alpha,plan)
% Conditional on fixed fractions, the free shares of the remaining mass are
% Dirichlet. Additional physical interval bounds truncate this joint law.
free = fractions(ismember(fractions,plan.freeFractions));
if numel(free)<2, return; end
a = alpha(ismember(fractions,free));
lo = arrayfun(@(k) spec.(k)(1),free);
hi = arrayfun(@(k) spec.(k)(2),free);
for attempt = 1:10000
    g = randg(a);
    f = plan.remaining*g/sum(g);
    if all(f>=lo & f<=hi)
        for i = 1:numel(free), p.(free(i)) = f(i); end
        return
    end
end
error('TMImix:FractionSampling','Dirichlet draws missed the fraction bounds 10000 times. Check the bounds or concentrations.');
end
