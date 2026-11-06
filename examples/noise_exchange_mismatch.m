function result = noise_exchange_mismatch(c,exchangeTime,varargin)
% Fit the exchanging-tissue observations with impermeable sticks, ECS and spheres.
assert(isequal(c.parts,["stick","ecs","sphere"]),'Expected impermeable sticks, ECS and spheres.');
result = noise_exchange_population(c,exchangeTime,varargin{:});
end
