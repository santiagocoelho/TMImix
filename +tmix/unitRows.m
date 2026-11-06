function unitRows(d)
validateattributes(d,{'numeric'},{'real','finite','2d','ncols',3,'nonempty'});
assert(all(abs(sum(d.^2,2)-1)<1e-6),'Expected unit directions as N x 3 rows.');
end
