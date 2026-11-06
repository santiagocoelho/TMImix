function [axes,mass]=odf(axes,mass)
tmix.unitRows(axes);
validateattributes(mass,{'numeric'},{'vector','finite','real','nonnegative','numel',size(axes,1)});
mass=mass(:); assert(sum(mass)>0,'The fODF must have positive total mass.');
mass=mass/sum(mass);
end
