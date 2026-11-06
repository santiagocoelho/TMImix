classdef FitPlan
    % Map a box of independent unit coordinates onto physical parameters.
    % Fixed fractions stay fixed. Free fractions share only the remaining mass.
    properties (SetAccess=private)
        spec
        template
        freeFractions
        names
        lower
        upper
        remaining
        dimension
    end
    methods
        function obj=FitPlan(model,spec,fractionOrder)
            tmix.validate(model,spec,true);
            obj.spec=spec; obj.template=spec; f=model.fractions;
            if nargin==3
                assert(isequal(sort(fractionOrder),sort(f)), ...
                    'Fraction order must contain each active pool once.');
                f=fractionOrder;
            end
            vary=arrayfun(@(k) ~isscalar(spec.(k)),f);
            obj.freeFractions=f(vary); obj.remaining=1;
            for k=f(~vary), obj.remaining=obj.remaining-spec.(k); end
            lo=arrayfun(@(k) spec.(k)(1),f(vary));
            hi=arrayfun(@(k) spec.(k)(end),f(vary));
            % A face of the simplex may collapse to one point.
            if abs(sum(lo)-obj.remaining)<1e-12 || abs(sum(hi)-obj.remaining)<1e-12
                z=lo; if abs(sum(hi)-obj.remaining)<1e-12, z=hi; end
                for k=1:numel(z), obj.template.(obj.freeFractions(k))=z(k); end
                obj.freeFractions=strings(1,0); lo=[]; hi=[];
            end
            allNames=string(fieldnames(spec))';
            obj.names=allNames(~ismember(allNames,f) & arrayfun(@(k) ~isscalar(spec.(k)),allNames));
            obj.lower=lo; obj.upper=hi;
            obj.dimension=max(0,numel(obj.freeFractions)-1)+numel(obj.names);
        end

        function p=decode(obj,x)
            x=x(:)'; assert(numel(x)==obj.dimension && all(isfinite(x)),'Invalid fit coordinates.');
            x=min(1,max(0,x)); p=obj.template; left=obj.remaining; j=0;
            for i=1:numel(obj.freeFractions)
                if i==numel(obj.freeFractions), value=left;
                else
                    [lo,hi]=obj.interval(i,left); j=j+1; value=lo+(hi-lo)*x(j);
                end
                p.(obj.freeFractions(i))=max(0,min(1,value)); left=left-value;
            end
            for k=obj.names
                a=obj.spec.(k); j=j+1; p.(k)=a(1)+(a(2)-a(1))*x(j);
            end
        end

        function x=encode(obj,p)
            % Reconcile an incompatible start with the imposed bounds. On
            % matching truth this is the identity (up to roundoff).
            x=zeros(1,obj.dimension); left=obj.remaining; j=0;
            for i=1:numel(obj.freeFractions)-1
                [lo,hi]=obj.interval(i,left);
                value=min(hi,max(lo,p.(obj.freeFractions(i)))); j=j+1;
                if hi>lo, x(j)=(value-lo)/(hi-lo); end
                left=left-value;
            end
            for k=obj.names
                a=obj.spec.(k); j=j+1; x(j)=min(1,max(0,(p.(k)-a(1))/(a(2)-a(1))));
            end
        end
    end
    methods (Access=private)
        function [lo,hi]=interval(obj,i,left)
            lo=max(obj.lower(i),left-sum(obj.upper(i+1:end)));
            hi=max(lo,min(obj.upper(i),left-sum(obj.lower(i+1:end))));
        end
    end
end
