function thetavecnew = BAwtcell2vec(modelCell, coefCell)
%  Extracts the ESTIMATED weight coefficients only from MODELCELL
%  and assembles them into a vector.

%  Last revised 10 March 2016

if nargin < 1
    error('Number of arguments is less than one.');
end

if ~iscell(modelCell)
    error('MODELCELL is not a cell array.');
end

[coefCell, ntheta] = coefcheck(coefCell);

nvar = length(modelCell);
thetavecnew = zeros(ntheta,1);
for ivar=1:nvar
    modelStructi = modelCell{ivar};
    for iw=1:modelStructi.nallXterm
        modelStructiw = modelStructi.XCell{iw};
        ncoefw    = modelStructiw.ncoef;
        coefStrw  = coefCell{ncoefw};
        if coefStrw.estimate
            coefw = coefStrw.parvec;
            indw  = coefStrw.index;
            thetavecnew(indw) = coefw;
        end
    end
    if modelStructi.nallFterm > 0
        for jforce=1:modelStructi.nallFterm
            modelStructj = modelStructi.FCell{jforce};
            ncoefj    = modelStructj.ncoef;
            coefStrj  = coefCell{ncoefj};
            if coefStrj.estimate
                coefj = coefStrj.parvec;
                indj  = coefStrj.index;
                thetavecnew(indj) = coefj;
            end
        end
    end
end

