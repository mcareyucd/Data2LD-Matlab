function coefCellnew = BAwtvec2cell(thetavec, coefCell)
%  Places estimated weight coefficients only in THETAVEC into pre-existing 
%    cell array COEFCELL.

%  Last modified 10 March 2016

if nargin < 2
    error('Number of arguments is less than 2.');
end

if ~iscell(coefCell)
    error('coefCELL is not a cell array.');
end

ncoef  = length(coefCell);
ntheta = length(thetavec);
coefCellnew = cell(ncoef,1);
for icoef=1:ncoef
    coefStructi = coefCell{icoef};
    if any(coefStructi.estimate)
        coefCellnew{icoef} = coefStructi;
        index = coefStructi.index;
        if max(index) > ntheta
            error('Coefficient index out of range.');
        end
        coefStructi.parvec = thetavec(index);
    end
    coefCellnew{icoef} = coefStructi;
end
