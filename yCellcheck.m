function [nrep, nvec, dataWrd] = yCellcheck(yCell, nvar)

%  Last modified 14 July 2015

if ~iscell(yCell)
    error('YCELL is not a cell array.');
end
errwrd = 0;

nvec    = zeros(nvar,1);
dataWrd = zeros(nvar,1);
ysize   = zeros(nvar,1);
nobs    = 0;
%  find number of replications for first non-empty cell
for ivar=1:nvar
    if ~isempty(yCell{ivar})
        yStructi = yCell{ivar};
        nrep = size(yStructi.y,2);
        break
    end
end
%  loop through variables
for ivar=1:nvar
    if ~isempty(yCell{ivar}) && ~isempty(yCell{ivar}.y)
        dataWrd(ivar) = 1;
        yStructi = yCell{ivar};
        if ~isstruct(yStructi)
            warning(['YCELL{',num2str(ivar),'} is not a struct object.']);
            errwrd = 1;
        end
        if ~isfield(yStructi,'argvals')
            warning(['argvals is not a field for YCELL{', ...
                num2str(ivar),'}.']);
            errwrd = 1;
        end
        ni = length(yStructi.argvals);
        nvec(ivar) = ni;
        if ~isfield(yStructi,'y')
            warning(['y is not a field for YCELL{',num2str(ivar),'}.']);
            errwrd = 1;
        end
        ysizei = size(yStructi.y);
        if length(ysizei) > 2
            warning(['More than two dimensions for y in YCELL{', ...
                num2str(ivar),'}.']);
            errwrd = 1;
        else
            ysize(ivar) = ysizei(1);
        end
        %  set up and check NREP
        nrepi = ysizei(2);
        if nrepi ~= nrep
            warning('Second dimensions of YStruct.y are not equal.');
            errwrd = 1;
        end
        nobs = nobs + 1;
        if ni ~= ysizei(1)
            warning(['Length of ARGVALS and first dimension of Y', ...
                'are not equal.']);
            errwrd = 1;
        end
    else
        dataWrd(ivar) = 0;
    end
end

if nobs == 0
    warning('No variables have observations.');
    errwrd = 1;
end

if errwrd
    error('One or more terminal error encountered in YCELL.');
end

