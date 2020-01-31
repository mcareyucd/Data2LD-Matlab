function modelCellnew = modelcheck(modelCell, coefCell)
% Check the modelCell structure, requiring at a minimum:
%   cellarray with a member for each variable, containing:
%     a cell array of length equal number of homogeneous terms present.  
%     Each cell has:
%       variable index i and/or tag
%       order of derivative in the  term
%       fdPar object for coefficient function
%   cell array of length either number of forcing functions or empty.  
%     If not empty, cells have:
%       fd object for forcing function(s)
%       fdPar object for coefficient function

%  Last modified 22 January 2017

ncoef = length(coefCell);

if isempty(modelCell)
    error('Argument MODELCELL is empty.')
end
if ~iscell(modelCell)
    error('Argument MODELCELL is not a cell object.');
end

%  check one-dimensionality of cell array

dim = size(modelCell);
if dim(1) > 1 && dim(2) > 1
    error('Argument MODELCELL is not a one-dimensional cell array.');
end

%  number of variables is length of cell array

nvar = length(modelCell);

%  check that each cell is a struct object

errwrd = 0;
for ivar=1:nvar
    modelStructi = modelCell{ivar};
    if ~isstruct(modelStructi)
        warning(['Object in modelCell{',num2str(ivar), ...
            '} is not a struct object.']);
        errwrd = 1;
    end
end
if errwrd
    error('one or more ModelCells do not contain struct objects.');
end

%  check that each cell has the necessary fields or, if not, 
%  set the fields to default values

for ivar=1:nvar
    modelStructi = modelCell{ivar};
    %  if order field not present, it defaults to 1
    if ~isfield(modelStructi, 'order')
        modelStructi.order = 1;
    end
    %  if XCell field not present, it defaults to {};
    if ~isfield(modelStructi,'XCell') || isempty(modelStructi.XCell)
        modelStructi.Xcell = {};
    end
    %  if FCell field not present, it defaults to {};
    if ~isfield(modelStructi,'FCell') || isempty(modelStructi.FCell)
        modelStructi.FCell = {};
    end
    %  if name field is not present, it defaults to 'x' + ivar
    if ~isfield(modelStructi, 'name')
        modelStructi.name = ['x',num2str(ivar)];
    end
    %  if nallXterm field is not present, it defaults length(XCell)
    if ~isfield(modelStructi, '')
        modelStructi.nallXterm = length(modelStructi.XCell);
    end
    %  if nallFterm field is not present, it defaults length(FCell)
    if ~isfield(modelStructi, 'nallFterm')
        modelStructi.nallFterm = length(modelStructi.FCell);
    end
    %  if weight field not present, it defaults to 1
    if ~isfield(modelStructi, 'weight')
        modelStructi.weight = 1;
    end
    %  replace original modelStruct
    modelCell{ivar} = modelStructi;
end

%  stop if errors found at this point

if errwrd
    error('One or more terminal errors encountered.');
end

%  now check the class of each field

% disp('second loop')
for ivar=1:nvar
    modelStructi = modelCell{ivar};
%     disp(modelStructi)
    %  check order field
    order = modelStructi.order;
    if ~isnumeric(order)
        disp(['order field in modelCell{',num2str(ivar), ...
            '} does not contain a numeric object.']);
        errwrd = 1;
    else
        if order ~= round(order)
            disp(['Field order in modelCell{',num2str(ivar), ...
                '} is not an integer.']);
            errwrd = 1;
        end
    end
    %  check XCell
    XCell = modelStructi.XCell;
    if ~iscell(XCell)
        disp(['XCell field in modelCell{',num2str(ivar), ...
            '} does not contain a cell object.']);
        errwrd = 1;
    end
    %  check FCell
    FCell = modelStructi.FCell;
    if ~iscell(FCell)
        disp(['XCell field in modelCell{',num2str(ivar), ...
            '} does not contain a cell object.']);
        errwrd = 1;
    end
    %  check nXterm
    nXterm = modelStructi.nallXterm;
    if ~isnumeric(nXterm)
        disp(['nXterm field in modelCell{',num2str(ivar), ...
            '} does not contain a numeric object.']);
        errwrd = 1;
    else
        if nXterm ~= round(nXterm)
            disp(['Field nXterm in modelCell{',num2str(ivar), ...
                '} is not an integer.']);
            errwrd = 1;
        end
    end
    %  check nFterm
    nFterm = modelStructi.nallFterm;
    if ~isnumeric(nFterm)
        disp(['nFterm field in modelCell{',num2str(ivar), ...
            '} does not contain a numeric object.']);
        errwrd = 1;
    else
        if nFterm ~= round(nFterm)
            disp(['Field nFterm in modelCell{',num2str(ivar), ...
                '} is not an integer.']);
            errwrd = 1;
        end
    end
    %  check weight
    weight = modelStructi.weight;
    if ~isnumeric(weight)
        disp(['weight field in modelCell{',num2str(ivar), ...
            '} does not contain a numeric object.']);
        errwrd = 1;
    else
        if weight < 0
            disp(['Field weight in modelCell{',num2str(ivar), ...
                '} is negative.']);
            errwrd = 1;
        end
    end
end

%  stop if errors found at this point

if errwrd
    error('One or more terminal errors encountered.');
end

%  Now check the contents of field Xterm

for ivar=1:nvar
    modelStructi = modelCell{ivar};
%     disp(modelStructi)
    if modelStructi.nallXterm ~= length(modelStructi.XCell)
        disp(['Length of Xterm in modelCell{',num2str(ivar), ...
            '} is not equal to field nXterm.']);
        errwrd = 1;
    else
        if ~isempty(modelStructi.XCell)
            for iterm=1:modelStructi.nallXterm
                modelStructij = modelStructi.XCell{iterm};
                if ~isstruct(modelStructij)
                    disp(['XCell{',num2str(iterm), ...
                        '} in modelCell{',num2str(ivar), ...
                        '} is not a struct object.']);
                    errwrd = 1;
                else
                    %  check ncoef field
                    if ~isfield(modelStructij,'ncoef')
                        disp(['Struct object in XCell{', ...
                            num2str(iterm), ...
                            '} in modelCell{', ...
                            num2str(ivar), ...
                            '} does not have a ncoef field.']);
                        errwrd = 1;
                    else
                        ncoefi     = modelStructij.ncoef;
                        if ncoefi < 1 || ncoefi > ncoef
                            disp(['Struct field ncoef in XCell{', ...
                                num2str(iterm), ...
                                '} in modelCell{', ...
                                num2str(ivar), ...
                                '} is not in coefCell.']);
                            errwrd = 1;
                        end
                    end
                    %  check variable field
                    if ~isfield(modelStructij,'variable')
                        disp(['Struct object in XCell{', ...
                            num2str(iterm), ...
                            '} in modelCell{', ...
                            num2str(ivar), ...
                            '} does not have an variable field.']);
                        errwrd = 1;
                    else
                        if ~isnumeric(modelStructij.variable) || ...
                            round(modelStructij.variable) ~= ...
                            modelStructij.variable
                            disp(['Struct object in XCell{', ...
                                num2str(iterm), ...
                                '} in modelCell{', num2str(ivar), ...
                                '} does not have an integer ', ...
                                'variable field.']);
                            errwrd = 1;
                        else
                            if modelStructij.variable > nvar
                                warning(['Struct object in XCell{', ...
                                    num2str(iterm), ...
                                    '} in modelCell{', num2str(ivar), ...
                                    '} has variable field > NVAR.']);
                                errwrd = 1;
                            end
                        end
                    end
                    %  check derivative field
                    if ~isfield(modelStructij,'derivative')
                        disp(['Struct object in XCell{', ...
                            num2str(iterm), ...
                            '} in modelCell{', ...
                            num2str(ivar), ...
                            '} does not have an derivative field.']);
                        errwrd = 1;
                    else
                        if ~isnumeric(modelStructij.derivative) || ...
                            round(modelStructij.derivative) ~= ...
                            modelStructij.derivative
                            disp(['Struct object in XCell{', ...
                                num2str(iterm), ...
                                '} in modelCell{', num2str(ivar), ...
                                '} does not have an integer ', ...
                                'derivative field.']);
                            errwrd = 1;
                        else
                            if modelStructij.derivative >= ...
                                modelCell{modelStructij.variable}.order
                                disp(['Struct object in XCell{', ...
                                    num2str(iterm), ...
                                    '} in modelCell{', num2str(ivar), ...
                                    '} has derivative field >= order.']);
                                errwrd = 1;
                            end
                        end
                    end
                end
            end
        end
    end
end

%  Now check the contents of field FCell

for ivar=1:nvar
    modelStructi = modelCell{ivar};
    if modelStructi.nallFterm ~= length(modelStructi.FCell)
        disp(['Length of Fterm in modelCell{',num2str(ivar), ...
            '} is not equal to field nFterm.']);
        errwrd = 1;
    else
        if ~isempty(modelStructi.FCell)
            for iterm=1:modelStructi.nallFterm
                modelStructij = FCell{iterm};
                if ~isstruct(modelStructij)
                    disp(['FCell{',num2str(iterm), ...
                        '} in modelCell{',num2str(ivar), ...
                        '} is not a struct object.']);
                    errwrd = 1;
                else
                    %  check ncoef field
                    if ~isfield(modelStructij,'ncoef')
                        disp(['Struct object in FCell{', ...
                            num2str(iterm), ...
                            '} in modelCell{', ...
                            num2str(ivar), ...
                            '} does not have a ncoef field.']);
                        errwrd = 1;
                    else
                        ncoefi = modelStructij.ncoef;
                        if ncoefi < 1 || ncoefi > ncoef
                            disp(['Struct field ncoef in FCell{', ...
                                num2str(iterm), ...
                                '} in modelCell{', ...
                                num2str(ivar), ...
                                '} is not in coefCell.']);
                            errwrd = 1;
                        end
                    end
                    %  check Ufd field
                    if ~isfield(modelStructij,'Ufd')
                        disp(['Struct object in FCell{', ...
                            num2str(iterm), ...
                            '} in modelCell{', ...
                            num2str(ivar), ...
                            '} does not have a Ufd field.']);
                        errwrd = 1;
                    else
                        Ufd = modelStructij.Ufd;
                        if ~isa_fd(Ufd)
                            disp(['Struct field Ufd in FCell{', ...
                                num2str(iterm), ...
                                '} in modelCell{', ...
                                num2str(ivar), ...
                                '} is not a fd object.']);
                            errwrd = 1;
                        end
                    end
                end
            end
        end
    end
end

%  stop if errors found at this point

if errwrd
    error('One or more terminal errors encountered.');
end

%  check that struct objects contain a field named 'factor', and, if
%  not, or if empty, replace by factor = 1.

for ivar=1:nvar
    modelStructi = modelCell{ivar};
    for iterm=1:modelStructi.nallXterm
        modelStructij = modelStructi.XCell{iterm};
        if ~isfield(modelStructij, 'factor') || ...
            isempty(modelStructij.factor)
            modelStructij.factor = 1;
            modelStructi.XCell{iterm} = modelStructij;
        end
    end
    for iterm=1:modelStructi.nallFterm
        modelStructij = modelStructi.FCell{iterm};
        if ~isfield(modelStructij, 'factor') || ...
            isempty(modelStructij.factor)
            modelStructij.factor = 1;
            modelStructi.FCell{iterm} = modelStructij;
        end
    end
    modelCell{ivar} = modelStructi;
end

modelCellnew = modelCell;
