function [modelCell,coefCell] = create_modelCell(coef_ODE,F,basisobj,estimate)

%Number of Equations
neqns     = size(estimate,1);

%Model Cell
modelCell = cell(neqns,1);

%# find empty cells
emptyCells = cellfun(@isempty,coef_ODE);

%Coef Cell
coefCell  = cell(sum(sum(~emptyCells)),1);

q        = 1;

for k = 1:neqns
    
    %% Set up modelCell
    if(~isempty(F{k}))
        nbeta = size(coef_ODE,2) - size(F,2);
    else
        nbeta = size(coef_ODE,2);
    end
    nzbeta = nbeta - sum(emptyCells(k,:));
    
    % Set up homogenous part of PDE
    modelStruct.order = nbeta;
    
    % Set up the derivative informaation
    derivative_info = 0';
    for j = 1:(modelStruct.order-1)
        derivative_info = [derivative_info,j];
    end
    
    % Set up terms of homogenous part of PDE
    modelStruct.XCell = cell(nzbeta,1);
    
    %Set up coefCell
    rng       = getbasisrange(basisobj{1});
    Wbasisobj = create_constant_basis([rng(1),rng(2)]);
    q1        = 1;
    q2        = 1;
    m2        = 0;
    m4        = nbeta;
    
    for i = 1:nbeta
        if(emptyCells(k,i)==0)
            %ModelCell
            XtermStruct.variable   = k;
            XtermStruct.derivative = derivative_info(i)';
            XtermStruct.ncoef      = q;
            XtermStruct.factor     = 1;
            modelStruct.XCell{q1}  = XtermStruct;
            
            %CoefCell
            if(~isa_fd(coef_ODE{k,i}) & ~isstruct(coef_ODE{k,i}))
                Wfd                    = fd(coef_ODE{k,i}, Wbasisobj);
                coefStr1.fun           = fdPar(Wfd, 0, 0, estimate(k,i));
                coefStr1.parvec        = coef_ODE{k,i};
                coefStr1.estimate      = estimate(k,i);
                coefStr1.index         = q2;
                q2                     = q2+1;
            elseif(isstruct(coef_ODE{k,i}))
                coefStr1.fun           = coef_ODE{k,i};
                if(~isfield(coefStr1.fun, 'more'))
                    coefStr1.fun.more      = [];
                end
                coefStr1.parvec        = coef_ODE{k,i}.theta0;
                m1                     = m2 + 1;
                m2                     = m2 + length(coef_ODE{k,i}.theta0);
                coefStr1.estimate      = any(estimate(k,m1:m2));
                coefStr1.index         = find(estimate(k,m1:m2));
                q2                     = q2 + length(coefStr1.parvec);
            else
                Wfd                    = coef_ODE{k,i};
                coefStr1.fun           = fdPar(Wfd, 0, 0, 1);
                coefStr1.parvec        = getcoef(coef_ODE{k,i});
                m1                     = m2 + 1;
                m2                     = m2 + length(coefStr1.parvec);
                coefStr1.estimate      = any(estimate(k,m1:m2));
                coefStr1.index         = find(estimate(k,m1:m2));
                q2                     = q2 + length(coefStr1.parvec);         
            end
            coefStr1.coeftype      = 'beta';
            coefCell{q}            = coefStr1;
            q                      = q+1;
            q1                     = q1+1;
        end
    end
    if(~isempty(F{k}))
        %Set up forcing function
        modelStruct.FCell = {};
        for j = 1:size(F,2)
            if(~isempty(F{k,j}))
                %CoefCell
                if(~isa_fd(coef_ODE{k,nbeta+j}) & ~isstruct(coef_ODE{k,nbeta+j}))
                    Wfd                    = fd(coef_ODE{k,nbeta+j}, Wbasisobj);
                    coefStr1.fun           = fdPar(Wfd, 0, 0, estimate(k,nbeta+j));
                    coefStr1.parvec        = coef_ODE{k,nbeta+j};
                    coefStr1.estimate      = estimate(k,nbeta+j);
                    coefStr1.index         = q2;
                    q2                     = q2+1;
                elseif(isstruct(coef_ODE{k,nbeta+j}))
                    coefStr1.fun           = coef_ODE{k,nbeta+j};
                    if(~isfield(coefStr1.fun, 'more'))
                        coefStr1.fun.more      = [];
                    end
                    coefStr1.parvec        = coef_ODE{k,nbeta+j}.theta0;
                    m3                     = m4 + 1;
                    m4                     = m4 + length(coef_ODE{k,nbeta+j}.theta0);
                    coefStr1.estimate      = any(estimate(k,nbeta-1+m3:m4));
                    coefStr1.index         = q2 + find(estimate(k,nbeta-1+m1:m2));
                    q2                     = q2 + length(coefStr1.parvec);
                else
                    Wfd                    = coef_ODE{k,nbeta+j};
                    coefStr1.fun           = fdPar(Wfd, 0, 0, 1);
                    coefStr1.parvec        = getcoef(coef_ODE{k,nbeta+j});
                    m3                     = m4 + 1;
                    m4                     = m4 + length(coefStr1.parvec);
                    coefStr1.estimate      = any(estimate(k,m3:m4));
                    coefStr1.index         = q2 + find(estimate(k,m3:m4))';
                    q2                     = q2 + length(coefStr1.parvec);
                end
                coefStr1.coeftype     = 'alpha';
                coefCell{q}           = coefStr1;
                q                     = q+1;
                
                %FCell
                if(isa_fd(F{k,j}))
                    FStruct.ncoef         = q-1;
                    FStruct.Ufd           = F{k,j};
                    FStruct.factor        = 1;
                    modelStruct.FCell{j,1}  = FStruct;
                elseif(isnumeric(F{k,j}))
                    FStruct.ncoef         = q-1;
                    FStruct.AfdPar        = fdPar(fd(F{k,j}, Wbasisobj), 0, 0, 1);
                    FStruct.Ufd           = fd(1, Wbasisobj);
                    FStruct.factor        = 1;
                    modelStruct.FCell{j,1}  = FStruct;
                end
            end
        end
    end
    
    %Complete model structure
    modelCell{k} = modelStruct;
    clear modelStruct
    
end

%Checks

[check, ntheta] = coefcheck(coefCell);

modelCell = modelcheck(modelCell, coefCell);

end