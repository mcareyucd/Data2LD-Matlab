function AtensorCell = Atensorfn(modelCell, coefCell)
%  Set up ATENSORCELL of length NVAR defining products of forcing terms.
%  Each celll AtensorCell{ivar} is a cell of size nforcei by nforcei.
%
%  This version was used prior to introducing coefCell as an argument.
%
%  Last modified 22 January 2017

nvar = length(modelCell);
AtensorCell = cell(nvar,1);
for ivar=1:nvar
    modelStructi = modelCell{ivar};
    if modelStructi.nallFterm > 0
        nforce = modelStructi.nallFterm;
        AtensorCell{ivar} = cell(nforce,nforce);
        %  Compute four-way tensors for cells BtensorCell{ivar}{iw,ix}
        %  unless one of the coefStruct.fun objects is a struct, in
        %  which case the cell is empty.
        for jv = 1:nforce
            modelStructiv = modelStructi.FCell{jv};
            ncoefv     = modelStructiv.ncoef;
            coefStrv   = coefCell{ncoefv};
            AfdParv    = coefStrv.fun;
            if ~isstruct(AfdParv)
                Ufdv    = modelStructiv.Ufd;
                if isa_basis(AfdParv)
                    Abasisv = AfdParv;
                elseif isa_fd(AfdParv)
                    Abasisv = getbasis(AfdParv);
                else
                    Abasisv = getbasis(getfd(AfdParv));
                end
                Ubasisv = getbasis(Ufdv);
                Atypev  = getbasistype(Abasisv);
                Utypev  = getbasistype(Ubasisv);
                for jx = 1:nforce
                    modelStructix = modelStructi.FCell{jx};
                    ncoefx     = modelStructix.ncoef;
                    coefStrx   = coefCell{ncoefx};
                    AfdParx    = coefStrx.fun;
                    if ~isstruct(AfdParx)
                        Ufdx    = modelStructix.Ufd;
                        if isa_basis(AfdParx)
                            Abasisx = AfdParx;
                        elseif isa_fd(AfdParx)
                            Abasisx = getbasis(AfdParx);
                        else
                            Abasisx = getbasis(getfd(AfdParx));
                        end
                        Ubasisx = getbasis(Ufdx);
                        Atypex  = getbasistype(Abasisx);
                        if strcmp(Atypev, 'const'  ) && ...
                                strcmp(Atypex, 'const'  )
                            XWXWmatij = inprod(Ubasisv, Ubasisx);
%                             tmp_fun = @(x)(eval_basis(x,...
%                             Ubasisv).*eval_basis(x,Ubasisx));
%                             rng = getbasisrange(Ubasisv);
%                             XWXWmatij = integral(tmp_fun,rng(1),...
%                                         rng(2),'ArrayValued',true);
                            XWXWmatij = reshape(XWXWmatij', ...
                                getnbasis(Ubasisv)*getnbasis(Ubasisx),1);
                        else
                            XWXWmatij = inprod_TPbasis( ...
                                Ubasisv, Abasisv, ...
                                Ubasisx, Abasisx, ...
                                0, 0, 0, 0);
                        end
                        AtensorCell{ivar}{jx,jv} = XWXWmatij;
                    end
                end
            end
        end
    else
        AtensorCell{ivar} = {};
    end
end

