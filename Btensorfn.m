function BtensorCell = Btensorfn(XbasisCell, modelCell, coefCell)
%  Set up BTENSORCELL defining homogeneous part of LX;  
%  This is a cell array of length NVAR
%    Each cell BTENSORCELL{i_1} contains a square cell array of order 
%      NALLXTERM+1.
%      Each cell BTENSORCELL{i_1}{j_1,j_2} contains the tensor product
%      of the coefficient basis functions for basis i in term j_1, 
%         the derivative order j of the X basis functions in term j_1,
%         the derivative order j of the X basis functions in term j_2,
%         the coefficient basis functions for basis i in term j_2,
%         j_1, j_2 = 1,...,NALLXTERM+1 where
%         j_1, j_2 = NALLXTERM+1 is for the mth derivative of the X basis
%         functions.

%  Last modified 22 January 2017

rng     = getbasisrange(XbasisCell{1});
Wbasism = create_constant_basis(rng);
Wfdm    = fd(1,Wbasism);
WfdParm = fdPar(Wfdm, 0, 0, 0);

%  set up the structure of BtensorCell

nvar = length(modelCell);
BtensorCell = cell(nvar,1);
for ivar=1:nvar
    modelStructi      = modelCell{ivar};
    BtensorCell{ivar} = cell(modelStructi.nallXterm+1);
    nderiv = modelStructi.nallXterm + 1;
    %  Compute four-way tensors for cells BtensorCell{ivar}{iw,ix}
    %  unless one of the coefStruct.fun objects is a struct, in
    %  which case the cell is empty.
    for iw = 1:nderiv        
%         disp(['iw = ',num2str(iw)])
        if iw < nderiv
            modelStructiw = modelStructi.XCell{iw};
            ncoefw   = modelStructiw.ncoef;
            coefStrw = coefCell{ncoefw};
            WfdParw  = coefStrw.fun;
            Xbasisw  = XbasisCell{modelStructiw.variable};
            jw       = modelStructiw.derivative;
        else
            WfdParw = WfdParm;
            Xbasisw = XbasisCell{ivar};
            jw = modelStructi.order;
        end
        if ~isstruct(WfdParw)
            if isa_basis(WfdParw)
                Wbasisw = WfdParw;
            elseif isa_fd(WfdParw)
                Wbasisw = getbasis(WfdParw);
            else
                Wbasisw = getbasis(getfd(WfdParw));
            end
            Wtypew  = getbasistype(Wbasisw);
            Xtypew  = getbasistype(Xbasisw);
            for ix = 1:nderiv
%                 disp(['ix = ',num2str(ix)])
                if ix < nderiv
                    modelStructix = modelStructi.XCell{ix};
                    ncoefx   = modelStructix.ncoef;
                    coefStrx = coefCell{ncoefx};
                    WfdParx  = coefStrx.fun;
                    Xbasisx  = XbasisCell{modelStructix.variable};
                    jx       = modelStructix.derivative;
                else
                    WfdParx = WfdParm;
                    Xbasisx = XbasisCell{ivar};
                    jx      = modelStructi.order;
                end
                if ~isstruct(WfdParx)
                    if isa_basis(WfdParx)
                        Wbasisx = WfdParx;
                    elseif isa_fd(WfdParx)
                        Wbasisx = getbasis(WfdParx);
                    else
                        Wbasisx = getbasis(getfd(WfdParx));
                    end
                    Wtypex  = getbasistype(Wbasisx);
                    Xtypex  = getbasistype(Xbasisx);
                    if      strcmp(Wtypew, 'const'  ) && ...
                            strcmp(Wtypex, 'const'  ) && ...
                            strcmp(Xtypew, 'bspline') && ...
                            strcmp(Xtypex, 'bspline')
                        XWXWmatij = inprod(Xbasisw, Xbasisx, jw, jx);
                        XWXWmatij = reshape(XWXWmatij', ...
                            getnbasis(Xbasisw)*getnbasis(Xbasisx),1);
                    else
%                         disp([iw,ix,jw,jx])
                        XWXWmatij = inprod_TPbasis( ...
                            Xbasisw, Wbasisw, Xbasisx, Wbasisx, ...
                            jw, 0, jx, 0);
                    end
                    BtensorCell{ivar}{iw,ix} = XWXWmatij;
                end
            end
        end
    end
end
