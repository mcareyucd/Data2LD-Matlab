function BAtensorCell = BAtensorfn(XbasisCell, modelCell, coefCell)
%  Set up BATENSORCELL of size NVAR by NVAR defining products of   
%  derivative terms and forcing terms in LX.
%  Each cell BATENSORCELL{i_1,i_2} contains a cell array of dimensions
%        NFORCE_{j_1} and NDERIVVEC_{l} + 1.
%  Each cell BATENSORCELL{i_1,i_2}{j_1,l_2} contains the inner product 
%  X basis functions for variable i_2,
%  variable weight basis functions for variable i_2, 
%  basis functions for forcing function j_1 and variable i_1
%  forcing weight basis for forcing function j_1 and variable i_1.
%  This version assumes that all forcing functions U have the same basis.

%  Last modified 22 January 2017

rng     = getbasisrange(XbasisCell{1});
Wbasism = create_constant_basis(rng);

%  set up the structure of BAtensorCell

nvar = length(modelCell);
BAtensorCell = cell(nvar,1);
for ivar=1:nvar
    modelStructi = modelCell{ivar};
    nallXterm = modelStructi.nallXterm;
    nallFterm = modelStructi.nallFterm;
    Xbasisi   = XbasisCell{ivar};
    Xtypei    = getbasistype(Xbasisi);
    if nallFterm > 0
        BAtensorCell{ivar} = cell(modelStructi.nallXterm+1, nallFterm);
        order = modelStructi.order;
        %  Compute four-way tensors for cells BtensorCell{ivar}{iw,ix}
        %  unless one of the coefStruct.fun objects is a struct, in
        %  which case the cell is empty.
        for jforce = 1:modelStructi.nallFterm
            modelStructij = modelStructi.FCell{jforce};
            ncoefj    = modelStructij.ncoef;
            coefStrj  = coefCell{ncoefj};
            AfdParj   = coefStrj.fun;
            if ~isstruct(AfdParj)
                if isa_basis(AfdParj)
                    Abasisj = AfdParj;
                elseif isa_fd(AfdParj)
                    Abasisj = getbasis(AfdParj);
                else
                    Abasisj = getbasis(getfd(AfdParj));
                end
                Atypej    = getbasistype(Abasisj);
                Ubasisj   = getbasis(modelStructij.Ufd);
                Utypej    = getbasistype(Ubasisj);
                if      strcmp(Atypej, 'const'  ) && ...
                        strcmp(Utypej, 'bspline') && ...
                        strcmp(Xtypei, 'bspline')
                    XWXWmatij = inprod(Xbasisi, Ubasisj, order, 0);
                    XWXWmatij = reshape(XWXWmatij', ...
                        getnbasis(Xbasisi)*getnbasis(Ubasisj),1);
                else
                    XWXWmatij = inprod_TPbasis(Xbasisi, Wbasism, ...
                                                Ubasisj, Abasisj, ...
                                                order, 0, 0, 0);
                end
                BAtensorCell{ivar}{modelStructi.nallXterm+1,jforce} = ...
                    XWXWmatij;
                for iw = 1:nallXterm
                    modelStructiw = modelStructi.XCell{iw};
                    ncoefw     = modelStructiw.ncoef;
                    coefStrw   = coefCell{ncoefw};
                    WfdParw    = coefStrw.fun;
                    if ~isstruct(WfdParw)
                        if isa_basis(WfdParw)
                            Wbasisw = WfdParw;
                        elseif isa_fd(WfdParw)
                            Wbasisw = getbasis(WfdParw);
                        else
                            Wbasisw = getbasis(getfd(WfdParw));
                        end
                        Wtypew     = getbasistype(Wbasisw);
                        derivative = modelStructiw.derivative;
                        Xbasisw    = XbasisCell{modelStructiw.variable};
                        Xtypew     = getbasistype(Xbasisw);
                        if  strcmp(Wtypew, 'const'  ) && ...
                                strcmp(Atypej, 'const'  ) && ...
                                strcmp(Utypej, 'bspline') && ...
                                strcmp(Xtypew, 'bspline')
                            XWXWmatiw = inprod(Xbasisw, Ubasisj, ...
                                               derivative, 0);
                            XWXWmatiw = reshape(XWXWmatiw', ...
                                getnbasis(Xbasisw)*getnbasis(Ubasisj),1);
                        else
                            XWXWmatiw  = inprod_TPbasis( ...
                                                    Xbasisw, Wbasisw, ...
                                                    Ubasisj, Abasisj, ...
                                                    derivative, 0, 0, 0);
                        end
                        BAtensorCell{ivar}{iw,jforce} = XWXWmatiw;
                    end
                end
            end
        end
    end
end

