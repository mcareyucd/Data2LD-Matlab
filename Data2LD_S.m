function [Smat, DSarray] = Data2LD_S(XbasisCell, modelCell, coefCell, ...
                                      rhoVec, ntheta, BAtensorCell, ...
                                      nrep, nforce)
%  Data2LD ... stands for "Data to Linear Dynamics"
%  Data2LD_S computes the penalty matrix S associated with the forcing
%  portion of a linear differential operator L, as well as its partial
%  derivatives with  respect to the parameter vector.
%  For a single variable whose approximated in terms of an exansion in
%  terms of a vector \phi of basis functions, S is
%                 S = \int [L \phi(t)] U' dt.
%  S has dimensions K and NREP, where K is the number of basis
%  functions in the expansion of the variable, NFORCE is the number of
%  forcing functions, and NREP is the number of replications.  The
%  forcing functions are assumed to vary from one replication to another.
%  This version loads BAtensorCell as an argument.
%
%  If multiple variables are involved, then S is a composite matrix
%  constructed from inner products and cross-products of the basis 
%  function vectors associate with each variable.  It's dimension will be
%  \sum K_i by NFORCE*NREP.
%
%  This version approximates the integrals in the penalty terms by using 
%  inprod_basis to compute the cross-product matrices for the  
%  \beta-coefficient basis functions and the corresponding derivative of 
%  the x-basis functions,and the cross-product matrices for the 
%  \alpha-coefficients and the corresponding U functions.  
%  These are computed upon the first call to Data2LD, and then retained 
%  for subsequent calls by using the persistent command.  See lines about 
%  560 to 590 for this code.
%
%  This version disassociates coefficient functions from equation 
%  definitions to allow some coefficients to be used repeatedly and for
%  both homogeneous and forcing terms.  It requires an extra argument
%  COEFCELL that contains the coefficients and the position of their
%  coefficient vectors in vector THETA.
%
%  Arguments:
%
%  XBASISCELL ... A functional data object or a BASIS object.  If so, the 
%               smoothing parameter LAMBDA is set to 0.
%  MODELCELL  ... A cell aray of length NVAR. Each cell contains a 
%                struct object with members:              
%                XCell ... cell array of length number of homogeneous terms
%                          Each cell contains a struct object with members:
%                          WfdPar     ... fdPar object for the coefficient
%                          variable   ... the index of the variable
%                          derivative ... the order of its derivative
%                          ncoef      ... if coefficient estimated, its 
%                                         location in the composite vector 
%                          factor     ... a scalar multiplier (def. 1)
%                FCell ... cell arrau of length number of forcing terms
%                          Each cell contains a struct object with members:
%                          AfdPar ... an fdPar object for the coefficient
%                          Ufd    ... an fd object for the forcing function
%                          ncoef  ... if coefficient estimated, its 
%                                     location in the composite vector 
%                          factor ... a scalar multiplier (def. 1)
%                order     ... the highest order of derivative
%                name      ... a  tag for the variable
%                nallXterm ... the number of homogeneous terms
%                nallFterm ... the number of forcing functions
%  COEFCELL  ... A cell array of length NCOEF containing struct objects
%                with fields:
%               parvec   ... a vector of parameters
%               estimate ... 0, held fixed, otherwise, estimated 
%               coeftype ... homogeneous or forcing
%               fun      ... functional basis, fd, or fdPar object, 
%                            or a struct object for a general function 
%                            with fields:
%                 fd      ... function handle for evaluating function
%                 Dfd     ... function handle for evaluating 
%                             partial derivative with respect to parameter
%                 more    ... object providing additional information for 
%                             evaluating coefficient function
%  RHOVEC    ... A vector of length NVAR containing values in [0,1].  
%                The data sums of squares are weighted by P and 
%                the roughness penalty by 1-P.
%  BAtensorCell ... A cell array of four-way tensors required for Smat
%  NREP       ... The number of replications of the system.  

%  Last modified 23 December 2016

%  ------------------------------------------------------------------------
%                         Set up analysis
%  ------------------------------------------------------------------------

if nargin < 4
    error('Less than four parameters supplied.');
end

%  compute number of variables

nvar = length(modelCell);

if nargin < 4,  rhoVec = 0.5*ones(nvar,1);  end
    
%  Set up a vector NCOEFVEC containing number of coefficients used
%  for the expansion of each variable

ncoefvec = zeros(nvar,1);
for ivar=1:nvar
    ncoefvec(ivar) = getnbasis(XbasisCell{ivar});
end
ncoefcum = cumsum([0;ncoefvec]);

%  get the width of the time domain

Xrange = getbasisrange(XbasisCell{1});
T      = Xrange(2) - Xrange(1);

%--------------------------------------------------------------------------
%                 Compute the penalty vector S(theta)
%--------------------------------------------------------------------------

ncoefsum = sum(ncoefvec);
if sum(nforce) == 0
    Smat = [];
    DSarray = [];
    return
end

Smat = zeros(ncoefsum,nrep);
m2 = 0;
for ivar=1:nvar
    modelStructi = modelCell{ivar};
    weighti = modelStructi.weight;
    m1  = m2 + 1;
    m2  = m2 + ncoefvec(ivar);
    ind = m1:m2;
    Xbasisi  = XbasisCell{ivar};
    order    = modelStructi.order;
    if nforce(ivar) > 0
        nXbasisi    = ncoefvec(ivar);
        nXtermi     = modelStructi.nallXterm;
        nFtermi     = modelStructi.nallFterm;
        for jforce = 1:nFtermi
            modelStructij = modelStructi.FCell{jforce};
            ncoefj   = modelStructij.ncoef;
            coefStrj = coefCell{ncoefj};
            Avecj    = coefStrj.parvec;
            factorj  = modelStructij.factor;
            Ufdj     = modelStructij.Ufd;
            Ubasisj  = getbasis(Ufdj);
            Ucoefj   = getcoef(Ufdj);
            nUbasisj = getnbasis(Ubasisj);
            nAbasisj = length(Avecj);
            funtypej = isstruct(coefStrj.fun);
            %  Crossproducts of homogeneous terms with forcing terms
            for iw=1:nXtermi
                modelStructiw = modelStructi.XCell{iw};
                ivw       = modelStructiw.variable;
                indw     = ncoefcum(ivw)+1:ncoefcum(ivw+1);
                nXbasisw = ncoefvec(ivw);
                Xbasisw  = XbasisCell{ivw};
                ncoefw   = modelStructiw.ncoef;
                coefStrw = coefCell{ncoefw};
                Bvecw    = coefStrw.parvec;
                factorw  = modelStructiw.factor;
                nWbasisw = length(Bvecw);
                derivw   = modelStructiw.derivative;
                funtypew = isstruct(coefStrw.fun);
                if funtypej || funtypew
                    [Smatjw,iter] = ...
                        inprod_basis_Data2LD(Xbasisw,  Ufdj,     ...
                                          coefStrw, coefStrj, ...
                                          derivw,   0);
                else
                    BAtenswj = BAtensorCell{ivar}{iw,jforce};
                    Smatjw = inprodijw(nXbasisw, nWbasisw, ...
                                      nUbasisj, nAbasisj, nrep, ...
                                      Bvecw, Avecj, Ucoefj, BAtenswj);
                end
                Smatjw = factorj*factorw*rhoVec(ivar).*Smatjw/T;
                Smat(indw,:) = Smat(indw,:) + weighti*Smatjw;
            end
            %  Crossproducts of D^m with forcing terms
            if funtypej
                [Smatji, iter] = inprod_basis_Data2LD(Xbasisi,  Ufdj, ...
                                                  1, coefStrj, ...
                                                  order,    0);
            else
                BAtenswj = BAtensorCell{ivar}{nXtermi+1,jforce};
                Smatji = inprodijw(nXbasisi, 1, ...
                                  nUbasisj, nAbasisj, nrep, ...
                                  1, Avecj, Ucoefj, BAtenswj);
            end
            Smatji = factorj.*rhoVec(ivar).*Smatji/T;
            Smat(ind,:) = Smat(ind,:) - weighti*Smatji;
        end
    end
end

if nargout > 1

DSarray = zeros(ncoefsum,nrep,ntheta);
m2 = 0;
for ivar=1:nvar
    m1   = m2 + 1;
    m2   = m2 + ncoefvec(ivar);
    indi = m1:m2;
    if nforce(ivar) > 0
        modelStructi = modelCell{ivar};
        weighti  = modelStructi.weight;
        nXbasisi = ncoefvec(ivar);
        nXtermi  = modelStructi.nallXterm;
        nFtermi  = modelStructi.nallFterm;
        Xbasisi  = XbasisCell{ivar};
        order    = modelStructi.order;
        %  partial derivatives of product of homogeneous terms
        %  and forcing terms with respect to homogeneous coefficients
        %  loop through all active forcing terms
        for jforce = 1:nFtermi
            modelStructij = modelStructi.FCell{jforce};
            ncoefj   = modelStructij.ncoef;
            coefStrj = coefCell{ncoefj};
            Avecj    = coefStrj.parvec;
            Aestimj  = coefStrj.estimate;
            factorj  = modelStructij.factor;
            nAbasisj = length(Avecj);
            Ufdj     = modelStructij.Ufd;
            Ucoefj   = getcoef(Ufdj);
            Ubasisj  = getbasis(Ufdj);
            nUbasisj = getnbasis(Ubasisj);
            funtypej = isstruct(coefStrj.fun);
            %  crossproducts of homogeneous terms with forcing terms
            for iw=1:nXtermi
                modelStructiw = modelStructi.XCell{iw};
                ncoefw   = modelStructiw.ncoef;
                coefStrw = coefCell{ncoefw};
                Westimw  = coefStrw.estimate;
                nderivw  = modelStructiw.derivative;
                funtypew = isstruct(coefStrw.fun);
                factorw  = modelStructiw.factor;
                if Westimw
                    ivw      = modelStructiw.variable;
                    indw     = ncoefcum(ivw)+1:ncoefcum(ivw+1);
                    indthw   = coefStrw.index;
                    nWbasisw = length(indthw);
                    nXbasisw = ncoefvec(ivw);
                    if funtypej || funtypew
                        [DBSarrayjw, iter] = ...
                            inprod_Dbasis_Data2LD(Xbasisw,  Ufdj, ...
                                                  coefStrw, coefStrj, ...
                                                  nderivw,  0);
                    else
                        BAtenswj   = BAtensorCell{ivar}{iw,jforce};
                        DBSarrayjw = inprodwj(nXbasisw, nWbasisw, ...
                                              nUbasisj, nAbasisj, ...
                                              nrep, Avecj, Ucoefj, BAtenswj);
                    end
                    DBSarrayjw = factorj*factorw*rhoVec(ivar)*DBSarrayjw/T;
                    DSarray(indw,:,indthw) = DSarray(indw,:,indthw) + ...
                                             weighti*DBSarrayjw;
                end
            end
            %  partial derivatives wrt forcing term coefficients for
            %  those forcing terms requiring estimation of their
            %  coefficients
            %  loop through all forcing terms with coefficients to be
            %  estimated
            if Aestimj
                indtha  = coefStrj.index;
                %  partial derivatives of products of homogeneous terms
                %  with forcing terms wrt forcing coefficients
                for iw=1:nXtermi
                    modelStructiw = modelStructi.XCell{iw};
                    ivw      = modelStructiw.variable;
                    indw     = ncoefcum(ivw)+1:ncoefcum(ivw+1);
                    nXbasisw = ncoefvec(ivw);
                    ncoefw   = modelStructiw.ncoef;
                    coefStrw = coefCell{ncoefw};
                    Bvecw    = coefStrw.parvec;
                    funtypew = isstruct(coefStrw.fun);
                    factorw  = modelStructiw.factor;
                    nderivw  = modelStructiw.derivative;
                    nWbasisw = length(Bvecw);
                    if funtypew || funtypej
                        [DASarrayjw, iter] = ...
                            inprod_Dbasis_Data2LD(Ufdj,     Xbasisw,  ...
                                               coefStrj, coefStrw, ...
                                               0,        nderivw);
                        DASarrayjw = permute(DASarrayjw,[2,1,3]);
                    else
                        BAtenswj = BAtensorCell{ivar}{iw,jforce};
                        DASarrayjw  = inprodjw(nXbasisw, nWbasisw, ...
                                            nUbasisj, nAbasisj, ...
                                            nrep, Bvecw, Ucoefj, BAtenswj);
                    end
                    DASarrayjw = factorj*factorw*rhoVec(ivar)*DASarrayjw/T;
                    DASarrayjw = weighti*DASarrayjw;
                    DSarray(indw,:,indtha) = DSarray(indw,:,indtha) + ...
                                             DASarrayjw;
                end
                %  partial derivatives of cross-products of D^m
                %  with forcing terms wrt forcing coefficients
                if funtypej
                    [DASarrayji, iter] = ...
                              inprod_Dbasis_Data2LD(Ufdj,  Xbasisi,     ...
                                                 coefStrj, coefStri, ...
                                                 0,        order);
                    DASarrayji = permute(DASarrayji,[2,1,3]);
               else
                    BAtensij = BAtensorCell{ivar}{nXtermi+1,jforce};
                    DASarrayji  = inprodjw(nXbasisi, 1, ...
                                        nUbasisj, nAbasisj, ...
                                        nrep, 1, Ucoefj, BAtensij);
                end
                DASarrayji = factorj*rhoVec(ivar)*DASarrayji/T;
                DSarray(indi,:,indtha) = DSarray(indi,:,indtha) - ...
                                         weighti*DASarrayji;
            end
        end
    end
end

end

%  ------------------------------------------------------------------------

function DASarray = inprodjw(nXbasisw, nWbasisw, nUbasisj, nAbasisj, ...
                            nrep, Bvecw, Ucoefj, BAtenswj)
ncum = cumprod([nAbasisj, nUbasisj, nWbasisw, nXbasisw]);
DASarray = zeros(nXbasisw,nrep,nAbasisj);
for irep=1:nrep
    for i=1:nXbasisw
        for j=1:nWbasisw
            for k=1:nUbasisj
                for l=1:nAbasisj
                    ijkl = (i-1)*ncum(3) + ...
                           (j-1)*ncum(2) + ...
                           (k-1)*ncum(1) + l;
                    DASarray(i,irep,l) = DASarray(i,irep,l) + ...
                        Bvecw(j)*BAtenswj(ijkl)* ...
                        Ucoefj(k,irep);
                end
            end
        end
    end
end

%  ------------------------------------------------------------------------

function DBSarray = inprodwj(nXbasisw, nWbasisw, nUbasisj, nAbasisj, ...
                            nrep, Avecj, Ucoefj, BAtenswj)
ncum = cumprod([nAbasisj, nUbasisj, nWbasisw, nXbasisw]);
DBSarray = zeros(nXbasisw,nrep,nWbasisw);
% BAtensor = zeros(nXbasisw, nWbasisw, nUbasisj, nAbasisj);
for irep=1:nrep
    for i=1:nXbasisw
        for j=1:nWbasisw
            for k=1:nUbasisj
                for l=1:nAbasisj
                    ijkl = (i-1)*ncum(3) + ...
                        (j-1)*ncum(2) + ...
                        (k-1)*ncum(1) + l;
                    DBSarray(i,irep,j) = DBSarray(i,irep,j) + ...
                        Avecj(l)*BAtenswj(ijkl)*Ucoefj(k,irep);
                end
            end
        end
    end
end

%  ------------------------------------------------------------------------

function Smatjw = inprodijw(nXbasisw, nWbasisw, nUbasisj, nAbasisj, ...
                            nrep, Bvecw, Avecj, Ucoefj, BAtenswj)
ncum = cumprod([nAbasisj, nUbasisj, nWbasisw, nXbasisw]);
Smatjw = zeros(nXbasisw,nrep);
for irep=1:nrep
    for i=1:nXbasisw
        for j=1:nWbasisw
            for k=1:nUbasisj
                for l=1:nAbasisj
                    ijkl = (i-1)*ncum(3) + ...
                        (j-1)*ncum(2) + ...
                        (k-1)*ncum(1) + l;
                    Smatjw(i,irep) = Smatjw(i,irep) + ...
                        Bvecw(j)*Avecj(l)*BAtenswj(ijkl)* ...
                        Ucoefj(k,irep);
                end
            end
        end
    end
end

