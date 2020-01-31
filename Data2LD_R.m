% function [Rmat, DRarray] = Data2LD_R(XbasisCell, modelCell, coefCell, ...
%                                       rhoVec, ntheta, BtensorCell)
function [Rmat, DRarray] = Data2LD_R(XbasisCell, modelCell, coefCell, ...
                                      rhoVec, ntheta, BtensorCell)
%  Data2LD ... stands for "Data to Linear Dynamics"
%  Data2LD_R computes the penalty matrix R associated with the homogeneous
%  portion of a linear differential operator L as well as its partial
%  derivative with respect to parameters defining the homogeneous portion.
%  This version inputs BtensorCell as an argument.
%  For a single variable whose approximated in terms of an exansion in
%  terms of a vector \phi of basis functions, R is
%                 R = \int [L \phi(t)] [L \phi(t)]' dt.
%  R is of order K, the number of basis functions, symmetric, and of rank
%  K - m where m is the order of the largest derivative in the operator.
%  The eigenanalysis of R can be used to construct an alternative basis
%  expansion defined in terms an increasing order of complexity of shape.
%
%  If multiple variables are involved, then R is a composite matrix
%  constructed from inner products and cross-products of the basis 
%  function vectors associate with each variable.  It's order will be
%  \sum K_i.
%
%  This version approximates the integrals in the penalty terms by using 
%  inprod_basis to compute the cross-product matrices for the  
%  \beta-coefficient basis functions and the corresponding derivative of 
%  the x-basis functions,and the cross-product matrices for the 
%  \alpha-coefficients and the corresponding U functions.  
%  These are computed upon the first call to Data2LD4, and then retained 
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
%
%  MODELCELL...  A cell aray of length NVAR. Each cell contains a 
%                struct object with members:              
%                XCell ... cell array of length number of homogeneous terms
%                          Each cell contains a struct object with members:
%                          WfdPar     ... fdPar object for the coefficient
%                          variable   ... the index of the variable
%                          derivative ... the order of its derivative
%                          ncoef      ... if coefficient estimated, its 
%                                         location in the composite vector 
%                          factor     ... a scalar multiplier (def. 1)
%                FCell ... cell array of length number of forcing terms
%                          Each cell contains a struct object with members:
%                          AfdPar ... an fdPar object for the coefficient
%                          Ufd    ... an fd object for the forcing function
%                          ncoef... if coefficient estimated, its location
%                                   in the composite vector 
%                order     ... the highest order of derivative
%                name      ... a  tag for the variable
%                nallXterm ... the number of homogeneous terms
%                nallFterm ... the number of forcing functions
%  COEFCELL  ... A cell array of length NCOEF containing struct objects
%                with fields:
%               parvec   ... a vector of parameters
%               coeftype ... homogeneous or forcing
%               estimate ... 0, held fixed, otherwise, estimated 
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
%
%  BTENSORCELL 

%  Last modified 23 January 2017

%  ------------------------------------------------------------------------
%                         Set up analysis
%  ------------------------------------------------------------------------

if nargin < 2
    error('Less than two parameters supplied.');
end

%  set default arguments

if nargin <  4,  error('BtensorCell not supplied');  end

%  compute number of variables

nvar = length(modelCell);

%  Set up a vector NCOEFVEC containing number of coefficients used
%  for the expansion of each variable

ncoefvec = zeros(nvar,1);
for ivar=1:nvar
    ncoefvec(ivar) = getnbasis(XbasisCell{ivar});
end

%  get the width of the time domain

Xrange = getbasisrange(XbasisCell{1});
T      = Xrange(2) - Xrange(1);

%  ------------------------------------------------------------------------
%                  Compute penalty matrix Rmat(theta)
%  ------------------------------------------------------------------------

ncoefsum = sum(ncoefvec);
Rmat     = zeros(ncoefsum,ncoefsum);
ncoefcum = cumsum([0;ncoefvec]);
m2 = 0;
for ivar=1:nvar
    modelStructi = modelCell{ivar};
    weighti = modelStructi.weight;
    m1  = m2 + 1;
    m2  = m2 + ncoefvec(ivar);
    indi = m1:m2;
    nXbasisi = ncoefvec(ivar);
    Xbasisi  = XbasisCell{ivar};
    nXtermi  = modelStructi.nallXterm;
    order    = modelStructi.order;
    if isempty(BtensorCell)
        order  = modelStructi.order;
        Rmatii = inprod_basis(Xbasisi, Xbasisi, order, order);
    else
        Btensii  = BtensorCell{ivar}{nXtermi+1,nXtermi+1};
        Rmatii   = inprodwx(nXbasisi, 1, nXbasisi, 1, 1, 1, Btensii);
    end
    Rmatii   = rhoVec(ivar)*Rmatii/T;
    Rmat(indi,indi) = Rmat(indi,indi) + weighti*Rmatii;
    for iw=1:nXtermi
        modelStructiw = modelStructi.XCell{iw};
        derivw   = modelStructiw.derivative;
        ivw      = modelStructiw.variable;
        indw     = ncoefcum(ivw)+1:ncoefcum(ivw+1);
        nXbasisw = ncoefvec(ivw);
        Xbasisw  = XbasisCell{ivw};
        ncoefw   = modelStructiw.ncoef;
        coefStrw = coefCell{ncoefw};
        Bvecw    = coefStrw.parvec;
        factorw  = modelStructiw.factor;
        nWbasisw = length(Bvecw);
        funtypew = isstruct(coefStrw.fun);
        for ix=1:nXtermi
            modelStructix = modelStructi.XCell{ix};
            derivx   = modelStructix.derivative;
            ivx      = modelStructix.variable;
            indx     = ncoefcum(ivx)+1:ncoefcum(ivx+1);
            nXbasisx = ncoefvec(ivx);
            Xbasisx  = XbasisCell{ivx};
            ncoefx   = modelStructix.ncoef;
            coefStrx = coefCell{ncoefx};
            Bvecx    = coefStrx.parvec;
            factorx  = modelStructix.factor;
            nWbasisx = length(Bvecx);
            funtypex = isstruct(coefStrx.fun);
            if funtypew || funtypex
                Rmatwx = inprod_basis_Data2LD(Xbasisw,  Xbasisx, ...
                                           coefStrw, coefStrx, ...
                                           derivw,   derivx);
            else
                Btenswx  = BtensorCell{ivar}{iw,ix};
                Rmatwx   = inprodwx(nXbasisw, nWbasisw, ...
                                    nXbasisx, nWbasisx, ...
                                    Bvecw, Bvecx, Btenswx);
            end
            Rmatwx   = factorw*factorx*rhoVec(ivar)*Rmatwx/T;
            Rmat(indw,indx) = Rmat(indw,indx) + weighti*Rmatwx;
        end
        if funtypew
            Rmatiw = inprod_basis_Data2LD(Xbasisi, Xbasisw, ...
                                       [], coefStrw, ...
                                       order, derivw);
        else
            Btensiw = BtensorCell{ivar}{nXtermi+1,iw};
            Rmatiw  = inprodwx(nXbasisi, 1, nXbasisw, nWbasisw, ...
                1, Bvecw, Btensiw);
        end
        Rmatiw  = factorw*rhoVec(ivar)*Rmatiw/T;
        Rmat(indi,indw) = Rmat(indi,indw) - weighti*Rmatiw;
        Rmat(indw,indi) = Rmat(indw,indi) - weighti*Rmatiw';
    end
end

%  ------------------------------------------------------------------------
%  Compute partial derivatives of R with respect to theta 
%  if DRarray is required.
%  ------------------------------------------------------------------------

if nargout > 1
    
%  compute DRarray

DRarray = zeros(ncoefsum,ncoefsum,ntheta);

m2 = 0;
for ivar=1:nvar
    modelStructi = modelCell{ivar};
    weighti = modelStructi.weight;
    m1   = m2 + 1;
    m2   = m2 + ncoefvec(ivar);
    indi = m1:m2;
    nXtermi  = modelStructi.nallXterm;
    nXbasisi = ncoefvec(ivar);
    Xbasisi  = XbasisCell{ivar};
    nderivi  = modelStructi.order;
    for iw=1:nXtermi
        modelStructiw = modelStructi.XCell{iw};
        nderivw  = modelStructiw.derivative;
        ncoefw   = modelStructiw.ncoef;
        coefStrw = coefCell{ncoefw};
        Westimw  = coefStrw.estimate;
        factorw  = modelStructiw.factor;
        if Westimw
            %  define coefficient of estimated variable and
            %  it's derivative index
            ivw      = modelStructiw.variable;
            jvw      = modelStructiw.derivative;
            indw     = ncoefcum(ivw)+1:ncoefcum(ivw+1);
            indthw   = coefStrw.index;
            nXbasisw = ncoefvec(ivw); 
            Xbasisw  = XbasisCell{ivw};
            funtypew = isstruct(coefStrw.fun);
            %  loop through all active variables within equation ivar
            for ix=1:nXtermi
                modelStructix = modelStructi.XCell{ix};
                nderivx  = modelStructix.derivative;
                %  define coefficient of active variable and
                %  it's derivative index
                ivx       = modelStructix.variable;
                jvx       = modelStructix.derivative;
                indx      = ncoefcum(ivx)+1:ncoefcum(ivx+1);
                nXbasisx  = ncoefvec(ivx);
                ncoefx    = modelStructix.ncoef;
                coefStrx  = coefCell{ncoefx};
                Bvecx     = coefStrx.parvec;
                funtypex  = isstruct(coefStrx.fun);
                factorx   = modelStructix.factor;
                nWbasisx  = length(Bvecx);
                %  get the tensor vector for this pair of coefficients
                %  and derivatives
                if funtypew || funtypex
                    DRarraywx = ...
                        inprod_Dbasis_Data2LD(Xbasisw,  Xbasisx,  ...
                                              coefStrw, coefStrx, ...
                                              nderivw,  nderivx);
                else
                    Btenswx   = BtensorCell{ivar}{iw,ix};
                    DRarraywx = ...
                        inprodDwx(nXbasisw, nWbasisw, ...
                                  nXbasisx, nWbasisx, ...
                                  Bvecx, Btenswx);
                end
                %  rescale the inner product
                DRarraywx = factorw*factorx*rhoVec(ivar)*DRarraywx/T;
                %  increment the inner product and its transpose for
                %  the appropriate location in DRarray
                if ivw == ivx && jvw == jvx
                    DRarray(indw,indw,indthw) = ...
                        DRarray(indw,indw,indthw) + 2*weighti*DRarraywx;
                else
                    DRarray(indw,indx,indthw) = ...
                        DRarray(indw,indx,indthw) +   weighti*DRarraywx;
                    DRarray(indx,indw,indthw) = ...
                        DRarray(indx,indw,indthw) + ...
                                weighti*permute(DRarraywx,[2,1,3]);
                end
            end
            %  partial derivatives wrt Wcoef for cross-products with D^m
            %  here x = ivar, Wbasisx is the constant basis, and
            %  Bvecx = 1;
            %  get the tensor vector for this pair of coefficients
            %  and derivatives
            %  compute inner product with respect to the ix dimension
            if funtypew
                % user code
                DRarraywi = inprod_Dbasis_Data2LD(Xbasisw,  Xbasisi,  ...
                                               coefStrw, coefStri, ...
                                               nderivw,  nderivi);
            else
                %  fda code
                Btenswi   = BtensorCell{ivar}{iw,nXtermi+1};
                DRarraywi = inprodDwi(nXbasisi, nXbasisw, nWbasisw, Btenswi);               
           end
            %  rescale the inner product
            DRarraywi = factorw*rhoVec(ivar)*DRarraywi/T;
            %  decrement the inner product and its transpose for
            % the appropriate location in DRarray
            if nWbasisw == 1
                DRarray(indi,indw,indthw) = DRarray(indi,indw,indthw) - ...
                    weighti*DRarraywi';
                DRarray(indw,indi,indthw) = DRarray(indw,indi,indthw) - ...
                    weighti*DRarraywi;
            else
                DRarraywit = permute(DRarraywi,[2,1,3]);
                DRarray(indi,indw,indthw) = DRarray(indi,indw,indthw) - ...
                    weighti*DRarraywit;
                DRarray(indw,indi,indthw) = DRarray(indw,indi,indthw) - ...
                    weighti*DRarraywi;
            end
        end
    end
end

else
    
    DRarray = [];
end

%  ------------------------------------------------------------------------

function Rmatwx = inprodwx(nXbasisw, nWbasisw, nXbasisx, nWbasisx, ...
                           Bvecw, Bvecx, Btenswx)
Rmatwx = zeros(nXbasisw,nXbasisx);
ncum   = cumprod([nWbasisx, nXbasisx, nWbasisw, nXbasisw]);
for i=1:nXbasisw
    for k=1:nXbasisx
        Rmatwx(i,k) = 0;
        for j=1:nWbasisw
            for l=1:nWbasisx
                ijkl = (i-1)*ncum(3) + ...
                       (j-1)*ncum(2) + ...
                       (k-1)*ncum(1) + l;
                Rmatwx(i,k) = Rmatwx(i,k) + ...
                    Btenswx(ijkl)*Bvecw(j)*Bvecx(l);
            end
        end
    end
end

%  ------------------------------------------------------------------------

function  DRarraywx = inprodDwx(nXbasisw, nWbasisw, nXbasisx, nWbasisx, ...
                               Bvecx, Btenswx)
DRarraywx = zeros(nXbasisw, nXbasisx, nWbasisw);
ncum = cumprod([nWbasisx, nXbasisx, nWbasisw, nXbasisw]);
for i=1:nXbasisw
    for j=1:nWbasisw
        for k=1:nXbasisx
            for l=1:nWbasisx
                ijkl = (i-1)*ncum(3) + ...
                       (j-1)*ncum(2) + ...
                       (k-1)*ncum(1) + l;
                if abs(ijkl) > eps
                    DRarraywx(i,k,j) = DRarraywx(i,k,j) + ...
                        Bvecx(l)*Btenswx(ijkl);
                end
            end
        end
    end
end

%  ------------------------------------------------------------------------

function  DRarraywi = inprodDwi(nXbasisi, nXbasisw, nWbasisw, Btenswi)
DRarraywi = zeros(nXbasisw, nXbasisi, nWbasisw);
ncum = cumprod([1, nXbasisi, nWbasisw, nXbasisw]);
for i=1:nXbasisw
    for j=1:nWbasisw
        for k=1:nXbasisi
            for l=1
                ijkl = (i-1)*ncum(3) + ...
                    (j-1)*ncum(2) + ...
                    (k-1)*ncum(1) + l;
                DRarraywi(i,k,j) = DRarraywi(i,k,j) + Btenswi(ijkl);
            end
        end
    end
end


