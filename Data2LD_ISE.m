function ISE = Data2LD_ISE(XbasisCell, modelCell, coefCell, coef, ...
                            Rmat, Smat,  nrep, nforce, rhoVec, AtensorCell)
%  Data2LD ... stands for "Data to Linear Dynamics"
%  Data2LD_ISE computes the value of the penalty term, the integrated
%  squared difference between the right and left sides of a linear
%  differential equation of the form
%      D^m x_i = sum_k^d sum_j^{m_k} beta_{kj} D^{j-1} x_k + 
%                sum_{\ell}^{M_i} \alpha_{\ell,i} u_{\ell,i}, 
%      i=1,...,d
%  where
%  and where each coefficient is expanded in terms of its own number of
%  B-spline basis functions:
%      \beta_{ij}(t)      = \bbold_{ij}'     \phibold_{ij}(t),
%      \alpha_{\ell,i}(t) = \abold_{\ell,i}' \psibold_{\ell,i}(t)
%  This version inputs AtensorCell as an argument.
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
%  Arguments:
%
%  XBASISCELL ... A functional data object or a BASIS object.  If so, the 
%               smoothing parameter LAMBDA is set to 0.
%
%  XBASISCELL ... A functional data object or a BASIS object.  If so, the 
%               smoothing parameter LAMBDA is set to 0.
%
%  MODELCELL...  A cell aray of length NVAR. Each cell contains a 
%                struct object with members:              
%                XCell ... cell array of length number of homogeneous terms
%                          Each cell contains a struct object with members:
%                          WfdPar    ... a fdPar object for the coefficient
%                          variable   ... the index of the variable
%                          derivative ... the order of its derivative
%                          npar ... if coefficient estimated, its location
%                                   in the composite vector 
%                FCell ... cell array of length number of forcing terms
%                          Each cell contains a struct object with members:
%                          AfdPar ... an fdPar object for the coefficient
%                          Ufd    ... an fd object for the forcing function
%                          npar ... if coefficient estimated, its location
%                                   in the composite vector 
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
%
%  COEF      ... The coefficient matrix
%
%  RMAT      ... The penalty matrix for the homogeneous part of L
%
%  SMAT      ... The penalty matrix for the forcing part of L
%
%  NREP      ... The number of replications
%
%  NFORCE    ... The vector containing the number of forcing functions
%                per variable.
%
%  ATENSORCELL ...

%  Last modified 23 December 2016

%  ------------------------------------------------------------------------
%                         Set up analysis
%  ------------------------------------------------------------------------

if nargin < 9
    error('Less than nine parameters supplied.');
end

%  ------------------------------------------------------------------------
%                    Check modelCell 
%  ------------------------------------------------------------------------

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
%  Compute unpenalized error integrated squares, ISE, the penalty term
%  loop through number of replications NREP
%  ------------------------------------------------------------------------

ISE1 = zeros(nvar,1);
ISE2 = zeros(nvar,1);
ISE3 = zeros(nvar,1);
for irep=1:nrep
    coefi = coef(:,irep);
    for ivar=1:nvar
        %  the homogeneous component
        ISE1i = coefi'*Rmat*coefi;
        modelStructi = modelCell{ivar};
        weighti = modelStructi.weight;
        nforcei  = nforce(ivar);
        if nforcei > 0
            %  the homogeneous/forcing component
            ISE2i = 2*coefi'*Smat(:,irep);
            ISE3i = 0;
            for jforce=1:nforcei
                % the forcing component
                modelStructij = modelStructi.FCell{jforce};
                ncoefj   = modelStructij.ncoef;
                coefStrj = coefCell{ncoefj};
                Avecj    = coefStrj.parvec;
                nAbasisj = length(Avecj);
                factorj  = modelStructij.factor;
                Ufdj     = modelStructij.Ufd;
                Ubasisj  = getbasis(Ufdj);
                Ucoefj   = getcoef(Ufdj);
                nUbasisj = getnbasis(Ubasisj);
                funtypej = isstruct(coefStrj.fun);
                for kforce = 1:nforcei
                    modelStructik = modelStructi.FCell{kforce};
                    ncoefk   = modelStructik.ncoef;
                    coefStrk = coefCell{ncoefk};
                    Aveck    = coefStrk.parvec;
                    factork  = modelStructik.factor;
                    nAbasisk = length(Aveck);
                    Ufdk     = modelStructik.Ufd;
                    Ubasisk  = getbasis(Ufdk);
                    Ucoefk   = getcoef(Ufdk);
                    nUbasisk = getnbasis(Ubasisk);
                    funtypek = isstruct(coefStrk.fun);
                    if funtypej || funtypek
                        ISE3i = inprod_basis_Data2LD(Ufdk,  Ufdj, ...
                                          coefStrk, coefStrj, ...
                                          0,   0);
                    else
                        ISE3i = 0;
                        Atensijk = AtensorCell{ivar}{jforce,kforce};
                        ncum     = cumprod([nAbasisk, nUbasisk, ...
                            nAbasisj, nUbasisj]);
                        for i=1:nUbasisj
                            for j=1:nAbasisj
                                for k=1:nUbasisk
                                    for l=1:nAbasisk
                                        ijkl = (i-1)*ncum(3) + ...
                                               (j-1)*ncum(2) + ...
                                               (k-1)*ncum(1) + l;
                                        ISE3i = ISE3i + ...
                                            Ucoefj(i,irep)* ...
                                            Avecj(j)*       ...
                                            Ucoefk(k,irep)* ...
                                            Aveck(l)*Atensijk(ijkl);
                                    end
                                end
                            end
                        end
                    end
                    ISE3i = factorj*factork*rhoVec(ivar)*ISE3i/T;
                end
            end
        else
            ISE2i = 0;
            ISE3i = 0;
        end
        ISE1(ivar) = ISE1(ivar) + weighti*ISE1i;
        ISE2(ivar) = ISE2(ivar) + weighti*ISE2i;
        ISE3(ivar) = ISE3(ivar) + weighti*ISE3i;
    end
end
ISE  = (ISE1 + ISE2 + ISE3)/nrep;

