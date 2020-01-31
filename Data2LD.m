function [MSE, DpMSE, D2ppMSE, XfdCell, ...
          df, gcv, ISE, Var_theta, Dcoef, Rmat, Smat, y2cMap] = ...
      Data2LD(yCell, XbasisCell, modelCell, coefCell, rhoVec, loadTensor)
%  Data2LD ... stands for "Data to Linear Dynamics"
%  It approximates the data in argument YCELL by one or smooth 
%  functions x_i, i=1,...,d.  This approximation is defined by a set
%  of linear differential or algebraic equations defined by a set of
%  parameters some of which require estimation from the data. 
%
%  Because Data2LD requires the computation of integrals of four-way 
%  tensors, it is not appropriate for very large numbers of basis functions
%  for either variables or coefficient functions multiplying variables.
%  It will work successfully for numbers of basis functions in the low 
%  hundreds, but beyond that will require more memory that Matlab can
%  provide.  If Data2LD can't handle the job, CollocInfer, which does
%  not work with these tensors, can.
%
%  The approximation minimizes the sum of squared residuals, expressed as 
%  follows in Latex notation:
%
%    H(\theta) = \sum_i^d \sum_j^n \sum_\ell^N [y_{ij \ell} - x_i(t_j)]^2
%
%  where:
%  i    = 1,...,d indexes equations in a system differential and algebraic
%                 equations.
%  j    = 1,...,n indexes times of observation of a variable
%  \ell = 1,...,N indexes replications of observations
%
%  But there is additional flexibility not captured in this expression:
%  1.  Only a subset of the variables may be observed, so that not all 
%      values of index i are actually used.
%  2.  The number and location of times of observation t_j  can vary
%      from one observed variable to another.
%  using a roughness penaltylinear differential operator that depends
%  on unknown parameters in cell array COEFCELL, which is described below.
%
%  The fitting functions x_i(t) in turn, here assumed to be defined over
%  the interval [0,T], are defined by basis function expansions:
%
%         x_i(t_j) = \sum_k^K_i c_{ik}(\theta|\rho) \phi_{ik}(t_j)
%
%  where \phi_{ik}(t) is the kth function in a system of K_i basis 
%  functions used to approximate the ith variable.
%  The number of K_i basis functions and the type of basis function system
%  can vary from one variable to another.  This information is contained
%  in argument XBASISCELL described below.
%
%  The coefficients c_{ik}(\theta|\rho) defining the smoothing functions  
%  are functions of the unknown parameters in vector \theta that define
%  the differential and algebraic equations and that require estimation
%  from the data.  The smoothing parameter $\rho is a value in the interval
%  [0,1).
%  The coefficient functions c_{ik}(\theta|\rho) minimize the inner 
%  least squares criterion, expressed here for simplicity for a single  
%  variable or equation:
%
%    J(c|\theta) = (1-\rho) \sum_j^n [y_j - c_k^K \phi_k(t_j)]^2/n +
%                  \rho \int_0^T {L(\theta)x(t)]^2 dt/T
%
%  Linear differential operator L(\theta) is a linear differential or
%  algebraic equation rewritten as an operator by subtracting the right
%  side of the equation from the left, so that a solution x of the 
%  equation satisfies Lx = 0.
%  Each L in a system of equation depends on one or more parameters that 
%  need to be estimated and are contained in parameter vector \theta.
%
%  The linear differential equation is of the form
%      D^m x_i = sum_k^d sum_j^{m_k} beta_{kj}(t) D^{j-1} x_k + 
%                sum_f^{F_i} \alpha_{fi}(t) u_{f,i}, 
%      i=1,...,d,  f=1,...,F_i
%  where
%  and where each coefficient is expanded in terms of its own number of
%  B-spline basis functions:
%      \beta_{ij}(t)  = \bbold_{ij}' \phibold_{ij}(t),
%      \alpha_{fi}(t) = \abold_{fi}' \psibold_{fi}(t)
%
%  As smoothing parameter \rho increases toward its upper limit of 1,
%  the second roughness penalty term in J is more and more emphasized
%  and the fit to the data less and less emphasized, so that x_i is
%  required to approach a solution to its respective equation.
%
%  The highest order of derivative can vary from one equation to another,
%  and in particular can be larger than 1.
%  Each right side may contain contributions from all variables in the
%  equation.  Morever, these contributions may be from all permissible
%  derivatives of each equation, wher "permissible" means up to one less 
%  than the highest order in the variable.
%  In addition, each equation may be forced by a number of forcing terms
%  that can vary from equation to equation.
%
%  This version approximates the integrals in the penalty terms by using 
%  inprod_basis to compute the cross-product matrices for the  
%  \beta-coefficient basis functions and the corresponding derivative of 
%  the x-basis functions,and the cross-product matrices for the 
%  \alpha-coefficients and the corresponding U functions.  
%  These are computed upon the first call to Data2LD, and then retained 
%  for subsequent calls by using the persistent command.  This technique
%  for speeding up computaton is called memoization.  
%  See lines 224 and 352 to 388 for this code.
%
%  The structure of the model is defined in cell array MODELCELL, which is
%  described below.
%
%  This version disassociates coefficient functions from equation 
%  definitions to allow some coefficients to be used repeatedly and for
%  both homogeneous and forcing terms.  It requires an extra argument
%  COEFCELL that contains the coefficients and the position of their
%  coefficient vectors in vector THETA.
%
%  ------------------------------------------------------------------------
%
%  Arguments:
%
%  YCELL     ... A cell array of length NVAR.  Each cell contains in turn
%                a struct object with fields:
%                  "argvals" is a vector of length n_i of observation times
%                  "y" contains a matrix with n_i rows and NREP columns.
%                The number of columns must be the same for all variables,
%                except that, if a cell is empty, that variable is taken to 
%                be not observed.
%
%  XBASISCELL ... A cell array of length NVAR.  Each cell contains in turn
%                 a functional data object or a BASIS object.  
%
%  MODELCELL...  A cell aray of length NVAR. Each cell contains a 
%                struct object with members:              
%                XCell ... cell array of length number of homogeneous terms
%                          Each cell contains a struct object with members:
%                          variable   ... the index of the variable
%                          derivative ... the order of its derivative
%                          ncoef      ... if coefficient estimated, its 
%                                         location in the composite vector 
%                          factor     ... a scalar multiplier (def. 1)
%                FCell ... cell array of length number of forcing terms
%                          Each cell contains a struct object with members:
%                          Ufd    ... an fd object for the forcing function
%                          ncoef  ... if coefficient estimated, its 
%                                     location in the composite vector 
%                          factor ... a scalar multiplier (def. 1)
%                order     ... the highest order of derivative..
%                name      ... a  tag for the variable
%                weight    ... a positive weight 
%                nallXterm ... the number of homogeneous terms
%                nallFterm ... the number of forcing functions
%
%  COEFCELL  ... A cell array of length NCOEF containing struct objects
%                with fields:
%               parvec   ... a vector of parameters
%               estimate ... 0, held fixed, otherwise, estimated 
%               coeftype ... homogeneous or forcing
%               fun      ... functional basis, fd, or fdPar object, 
%                            or a struct object for a general function 
%                            with fields:
%                 fd     ... function handle for evaluating function
%                 Dfd    ... function handle for evaluating 
%                            partial derivative with respect to parameter
%                 more   ... object providing additional information for 
%                            evaluating coefficient function
%
%  RHOVEC    ... A vector of length NVAR containing values in [0,1].  
%                The data sums of squares are weighted by P and 
%                the roughness penalty by 1-P.
%
%  LOAD_TENSOR ... If nonzero, attempt to load the cell arrays
%                BtensorCell, BAtensorCell and AtensorCell.  These must
%                have set up before any call to Data2LD and saves as 
%                .mat files with names BtensorCell.mat, BAtensorCell.mat
%                and AtensorCell.mat.  
%                For information on how these are set up, see the functions 
%                Btensorfn, BAtensorfn and Atensorfn.
%
%  ------------------------------------------------------------------------
%
%  Output objects (d = number of equations, 
%                  NTHETA is the total number of estimated parameters):
%
%  MSE     ... A vector of length d containing the error sum of squares for
%              each variable divided by the number of replications N and by 
%              the number of observations n_i.  However, if a variable
%              is not measured, the value is 0.  
%              The objective function minimized by LSMOOTH is the sum of 
%              this vector.  
%  DMSE    ... A vector of length NTHETA that is the gradient of the
%              objective function minimized by LSMOOTH.
%  D2MSE   ... A square symmetric matrx of order NTHETA that contains
%              the expected values of the second partial derivatives of the
%              objective function minimized by LSMOOTH.
%  XFDCELL ... A cell array of length d containing functional parameter
%              objects of class fdPar for the estimated functions x_i(t).
%  DF      ... An equivalent degrees of freedom value 
%                   df = trace(2*YM - YM*YM) where YM is the matrix
%              Y2CMAP described below.
%  GCV     ... The generalized cross-validation measure.  The value of
%              \rho corresponding to the minimim of GCV across values of
%              smoothing parameter \rho is often chose for an automatic 
%              data-driven level of smoothing.
%  ISE     ... The sum across variables of the integrated squared value of
%              the differential operator.  This value multiplied by \rho
%              and divided by T, the width of the domain, is the second
%              term the objective function.
%  Note:  I decided to eliminate the output of variable ISE.  We can
%  easily approximate it after an analysis, and it requires Atensor,
%  which can take a long time to calculate if the basis for a forcing
%  term is large.  This version output ISE = [].
%  RMAT    ... A square symmetric matrix of order equal to the sum of the
%              coefficients in the basis function expansions of the 
%              variables.  This matrix defines the size of the second term
%              in the objective function.
%  SMAT    ... Either a vector of length equal to the order of RMAT if
%              there is only one replication, or a matrix with number of
%              columns equal to the number of replications NREP.
%  Y2CMAP  ... A matrix with number of rows equal to the total number of
%              coefficients in the basis expansions of variables and 
%              number of columns equal the total number of observations.
%              This matrix is the linear map from the data to the 
%              combined coefficients.

%  Last modified 10 February 2017

%  Set default arguments

if nargin <  6,  loadTensor = 0;  end

%  ------------------------------------------------------------------------
%     Determine whether tensors need computing: define these variables
%        as persistent, so that memoization is activated.
%  ------------------------------------------------------------------------

persistent BtensorCell BAtensorCell AtensorCell basismatCell
% persistent BtensorCell BAtensorCell basismatCell

%  ------------------------------------------------------------------------
%                    Check modelCell 
%  ------------------------------------------------------------------------

modelCell = modelcheck(modelCell, coefCell);

nvar = length(modelCell);

%  ------------------------------------------------------------------------
%                    Check coefCell 
%  ------------------------------------------------------------------------

[coefCell, ntheta] = coefcheck(coefCell);

if nargin < 5,  rhoVec = 0.5*ones(nvar,1);  end
    
%  ------------------------------------------------------------------------
%  Store the number of homogeneous and forcing functions for each variable.
%  along with number of estimated functions and number of coefficients
%  ------------------------------------------------------------------------

nhomog  = zeros(nvar,1);  %  vector of numbers of homogeneous terms
nforce  = zeros(nvar,1);  %  vector of numbers of forcing terms
nthetaH = 0;  %  total number of estimated homogeneous terms
nthetaF = 0;  %  total number of estimated forcing terms
for ivar=1:nvar
    modelStructivar = modelCell{ivar};
    nXterm = modelStructivar.nallXterm;
    %  process homogeneous terms
    if nXterm > 0
        nhomog(ivar) = nXterm;
        for iterm=1:nXterm
            XStri    = modelStructivar.XCell{iterm};
            ncoefi   = XStri.ncoef;
            coefStri = coefCell{ncoefi};
            if coefStri.estimate
                nthetaH = nthetaH + length(coefStri.parvec);
            end
        end
    end
    nFterm = modelStructivar.nallFterm;
    %  process forcing terms
    if nFterm > 0
        nforce(ivar) = nFterm;
        for iterm=1:nFterm
            FStri    = modelStructivar.FCell{iterm};
            ncoefi   = FStri.ncoef;
            coefStri = coefCell{ncoefi};
            if coefStri.estimate
                nthetaF = nthetaF + length(coefStri.parvec);
            end
        end
    end
end

%  ------------------------------------------------------------------------
%                      Check structure of YCELL
%  ------------------------------------------------------------------------

[nrep, nvec, dataWrd] = yCellcheck(yCell, nvar);

%  ------------------------------------------------------------------------
%                   Check structure of XBASISCELL
%  ------------------------------------------------------------------------

%  Retrieve basis object for each variable from XBASISCELL and install it
%  in Xbasis Cell.
%  And set up a vector NCOEFVEC containing number of coefficients used
%  for the expansion of each variable

%  check that XbasisCell is a cell array

if ~iscell(XbasisCell)
    error('XBASISCELL is not a cell array.');
end

%  check length of XbasisCell

if length(XbasisCell) ~= nvar
    error('XBASISCELL is not of length NVAR.');
end

%  check that cells contain basis objects

errwrd = 0;
ncoefvec = zeros(nvar,1);
for ivar=1:nvar
    Xbasisi = XbasisCell{ivar};
    if ~isa_basis(Xbasisi)
        warning(['XBASIS is not a BASIS object for variable ', ...
                 num2str(ivar),'.']);
        errwrd = 1;
    else
        ncoefvec(ivar) = getnbasis(Xbasisi);
        XbasisCell{ivar} = Xbasisi;
    end
end
if errwrd
    error('One or more terminal error encountered in XBASISCELL.');
end

%  check length of rhoVec

if length(rhoVec) ~= nvar
    error('RHOVEC not of length NVAR.');
end

%  check that rho values are within [0,1).

for ivar=1:nvar
    if rhoVec(ivar) >= 1 || rhoVec(ivar) < 0
        error(['P is not in [0,1) for variable ',num2str(ivar),'.']);
    end
end

%  ------------------------------------------------------------------------
%  Compute four-way tensor products of D^j X-basis functions and 
%  W-basis functions, j=0,...,nderivvec if not already set up.  
%  If already set up, load these tensor products.
%  ------------------------------------------------------------------------

%  tensor products involving homogeneous terms only

if isempty(BtensorCell)
    if loadTensor
        load BtensorCell
        BtensorCellcheck(BtensorCell, modelCell);
    else
        tic;
        BtensorCell = Btensorfn(XbasisCell, modelCell, coefCell);
        toc
        save BtensorCell BtensorCell
    end
end

%  tensor products involving homogeneous and forcing terms

if isempty(BAtensorCell) && any(nforce > 0);
    if loadTensor
        load BAtensorCell
        BAtensorCellcheck(BAtensorCell, modelCell);
    else
        tic;
        BAtensorCell = BAtensorfn(XbasisCell, modelCell, coefCell);
        toc
        save BAtensorCell BAtensorCell
    end
end

%  tensor products involving forcing terms only

if isempty(AtensorCell) && any(nforce > 0);
    if loadTensor
        load AtensorCell
        AtensorCellcheck(modelCell);
    else
        tic;
        AtensorCell = Atensorfn(modelCell, coefCell);
        toc
        save AtensorCell AtensorCell
    end
end

%  ------------------------------------------------------------------------
%        Compute cell array basismatCell  and coefficient matrix Bmat
%  ------------------------------------------------------------------------

%  Set up cell array of matrices of basis function values for variables
%  having observations.
%  Each matrix has number of rows equal to the number of observation times
%  and number of columns equal to the number K_i of basis functions

basismatCell = cell(nvar,1);
for ivar=1:nvar
    if dataWrd(ivar)
        yStructi  = yCell{ivar};
        Xbasisi   = XbasisCell{ivar};
        basismati = eval_basis(yStructi.argvals, Xbasisi);
        basismatCell{ivar} = basismati;
    end
end

%  Set up square matrices of order K_i of basis value inner products,
%  each multiplied by the respective rho value and divided by the number
%  of sampling points.  These are needed only for measured variables.

nsum     = sum(nvec);
ncoefsum = sum(ncoefvec);
Bmat     = zeros(ncoefsum,ncoefsum);
basismat = zeros(nsum,ncoefsum);
ymat     = zeros(nsum,nrep);
m2 = 0;
n2 = 0;
for ivar=1:nvar
    modelStructi = modelCell{ivar};
    weighti  = modelStructi.weight;
    m1  = m2 + 1;
    m2  = m2 + ncoefvec(ivar);
    ind = m1:m2;
    if dataWrd(ivar)
        n1   = n2 + 1;
        n2   = n2 + nvec(ivar);
        indn = n1:n2;
        yStructi      = yCell{ivar};
        ymat(indn,:)  = yStructi.y;
        basismati     = basismatCell{ivar};
        Bmat(ind,ind) = weighti*(1-rhoVec(ivar)).*basismati'*basismati./ ...
                        nvec(ivar);
        basismat(indn,ind) = basismati;
    end
end

%  ------------------------------------------------------------------------
%      Compute roughness penalty matrices Rmat(theta) and Smat(theta)
%      and their derivatives with respect to estimated parameters
%  ------------------------------------------------------------------------

%  Matrices R and DR

if nthetaH > 0
    [Rmat, DRarray] = Data2LD_R(XbasisCell, modelCell, coefCell, ...
                                rhoVec, ntheta, BtensorCell);
else
    Rmat            = Data2LD_R(XbasisCell, modelCell, coefCell, ...
                                rhoVec, ntheta, BtensorCell);
    DRarray = [];
end

%  Matrices S and DS for variables having forcing functions

if nthetaF > 0
    [Smat, DSarray] = Data2LD_S(XbasisCell, modelCell, coefCell, ...
                                rhoVec, ntheta, BAtensorCell, nrep, nforce);
else
    Smat            = Data2LD_S(XbasisCell, modelCell, coefCell, ...
                                rhoVec, ntheta, BAtensorCell, nrep, nforce);
    DSarray = [];
end

%  ------------------------------------------------------------------------
%      Compute symmetric coefficient matrix Cmat for the linear equation
%  ------------------------------------------------------------------------

Cmat = Bmat + Rmat;

%  ------------------------------------------------------------------------
%                     Set up right side of equation
%  ------------------------------------------------------------------------

Dmat = zeros(ncoefsum,nrep);
m2 = 0;
for ivar=1:nvar
    modelStructi = modelCell{ivar};
    weighti  = modelStructi.weight;
    m1  = m2 + 1;
    m2  = m2 + ncoefvec(ivar);
    ind = m1:m2;
    if dataWrd(ivar)
        yStructi = yCell{ivar};
        basismati = basismatCell{ivar};
        yi = yStructi.y;
        ni = nvec(ivar);
        Dmat(ind,:) = weighti*(1-rhoVec(ivar)).*basismati'*yi/ni;
    end
end

%  ------------------------------------------------------------------------
%                  Compute coefficient matrix
%  ------------------------------------------------------------------------

[Cmatfac, rnkp1] = chol(Cmat);
if rnkp1 ~= 0
    disp(['The matrix linear equation coefficient matrix CMAT ', ...
          'is not positive definite.']);
      
    eigvals = eig(Cmat);
    disp('Eigenvalues of coefficient matrix:')
    disp(eigvals')
    
    error(['The coefficients for the variable approximations ', ...
           'cannot be computed.']);
end

%  compute coefficient matrix or array

if isempty(Smat)
    coef = Cmatfac\(Cmatfac'\Dmat);
else
    coef = Cmatfac\(Cmatfac'\(Dmat-Smat));
end

%  ------------------------------------------------------------------------
%                     Compute summary values
%  ------------------------------------------------------------------------

%  The composite matrix y2cMap mapping the composite data vector into
%  the composite coefficient vector

y2cFac = basismat';
m2 = 0;
n2 = 0;
for ivar=1:nvar
    modelStructi = modelCell{ivar};
    weighti  = modelStructi.weight;
    m1  = m2 + 1;
    m2  = m2 + ncoefvec(ivar);
    ind = m1:m2;
    if dataWrd(ivar)
        n1   = n2 + 1;
        n2   = n2 + nvec(ivar);
        indn = n1:n2;
        y2cFac(ind,indn) = ...
             weighti*(1-rhoVec(ivar))*y2cFac(ind,indn)/nvec(ivar);
    end
end
y2cMap = basismat*(Cmat\y2cFac);

%  Use y2cMap to compute a equivalent degrees of freedom measure

df = trace(2*y2cMap - y2cMap*y2cMap');

%  Set up the functional parameter (fdPar) data objects for the variables

XfdCell = cell(nvar,1);
m2 = 0;
for ivar=1:nvar
    m1  = m2 + 1;
    m2  = m2 + ncoefvec(ivar);
    ind = m1:m2;
    Xbasisi = XbasisCell{ivar};
    Xfdobji = fd(coef(ind,:),Xbasisi);
    XfdCell{ivar} = fdPar(Xfdobji);
end

%  ------------------------------------------------------------------------
%  Compute the vector of unpenalized mean square errors, MSE_i, 
%  the sum of which is the outer objective function H(\theta|\rho).
%  Each MSE_i is normalized by dividing by NREP and by the n_i's.
%  ------------------------------------------------------------------------

xmat = basismat*coef;

MSE  = zeros(nvar,1);  %  vector of MSE_i's
SSEtot = 0;  %  total un-normalized sum of squared errors
m2 = 0;
for ivar=1:nvar
    modelStructi = modelCell{ivar};
    weighti  = modelStructi.weight;
    if dataWrd(ivar)
        m1 = m2 + 1;
        m2 = m2 + nvec(ivar);
        %  fit to data for measured variable i
        xmati  = xmat(m1:m2,:);
        ymati  = yCell{ivar}.y;
        %  residual vector 
        rmati  = ymati - xmati;
        %  un-normalized error sum of squares
        SSEi   = sum(sum(rmati.^2));
        SSEtot = SSEtot + weighti*SSEi;
        %  mean squared error
        MSE(ivar) = SSEi/nrep/nvec(ivar);
    else
        MSE(ivar) = [];
    end
end

%  compute residual variance

Rvar = SSEtot/nsum;

%  compute GCV, generalized cross-validation index

if df < nsum
    gcv = Rvar/((nsum - df)/nsum)^2;
else
    gcv = NaN;
end

%  ------------------------------------------------------------------------
%  Compute unpenalized error integrated squares, ISE, the roughness penalty 
%  term.  ISE is a vector of length NVAR containing these terms for each
%  variable.
%  ------------------------------------------------------------------------

if nargout > 6
    ISE = Data2LD_ISE(XbasisCell, modelCell, coefCell, coef, ...
                       Rmat, Smat, nrep, nforce, rhoVec, AtensorCell);
else
    ISE = [];
end

%  ------------------------------------------------------------------------
%       Compute total derivative of MSE wrt theta if required
%  ------------------------------------------------------------------------

if nargout == 1 || ntheta == 0 
    DpMSE     = [];
    D2ppMSE   = [];
    Var_theta = [];
    return;
end

%  ---------  Compute the total derivative  ----------------

%  Compute the partial derivatives of the coefficients with respect to the
%  estimated parameters,  dc/dtheta

Dcoef = zeros(ncoefsum,ntheta,nrep);
for itheta=1:ntheta
    if nthetaH > 0
        DRmati = squeeze(DRarray(:,:,itheta));
        if nthetaF > 0
            for irep=1:nrep
                DRi = -DRmati*coef(:,irep);
                DSi = -DSarray(:,irep,itheta);
                Dcoef(:,itheta,irep) = Cmat\(DRi + DSi);
            end
        else
            for irep=1:nrep
                DRi = -DRmati*coef(:,irep);
                Dcoef(:,itheta,irep) = Cmat\DRi;
            end
        end
    else
        if nthetaF > 0
            for irep=1:nrep
                DSi = -DSarray(:,irep,itheta);
                Dcoef(:,itheta,irep) = Cmat\DSi;
            end
        end
    end
end

%  ------------------------------------------------------------------------
%              Compute the total theta-gradient of H 
%  ------------------------------------------------------------------------

xmat    = basismat*coef;
DpMSE   = zeros(ntheta,1);
D2ppMSE = zeros(ntheta, ntheta);
D2pyMSE = zeros(ntheta,nsum);

m2 = 0;
for ivar=1:nvar
    if dataWrd(ivar)
        m1 = m2 + 1;
        m2 = m2 + nvec(ivar);
        weighti   = modelStructi.weight;
        basismati = basismat(m1:m2,:);
        xmati     = xmat(m1:m2,:);
        yVeci     = yCell{ivar}.y;
        rmati     = yVeci - xmati;
        for irep=1:nrep
            Dcoefi = squeeze(Dcoef(:,:,irep));
            BasDcoefi = basismati*Dcoefi;
            DpMSE   = DpMSE   - 2*weighti*(rmati(:,irep)'* ...
                 BasDcoefi)'/nrep/nvec(ivar);
            D2ppMSE = D2ppMSE + 2*weighti* ...
                (BasDcoefi'*BasDcoefi)/nrep/nvec(ivar);
            D2pyMSE(:,m1:m2) = D2pyMSE(:,m1:m2) - 2.*weighti.* ...
                 BasDcoefi'/nrep/nvec(ivar);
        end
    end
end

if nargout > 7
    DpDy      = -D2ppMSE\D2pyMSE;
    sigma_sq  = SSEtot./(nsum - df);
    Var_theta = sigma_sq.*DpDy*DpDy';
else
    Var_theta = [];
end

