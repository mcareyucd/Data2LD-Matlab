function  [ss, iter] = inprod_Dbasis_Data2LD(fdobj1, fdobj2,     ...
                                             coefStr1, coefStr2, ...
                                             Lfdobj1, Lfdobj2,   ...
                                             EPS, JMAX, JMIN)
%  INPROD_DBASIS  Computes the three-way tensor of inner products
%  where the first two dimensions are the number of basis functions 
%  of fd functions in FDOBJ1 and FDOBJ2 respectively; 
%  and the third dimension is the partial derivatives of the first 
%  coefficient function with respect to its defining parameter vector.
%  The  integration is approximated using Romberg integration with the 
%  trapezoidal rule.
%
%  The coefficient function arguments coefStr1, coefStr2 can be either 
%  a struct object or a functional basis object.
%
%  Arguments:
%  FDOBJ1 and FDOBJ2     ... these are either fd or basis objects.
%  COEFSTR1 and COEFSTR2 ... Struct objects defining the coefficient 
%             functions containing the following fields
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
%  LFDOBJ1 and LFDOBJ2 ...  differential operators for inner product for
%                    BASIS1 and BASIS2, respectively
%  EPS     ...  A convergence criterion, defaults to 1e-4.
%  JMAX    ...  Maximum number of Richardson extrapolation iterations.
%            Defaults to 15.
%  JMIN    ...  Minimum number of Richardson extrapolation iterations.
%            Defaults to 5.
%
%  Return:
%  A NBASIS1 by NBASIS2 matrix SS of inner products for each possible pair
%  of basis functions.

%  Last modified 9 December 2016

%  Check that first four arguments are supplied

if nargin < 4,  error('Less than four arguments are supplied.');  end

%  set up default values of arguments

if nargin < 9 || isempty(JMIN),    JMIN = 5;             end
if nargin < 8 || isempty(JMAX),    JMAX = 16;            end
if nargin < 7 || isempty(EPS),     EPS  = 1E-5;          end
if nargin < 6 || isempty(Lfdobj2), Lfdobj2 = int2Lfd(0); end
if nargin < 5 || isempty(Lfdobj1), Lfdobj1 = int2Lfd(0); end

%  check LFDOBJ1 and LFDOBJ2

Lfdobj1 = int2Lfd(Lfdobj1);
Lfdobj2 = int2Lfd(Lfdobj2);

%  Determine where FDOBJ1 and FDOBJ2 are basis or fd objects, define
%  BASIS1 and BASIS2, and check for common range

errwrd = 0;
if     isa_basis(fdobj1)
    basis1  = fdobj1;
    nbasis1 = getnbasis(basis1) - length(getdropind(basis1));
    ndim1   = nbasis1;
    fdtype1 = 0;
elseif isa_fd(fdobj1)
    basis1  = getbasis(fdobj1);
    nbasis1 = getnbasis(basis1) - length(getdropind(basis1));
    ndim1   = size(getcoef(fdobj1),2);
    fdtype1 = 1;
else
    errwrd = 1;
    disp('First argument is neither a basis object nor an fd object.')
end

if     isa_basis(fdobj2)
    basis2  = fdobj2;
    nbasis2 = getnbasis(basis2) - length(getdropind(basis2));
    ndim2   = nbasis2;
    fdtype2 = 0;
elseif isa_fd(fdobj2)
    basis2  = getbasis(fdobj2);
    nbasis2 = getnbasis(basis2) - length(getdropind(basis2));
    ndim2   = size(getcoef(fdobj2),2);
    fdtype2 = 1;
else
    errwrd = 1;
    disp('Second argument is neither a basis object nor an fd object.')
end

if errwrd
    err('Terminal error encountered.');
end

%  get coefficient vectors

bvec1 = coefStr1.parvec;
bvec2 = coefStr2.parvec;

%  set up beta functions

if isstruct(coefStr1.fun)
    type1 = 1;
    betaDfd1 = coefStr1.fun.Dfd;
    more1    = coefStr1.fun.more;
else
    type1   = 0;
    fdobj = coefStr1.fun;
    if isa_basis(fdobj) 
        betabasis1 = fdobj; 
    end
    if isa_fd(fdobj)    
        betabasis1 = getbasis(fdobj);  
    end
    if isa_fdPar(fdobj) 
        betabasis1 = getbasis(getfd(fdobj));
    end
end

if isstruct(coefStr2.fun)
    type2 = 1;
    betafd2 = coefStr2.fun.fd;
    more2   = coefStr2.fun.more;
else
    type2   = 0;
    fdobj = coefStr2.fun;
    if isa_basis(fdobj) 
        betafd2 = fd(bvec2,fdobj); 
    end
    if isa_fd(fdobj)    
        betafd2 = putcoef(fdobj, coefStr2.parvec);  
    end
    if isa_fdPar(fdobj) 
        betafd2 = putcoef(getfd(fdobj),coefStr2.parvec);
    end
end

npar = length(bvec1);

%  check for any knot multiplicities in either argument

knotmult = [];

%  check first functional object for knot multiplicities

if strcmp(getbasistype(basis1),'bspline')
    % Look for knot multiplicities in first basis
    params1  = getbasispar(basis1);
    nparams1 = length(params1);
    for i=2:nparams1
        if params1(i) == params1(i-1) || ...
           getnbasis(basis1) == length(nparams1) + 1
            knotmult = [knotmult, params1(i)];
        end
    end
end

%  check second functional object for knot multiplicities

%  univariate case
if strcmp(getbasistype(basis2),'bspline')
    % Look for knot multiplicities in first basis
    params2  = getbasispar(basis2);
    nparams2 = length(params2);
    for i=2:nparams2
        if params2(i) == params2(i-1) || ...
           getnbasis(basis2) == length(nparams2) + 1
            knotmult = [knotmult, params2(i)];
        end
    end
end
   
%  Set up RNGVEC defining subinvervals if there are any
%  knot multiplicities.

rng = getbasisrange(basis1);
if ~isempty(knotmult)
    knotmult = sort(unique(knotmult));
    knotmult = knotmult(knotmult > rng(1) & knotmult < rng(2));
    rngvec = [rng(1), knotmult, rng(2)];
else
    rngvec = rng;
end

%  -----------------------------------------------------------------
%                   loop through sub-intervals
%  -----------------------------------------------------------------

%  outer loop is over inter-multiple-knot ntervals

nrng = length(rngvec);
for irng = 2:nrng
    rngi = [rngvec(irng-1),rngvec(irng)];
    %  change range so as to avoid being exactly on
    %  multiple knot values
    if irng > 2
        rngi(1) = rngi(1) + 1e-10;
    end
    if irng < nrng
        rngi(2) = rngi(2) - 1e-10;
    end
    
  %  set up first iteration

    width = rng(2) - rng(1);
    JMAXP = JMAX + 1;
    h     = ones(JMAXP,1);
    h(2)  = 0.25;
    tnm   = 0.5;
    s = reshape(zeros(JMAXP*ndim1*ndim2*npar,1), ...
               [JMAXP,ndim1,ndim2,npar]);
    x = rngi;
    %  the first iteration uses just the endpoints
    %  For first argument:
    %  matrix of partial derivative values of first coefficient
    if type1
        betamat1 = betaDfd1(x, bvec1, more1);
    else
        betamat1  = eval_basis(x, betabasis1);
    end
    %  basis or fd matrix for first argument
    if isa_basis(fdobj1)
        basismat1 = eval_basis(x, basis1, Lfdobj1);
    else
        basismat1 = eval_fd(x, fdobj1, Lfdobj1);
    end
    %  vector of values of second coefficient
    if type2
        betavec2 = betafd2(x, bvec2, more2);
    else
        betavec2 = eval_fd(x, betafd2);
    end
    %  For second argument:
    %  basis or fd matrix for second argument
    if isa_basis(fdobj2)
        basismat2 = eval_basis(x, basis2, Lfdobj2);
    else
        basismat2 = eval_fd(x, fdobj2, Lfdobj2);
    end
    
    %  loop through parameters to define three-way tensor
    temp2 = basismat2.*repmat(betavec2,1,nbasis2);
    for k=1:npar
        temp1 = basismat1.*repmat(betamat1(:,k),1,nbasis1);
        chs = width.*temp1'*temp2./2;
        s(1,:,:,k) = chs;
    end

    %  now iterate to convergence

    for iter = 2:JMAX
        tnm = tnm.*2;
        del = width./tnm;
        x   = rng(1)+del/2:del:rng(2);
        %  For first argument:
         %  matrix of partial derivative values of first coefficient
        if type1
            betamat1 = betaDfd1(x, bvec1, more1);
        else
            betamat1  = eval_basis(x, betabasis1);
        end
        %  basis or fd matrix for first argument
        if isa_basis(fdobj1)
            basismat1 = eval_basis(x, basis1, Lfdobj1);
        else
            basismat1 = eval_fd(x, fdobj1, Lfdobj1);
        end
        %  vector of values of second coefficient
        if type2
            betavec2 = betafd2(x, bvec2, more2);
        else
            betavec2 = eval_fd(x, betafd2);
        end
        %  For second argument:
        %  basis or fd matrix for second argument
        if isa_basis(fdobj2)
            basismat2 = eval_basis(x, basis2, Lfdobj2);
        else
            basismat2 = eval_fd(x, fdobj2, Lfdobj2);
        end       
        %  loop through parameters to define three-way tensor
        temp2 = basismat2.*repmat(betavec2,1,nbasis2);
        for k=1:npar
            temp1 = basismat1.*repmat(betamat1(:,k),1,nbasis1);
            chs = width.*temp1'*temp2./tnm;
            chsold = reshape(squeeze(s(iter-1,:,:,k)),size(chs));
            s(iter,:,:,k) = (chsold + chs)./2;
        end
        if iter >= 5
            ind = (iter-4):iter;
            ya = s(ind,:,:,:);
            xa = h(ind);
            absxa = abs(xa);
            [absxamin, ns] = min(absxa);
            cs = ya;
            ds = ya;
            y  = squeeze(ya(ns,:,:,:));
            ns = ns - 1;
            for m = 1:4
                for i = 1:(5-m)
                    ho      = xa(i);
                    hp      = xa(i+m);
                    w       = (cs(i+1,:,:,:) - ds(i,:,:,:))./(ho - hp);
                    ds(i,:,:,:) = hp.*w;
                    cs(i,:,:,:) = ho.*w;
                end
                if 2*ns < 5-m
                    dy = squeeze(cs(ns+1,:,:,:));
                else
                    dy = squeeze(ds(ns,:,:,:));
                    ns = ns - 1;
                end
                y = y + dy;
            end
            ss = reshape(y, ndim1, ndim2, npar);
            errval = max(abs(reshape(dy,ndim1*ndim2*npar,1)));
            ssqval = max(abs(reshape(ss,ndim1*ndim2*npar,1)));
            if all(ssqval > 0)
                crit = errval./ssqval;
            else
                crit = errval;
            end
            if crit < EPS && iter >= JMIN
                return
            end
        end
        s(iter+1,:,:,:) = s(iter,:,:,:);
        h(iter+1)   = 0.25.*h(iter);
    end
%     disp(['No convergence after ',num2str(JMAX), ...
%           ' steps in INPROD_DBASIS.']);

end

