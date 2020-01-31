function [thetastore, dfstore, gcvstore, coefCell_optCell] = ...
         Data2LD_Opt(yCell, XbasisCell, modelCell, coefCell, rhoMat, ...
                     conv, iterlim, dbglev, active, parMap, load_tensor)
%  Data2LD ... stands for "Data to Linear Dynamics"
%  Data2LD_OPT optimizes a parameter vector theta defining a 
%    linear differential operator object used to smooth a set of data.
%
%  Arguments:
%  YCELL    ... an array containing values of curves
%               If the array is a matrix, rows must correspond to argument
%               values and columns to replications, and it will be assumed
%               that there is only one variable per observation.
%               If Y is a three-dimensional array, the first dimension
%               corresponds to argument values, the second to replications,
%               and the third to variables within replications.
%               If Y is a vector, only one replicate and variable are 
%               assumed.
%  XBASISCELL ... A functional data object or a BASIS object.  If so, the 
%               smoothing parameter LAMBDA is set to 0.
%  MODELCELL...  A cell aray of length NVAR. Each cell contains a 
%                struct object with members:              
%                Xcell ... cell array of length number of homogeneous terms
%                          Each cell contains a struct object with members:
%                          WfdPar ... a fdPar object for the coefficient
%                          variable   ... the index of the variable
%                          derivative ... the order of its derivative
%                          npar ... if coefficient estimated, its location
%                                   in the composite vector 
%                Fcell ... cell arrau of length number of forcing terms
%                          Each cell contains a struct object with members:
%                          AfdPar ... an fdPar object for the coefficient
%                          Ufd    ... an fd object for the forcing function
%                          npar ... if coefficient estimated, its location
%                                   in the composite vector 
%                order     ... the highest order of derivative
%                weight    ... a positive weight 
%                name      ... a  tag for the variable
%                nallXterm ... the number of homogeneous terms
%                nallFterm ... the number of forcing functions
%  RHOMAT   ... A vector or matrix of values in [0,1].  
%               The data are weighted by rho and ...
%               the roughness penalty by 1-rho.
%               the length of the vector of the number of rows of the 
%               matrix must be equal to the number of variables.
%               If RHOMAT is a matrix, the number of columns corresponds
%               to the number of sequential optimization to perform,
%               each optimization using the parameter vector from the last
%               optimization.
%               This is designed to allow the user to conveniently step up
%               the P values from low levels where the optimization is 
%               rapid to higher values where the optimization is more
%               difficult unless the initial parameter vector is close to
%               the optimal vector.
%  CONV     ... One convergence criterion, or a vector of two criteria.
%               The first criterion is applied to the function change,
%               and the second is applied to the gradient norm.
%  ITERLIM  ... Maximum number of iterations allowed.
%  DBGLEV   ... An integer controlling amount of output per iteration.
%               Defaults to 1, which prints summary results for 
%               each iteration.
%  ACTIVE   ... Indices of constrained parameters that are to be estimated.
%               the others will be held fixed at their initial values.
%               Default:  all parameters
%  PARMAP   ... A rectangular matrix with number of rows equal to the
%               number of parameters to be estimated as defined in
%               BWTCELL, and number of columns equal to the number of 
%               parameters less the number of linear constraints on the
%               estimated parameters.  The columns of PARMAP must be
%               orthnormal so that PARMAP'*PARMAP is an identity matrix.
%               PARMAP'*THETA maps unconstrained parameters and the 
%               corresponding gradient into   constrained parameter space.
%               PARMAP*THETA  maps constrained parameters and the 
%               corresponding gradient into unconstrained parameter space.
%               PARMAP will usually be set up using the full QR 
%               decomposition of a linear constraint coefficient matrix A'
%               where the constraints are of the form A P = B, A and B 
%               being known matrices.  An example of such a constraint
%               that arises often is one where two estimated coefficients
%               are constrained to be equal.  For example, if a vvariable 
%               X involved in an equation the form a(x - x_0), where x_0
%               is a fixed set point or defined target level for variable 
%               X, then this would be set up as a_1 x + a_2 x_0, where
%               coefficients a_1 and a_2 are constrained to be equal in
%               magnitude but opposite in sign, or a_1 + a_2 = 0.  
%  LOAD_TENSOR ... If nonzero, attempt to load the cell arrays
%                BtensorCell, BAtensorCell and AtensorCell.  These must
%                have set up before any call to Data2LD and saves as 
%                .mat files with names BtensorCell.mat, BAtensorCell.mat
%                and AtensorCell.mat.  
%                For information on how these are set up, see the functions 
%                Btensorfn, BAtensorfn and Atensorfn.

%  Returns:
%
%  THETA_OPT    ... The optimal parameter values.
%  COEFCELL_OPT ... The parameter values installed in a cell array
%                   containing the coefficient fdPar objects

%  Last modified 15 June 2017 by Jim Ramsay

%  ------------------------------------------------------------------------
%                     Check input parameters
%  ------------------------------------------------------------------------

%  check that first five arguments are present.

if nargin <   5
    error('Less than five arguments supplied.');
end

%  Get the vector of parameters

theta  = BAwtcell2vec(modelCell, coefCell);
ntheta = length(theta);
npar   = length(theta);

%  the number of variables

nvar  = length(yCell);

%  Set default arguments

if nargin <  11 || isempty(load_tensor),  load_tensor  = 0;       end
if nargin <  10 || isempty(parMap),       parMap  = eye(npar);    end
if nargin <   8 || isempty(dbglev),       dbglev  =  0;           end
if nargin <   7 || isempty(iterlim),      iterlim = 20;           end
if nargin <   6 || isempty(conv),         conv    = 1e-6;         end

thetaCon = parMap'*theta;  

nparCon = length(thetaCon);
climit  = 1000*([-1; 1]*ones(1,nparCon));
if nargin <   9 || isempty(active)       
    active  =  1:nparCon;   
end

%  check rhoMat and get nopt

[nvartmp, nopt] = size(rhoMat);

if nvartmp ~= nvar
    error(['The first dimension of RHOMAT is not equal to ', ...
           'the number of variables'])
end

%   Cell arrays to contain results if more than a single set of rho's
%   are involved

thetastore = zeros(ntheta,nopt);
dfstore    = zeros(nopt,1);
gcvstore   = zeros(nopt,1);
coefCell_optCell = cell(nopt,1);

coefCell_opti = coefCell;

%  ------------------------------------------------------------------------
%                  loop through rho vectors
%  ------------------------------------------------------------------------

for iopt = 1:nopt
        
rhoVeci = rhoMat(:,iopt);

disp(['Rho = ',num2str(rhoVeci')])
if nopt > 1
    disp(' ')
    disp(['Rho = ',num2str(rhoVeci')])
end

%  compute initial criterion value, gradient and hessian
tic;
[fvec, grad, hessmat] = ...
       Data2LD(yCell, XbasisCell, modelCell, coefCell_opti, rhoVeci);
gradCon    = parMap'*grad;
hessmatCon = parMap'*hessmat*parMap;
toc;
f = sum(fvec);

norm = sqrt(mean(gradCon.^2));

if ntheta > 0
    
%  evaluate the initial update vector for correcting the initial theta

deltac = -gradCon;

%  initialize iteration status arrays

iternum = 0;
status = [iternum, f, norm];
if dbglev >= 1
    fprintf('\nIter.        Criterion   Grad Length\n')
    fprintf('%3.f %18.6f %12.6f\n', status);
end
iterhist = zeros(iterlim+1,length(status));
iterhist(1,:)  = status;
if iterlim == 0, return;  end

%  -------  Begin main iterations  -----------

MAXSTEPITER = 5;
MAXSTEP     = 1000;
trial       = 1;
reset       = 0;
linemat     = zeros(3,5);
thetaoldCon = thetaCon;
fold        = f;
gradoldCon  = gradCon;
dbgwrd      = dbglev >= 2;

%  ---------------  beginning of optimization loop  -----------

for iter = 1:iterlim
    iternum = iternum + 1;
    %  set logical parameters
    dblwrd = [0,0]; limwrd = [0,0]; ind = 0;  ips = 0;
    %  compute slope
    linemat(2,1) = sum(deltac.*gradoldCon);
    %  normalize search direction vector
    sdg          = sqrt(sum(deltac.^2));
    deltac       = deltac./sdg;
    linemat(2,1) = linemat(2,1)/sdg;
    % initialize line search vectors
    linemat(:,1:4) = [0; linemat(2,1); f]*ones(1,4);
    stepiter  = 0;
    if dbglev >= 2
        fprintf('                 %3.f %10.4f %12.6f %12.6f\n', ...
            [stepiter, linemat(:,1)']);
%  ------------------  start of code change August 2016  ------------------
        if dbglev > 2
            linestore = linemat(:,1)';
        end
%  -------------------  end of code change August 2016  -------------------
    end
    %  return with error condition if initial slope is nonnegative
    if linemat(2,1) >= 0
        if dbglev >= 2
            disp(['Initial slope ',num2str(linemat(2,1)),' nonnegative.']); 
        end
%  ------------------  start of code change August 2016  ------------------
        deltac       = -gradCon;
        linemat(2,1) = -sum(gradCon.^2);
        sdg          = sqrt(sum(deltac.^2));
        deltac       = deltac./sdg;
        linemat(2,1) = linemat(2,1)/sdg;
        linemat(:,1:4) = [0; linemat(2,1); f]*ones(1,4);
        linestore = linemat(:,1)';
        if dbglev > 2
            fprintf('                 %3.f %10.4f %12.6f %12.6f\n', ...
                    [stepiter, linemat(:,1)']);
        end
%  -------------------  end of code change August 2016  -------------------
    end
    %  return successfully if initial slope is very small
    if linemat(2,1) >= -min([1e-3,conv])
        if dbglev >= 2, disp('Initial slope too small'); end
        status = [iternum, f, norm];
        if dbglev >= 1
            fprintf('%3.f %18.6f %12.6f\n', status);
        end
        coefCellnew = coefCell;
        break;
    end
    %  first step set to trial
    linemat(1,5)  = trial;
    %  ------------  begin line search iteration loop  ----------
    thetanewCon = thetaCon;
    gradnewCon  = gradCon;
    hessmatnewCon = hessmatCon;
    % gradnewCon  = gradCon;
    for stepiter = 1:MAXSTEPITER
        %  check the step size and modify if limits exceeded
        [linemat(1,5), ind, limwrd] = ...
            stepchk(linemat(1,5), thetaCon, deltac, limwrd, ind, ...
                    climit, active, dbgwrd);
        if ind == 1 
            if dbglev >= 2
                fprintf('ind == 1\n');
            end
            thetanew    = theta;
            fnew        = f;
            thetanewCon = thetaCon;
            gradnewCon  = gradCon;
            coefCellnew = coefCell;
            break; 
        end % break of limit hit twice in a row
        if linemat(1,5) <= 1e-7
            %  Current step size too small ... terminate
            if dbglev >= 2
                fprintf('Stepsize too small: %15.7f\n', linemat(1,5));
            end
            thetanew   = theta;
            fnew       = f;
            hessmatnew = hessmat;
            break;
        end
        %  update parameter vector
        thetanewCon = thetaCon + linemat(1,5).*deltac;
        %  ---------  update function, gradient and hessian  -----------
        thetanew = parMap*thetanewCon;
        coefCellnew = BAwtvec2cell(thetanew, coefCell);
        [fvecnew, gradnew, hessmatnew] = ...
            Data2LD(yCell, XbasisCell, modelCell, coefCellnew, rhoVeci);
        fnew          = sum(fvecnew);
        gradnewCon    = parMap'*gradnew;
        hessmatnewCon = parMap'*hessmatnew*parMap;
        %  -------------------------------------------------------------
        linemat(3,5) = fnew;
        %  compute new directional derivative
        linemat(2,5) = sum(deltac.*gradnewCon);
        if dbglev >= 2
            fprintf('                 %3.f %10.4f %12.6f %12.6f\n', ...
                [stepiter, linemat(:,5)']);
%  ------------------  start of code change August 2016  ------------------
            if dbglev > 2
                linestore = [linestore; linemat(:,5)'];
            end
%  -------------------  end of code change August 2016  -------------------
        end
        %  compute next line search step, also testing for convergence
        [linemat, ips, ind, dblwrd] = ...
            stepit(linemat, ips, dblwrd, MAXSTEP);
        trial  = linemat(1,5);
        %  ind == 0 implies convergence
        if ind == 0 || ind == 1 || ind == 5, break; end
    end
%  ------------------  start of code change August 2016  ------------------
    if dbglev > 2
        figure(20)
        linesrchplot(linestore(:,1),linestore(:,2),linestore(:,3))
        title(['Iteration ',num2str(iter)])
        pause
    end
%  -------------------  end of code change August 2016  -------------------
    %  ------------  end line search iteration loop  ----------
    thetaCon   = thetanewCon;
    theta      = thetanew;
    f          = fnew;
    gradCon    = gradnewCon;
    hessmatCon = hessmatnewCon;
    %  test for function value made worse
    if f > fold
        %  Function value worse ... warn and terminate
        if dbglev >= 2
            fprintf('Criterion increased, terminating iterations.\n');
            fprintf('%10.4f %10.4f\n',[fold, f]);
        end
        %  reset parameters and fit
        thetaCon  = thetaoldCon;
        f      = fold;
%  ------------------  start of code change August 2016  ------------------
        deltac = -gradCon;
%  -------------------  end of code change August 2016  -------------------
        if dbglev > 2
            for i = 1:npar, fprintf('%10.4f%', theta(i)); end
            fprintf('\n');
        end
        if reset == 1
            %  This is the second time in a row that this
            %     has happened ...  quit
            if dbglev >= 2
                fprintf('Reset twice, terminating.\n');
                fprintf('Convergence not attained.\n')
            end
            %  return current status of optimization
            break;
        else
            reset = 1;
        end
    else
        %  function value has not increased,  check for convergence
        RMSgrad = sqrt(mean(gradCon.^2));
        if length(conv) > 1
            convtest = fold-f < conv(1) && RMSgrad < conv(2);
        else
            convtest = fold-f < conv;
        end
        if convtest
            norm   = sqrt(mean(gradCon.^2));
            status = [iternum, f, norm];
            if dbglev >= 1
                fprintf('%3.f %18.6f %12.6f\n', status);
                fprintf('Convergence reached\n')
            end
            %  return current status of optimization
            break;
        end
        %  update old parameter vectors and fit structure
        thetaoldCon = thetaCon;
        fold        = f;
        gradoldCon  = gradCon;
        hessmatoldCon  = parMap'*hessmatnew*parMap;
        %  update the line search direction vector
        deltac      = -hessmatCon\gradCon;
        reset       = 0;
    end
    norm   = sqrt(mean(gradCon.^2));
    status = [iternum, f, norm];
    iterhist(iter+1,:) = status;
    if dbglev >= 1
        fprintf('%3.f %18.6f %12.6f\n', status);
    end
    if dbglev >= 1 && iter == iterlim
        fprintf(...
            '\nMaximum iterations reached but convergence not attained.\n')
    end
end

%  ---------------  end of optimization loop  -----------

%  return current status of optimization

[SSEi, DSSEi, D2SSEi, XfdParCelli, dfi, gcvi] = ...
             Data2LD(yCell, XbasisCell, modelCell, coefCellnew, rhoVeci);

thetastore(:,iopt) = theta;
dfstore(iopt)      = dfi;
gcvstore(iopt)     = gcvi;
coefCell_optCell{iopt} = coefCellnew;

else
    
[SSEi, DSSEi, D2SSEi, XfdParCelli, dfi, gcvi] = ...
             Data2LD(yCell, XbasisCell, modelCell, coefCell, rhoVeci);

dfstore(iopt)      = dfi;
gcvstore(iopt)     = gcvi;
    
end


end
