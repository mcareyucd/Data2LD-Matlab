function [Res,Results_cell] = D2LD_ODE(data,XbasisCell,coef_ODE,...
    estimate,FCell,plots)
%  D2LD ... stands for "Data to Linear Dynamics"
%
%  It approximates the data in argument data by one or smooth
%  functions x_i, i=1,...,d.  This approximation is defined by a set
%  of linear differential equations defined by a set of
%  parameters some of which require estimation from the data.
%
%  The approximation minimizes the sum of squared residuals, expressed as
%  follows:
%
%    H(\theta) = \sum_i^d \sum_j^n \sum_\ell^N [y_{ij \ell} - x_i(t_j)]^2
%
%  where:
%  i    = 1,...,d indexes equations in a system differential
%                 equations.
%  j    = 1,...,n indexes times of observation of a variable
%  \ell = 1,...,N indexes replications of observations
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
%  in argument XbasisCell described below.
%
%  The coefficients c_{ik}(\theta|\rho) defining the smoothing functions
%  are functions of the unknown parameters in vector \theta that define
%  the differential equation and that require estimation
%  from the data.
%
%  The smoothing parameter $\rho is a value in the interval [0,1).
%  The coefficient functions c_{ik}(\theta|\rho) minimize the inner
%  least squares criterion, expressed here for simplicity for a single
%  variable or equation:
%
%    J(c|\theta) = (1-\rho) \sum_j^n [y_j - c_k^K \phi_k(t_j)]^2/n +
%                  \rho \int_0^T {L(\theta)x(t)]^2 dt/T
%
%  Linear differential operator L(\theta) is a linear differential equation
%  rewritten as an operator by subtracting the right
%  side of the equation from the left, so that a solution x of the
%  equation satisfies Lx = 0.
%
%  Each L in a system of equation depends on one or more parameters that
%  need to be estimated and are contained in parameter vector \theta.
%
%  The linear differential equation is of the form
%      D^m x_i = sum_k^d sum_j^{m_k} beta_{kj}(t) D^{j-1} x_k +
%                sum_f^{F_i} \alpha_{fi}(t) u_{f,i},
%      i=1,...,d,  f=1,...,F_i
%
%  As smoothing parameter \rho increases toward its upper limit of 1,
%  the second roughness penalty term in J is increasing emphasized
%  and the fit to the data decreasing emphasized, so that x_i is
%  required to approach a solution to its respective equation.
%
%  ------------------------------------------------------------------------
%
%  Arguments:
%
%  data           A cell array of length NVAR. Each cell contains in turn
%                 the x locations and the observations at the x locations for each curve.
%  XbasisCell     A cell array of length NVAR.  Each cell contains in turn
%                 a functional data object or a BASIS object for each curve.
%  coef_ODE       the intial contions for the parameters of the ODE.
%  estimate       An indictor variable denoting if the parameter is to be estimated (1) or
%                 kept fixed (0).
%  F              A cell array of length NVAR.  Each cell contains in turn
%                 a functional data object for the forcing function for term in the ODE.
%  plots          A indictor variable denoting if the plots are required (1) or
%                 not (0).
%
%  ------------------------------------------------------------------------
%
%  Output objects
%  Res     ... A matrix containing the estimated parameters, their standard deivation,
%          the mean squared error, the integrated squared error of the
%          differential equation, the degrees of freedom and the
%          generalized cross validation.

if nargin <  5,  FCell{1} = [];  end
if nargin <  6,  plots = 0;  end

%Variables in the dynamical system
neqns   = length(data);

%Set up the rho values
gamvec  = [-4:2,3:0.5:10];
rhoMat  = exp(gamvec)./(1+exp(gamvec));

%Set up the basis functions
basismat    = cell(neqns,1);
for i = 1:neqns
    basismat{i} = eval_basis(data{i}(:,1),XbasisCell{i});
end

%Data
yCell           = cell(neqns,1);
for i = 1:neqns
    yStruct.argvals = data{i}(:,1);
    yStruct.y       = data{i}(:,2);
    yCell{i}        = yStruct;
end

%Model Cell and CoefCell
[modelCell,coefCell] = create_modelCell(coef_ODE,FCell,...
    XbasisCell,estimate);

%Intilize BtensorCell
rho_int     = repmat(0.5,neqns,1);
[MSE, DSSE, D2SSE] = Data2LD(yCell, XbasisCell, modelCell, coefCell, rho_int);

%Set up Dynamical System
%# find empty cells
ntheta        = size(DSSE,1);

nP         = length(rhoMat);
thetastore = zeros(ntheta,nP);
std_theta  = zeros(ntheta,nP);
dfstore    = zeros(nP,neqns);
gcvstore   = zeros(nP,neqns);
MSEstore   = zeros(nP,neqns);
ISEstore   = zeros(nP,neqns);
fd_Par     = cell(nP,1);
CI         = cell(nP,2);
V          = cell(nP,1);

dbglev   = 1;      %  debugging level
iterlim  = 200;    %  maximum number of iterations
conv     = 1e-8;   %  convergence criterion

coefCell_opti = coefCell;

tic;
for irho = 1:nP
    rhoi = repmat(rhoMat(:,irho),neqns,1);
    theta_opti = Data2LD_Opt(yCell, XbasisCell, modelCell, ...
        coefCell_opti, rhoi, ...
        conv, iterlim, dbglev);
    coefCell_opti = BAwtvec2cell(theta_opti, coefCell_opti);
    [MSE, ~, ~, XfdParCell, df, gcv, ISE, Var_theta, Dcoef] = ...
        Data2LD(yCell, XbasisCell, modelCell, coefCell_opti, rhoi);
    CI{irho}           = {Var_theta, Dcoef};
    fd_Par{irho}       = XfdParCell;
    thetastore(:,irho) = theta_opti;
    V{irho}            = Var_theta;
    std_theta(:,irho)  = sqrt(diag(Var_theta));
    MSEstore(irho,:)     = MSE';
    ISEstore(irho,:)     = ISE';
    dfstore(irho,:)      = df;
    gcvstore(irho,:)     = gcv;
end
toc

%Optimal rho
if(size(thetastore,1)>1)
    ind_rho = find(rhoMat>0.99);
    ind_rho(end) = [];
    con = sum(abs(diff(thetastore'))');
    ind = find(con(ind_rho)==min(con(ind_rho)));
    ind = min(sum(rhoMat<0.99)+ind(1)+1,size(rhoMat,2));
else
    ind_rho = find(rhoMat>0.99);
    ind_rho(end) = [];
    con = abs(diff(thetastore'))';
    ind = find(con(ind_rho)==min(con(ind_rho)));
    ind = sum(rhoMat<0.99)+ind(1);
end




%Variance of the fitted curve
m2 = 0;
Results_cell = cell(4,neqns);
Results_cell{2,1} = [thetastore(:,ind),std_theta(:,ind),V{ind}];
for i = 1:neqns
m1        = m2 + 1;
m2        = m2 + size(basismat{i,1},2);
tmp       = CI{ind}{2}(m1:m2,:);
Var_FC{i} = diag(basismat{i,1}*tmp*CI{ind}{1}*tmp'*basismat{i,1}');
Results_cell{1,i} = fd_Par{ind}{i};
Results_cell{3,i} = [MSEstore(ind,i),ISEstore(ind,i),dfstore(ind,i),gcvstore(ind,i)];
fitted_curve = eval_fd(yStruct.argvals, getfd(fd_Par{ind}{i}));
LCI = fitted_curve-tinv(0.975,length(yStruct.argvals)-dfstore(ind,i)).*sqrt(Var_FC{i});
HCI = fitted_curve+tinv(0.975,length(yStruct.argvals)-dfstore(ind,i)).*sqrt(Var_FC{i});
LPI = fitted_curve-tinv(0.975,length(yStruct.argvals)-dfstore(ind,i)).*sqrt(Var_FC{i}+MSEstore(ind,i));
HPI = fitted_curve+tinv(0.975,length(yStruct.argvals)-dfstore(ind,i)).*sqrt(Var_FC{i}+MSEstore(ind,i));
Results_cell{4,i} = [fitted_curve,LCI,HCI,LPI,HPI];
Results_cell{5,i} = {modelCell,coefCell};
end

%Results for all rho
if(neqns==1)
Restab          = num2cell(round([rhoMat',thetastore',...
    std_theta', MSEstore,ISEstore,dfstore,gcvstore],3));
else
Restab          = num2cell(round([rhoMat',thetastore',...
    std_theta', sum(MSEstore')',sum(ISEstore')',dfstore(:,1),sum(gcvstore')'],3));
end
names(1)        = {'Rho  '};
names(2:(ntheta+1)) = cellstr(arrayfun(@(num) ['Theta',num2str(num)], ...
    1:ntheta , 'UniformOutput', false));
names((ntheta+2):(2*ntheta+1)) = cellstr(arrayfun(@(num) ['sdTheta',num2str(num)], ...
    1:ntheta , 'UniformOutput', false));
names((2.*ntheta+2):(2.*ntheta+5))={'MSE  ', 'ISE  ', 'Df   ', 'GCV  '};
Res     = cell2table(Restab,'VariableNames',names');

%Optimal rho
disp('Optimal rho values:')
disp(Res(ind,:));

disp('All rho values:')
if(size(Res,2)<=8)
    disp(Res)
else
    disp(Res(:,1:8))
    disp(Res(:,9:end))
end

if(plots==1)
    
    %Plot the data, fitted curve and confidence and prediction interval
    figure()
    hp = patch([yStruct.argvals; yStruct.argvals(end:-1:1);...
        yStruct.argvals(1)],[LPI; HPI(end:-1:1); LPI(1)],...
        [0.85 0.85 0.85]);
    hold on;
    plot(yStruct.argvals,[LCI,HCI],'k--')
    plot(yStruct.argvals,fitted_curve,'k-')
    phdl=plot(yStruct.argvals, yStruct.y, 'ko');
    xlabel('x'), ylabel('y')
    legend('Prediction Interval','Confidence Interval','Confidence Interval','Fitted Curve','Data')
    
    %Plot the Residuals
    figure()
    plot(yStruct.argvals,yStruct.y-fitted_curve,'bo')
    xlabel('x'), ylabel('y'), title('Residuals')
    
    %Plot the parameters as a function of rho
    figure()
    for i=1:ntheta
        subplot(ntheta,1,i)
        indexOfInterest = (rhoMat > 0.95);
        plot(rhoMat(indexOfInterest),thetastore(i,indexOfInterest),'-*k')
        hold on;
        errorbar(rhoMat(indexOfInterest),thetastore(i,indexOfInterest),1.96.*std_theta(i,indexOfInterest),'ko')
        xlim([0.9815,1])
        ylabel(['$\theta$',num2str(i)],'interpreter','latex')
        xlabel('$\rho$','interpreter','latex')
    end
    
    tilefigs
end

end