%% Demonstration of Data2LD: ﻿Incidence of cancer vrs age for American females in 2001
%
%
% This page provides a detailed description of the MATLAB calculations
% necessary to run Data2LD code for estimating the parameters and the
% solution of a linear differential equation.

%% Set-up
clear
clear functions

%%  Load the data

load('Cancer_D.mat')

% Create a cell containing the data
yCell = cell(1);
yCell{1}(:,1)  = Cancer_D(:,1); %age of the individuals
yCell{1}(:,2)  = Cancer_D(:,2); %cancer cases per 100,000 population

% Plot the data
figure(1)
phdl=plot(yCell{1}(:,1), yCell{1}(:,2), 'ko');
set(phdl, 'LineWidth', 2)
xlabel('Age')
ylabel('Cancer cases per 100,000 population')

%%
% Figure (1) illustrates $19$ observations of cancer cases 
% measured from age 0.5 to age 87. 

%% The first order linear differential equation
%
% The first order linear differential equation, as used in this example is
%
% $$Dx(t) =  \beta(t) x(t)$$
%
% The parameter $\beta(t)$ in the linear ODE conveys the age-varying rate of growth 
% (if positive) or decay (if negative) in the 
% cancer cases per 100,000 population as the individuals increase in age.
%
% We have no prior knowledge of the form of $\beta(t)$ so we will estimate it with
% a simple step function (B-spline basis of order 1) $\beta(t)=sum_{j=1}^{J}b_j\phi_j(t)$. 
%
% Our aim is to estimate the cancer cases  $x(t)$
% and the parameter $\theta=[b_1,\ldots,b_J]$
% from the data in Figure (1).

%% Set up basis functions for X
rng       = [yCell{1}(1,1),yCell{1}(end,1)];
knots     = linspace(rng(1),rng(end),17);
norder    = 4;
nbasis    = length(knots)+ (norder - 2);
Xbasisobj = create_bspline_basis(rng,nbasis,norder,knots);

%Plot the basis
figure(2)
plot(Xbasisobj)

%%
% For the expansion of $x(t)$ we used order four B-spline functions which by their nature have quadratic 1st derivatives.
% We used 17 equally spaced knots between the first age and last age. 

%Create a cell for the basis for each variable
XbasisCell = cell(1);
XbasisCell{1} = Xbasisobj;

%% Create a functional data object for the beta function
% For the expansion of $\beta(t)$ we used order one B-spline functions over 
% the following knot sequence [0.5,12,27,42,57,70,72,74,76,80,87]. 

rng       = [yCell{1}(1,1),yCell{1}(end,1)];
knots     = [0.5,12,40,45,50,55,60,65,70,75,80,85,87];
norder    = 1;
nbasis    = length(knots)+ (norder - 2);
basisobjC = create_bspline_basis(rng, nbasis, norder, knots);
fd_beta   = fd(0.01.*ones(nbasis,1),basisobjC);

%Plot the beta
figure(2)
plot(fd_beta)
hold on;
plot([knots; knots], [0; 1], '--r')


%% Run Data2LD

%Intial parameter estimates
intial_theta = {fd_beta};

%Which parameters would we like to estimate. In this case all parameters.
estimate = ones(1,nbasis);

%Run Data2LD
[Res,Results_cell] = D2LD_ODE(yCell,XbasisCell,intial_theta,...
    estimate);
 

%%
% The GCV indicates a model with a lower $\rho$ provides a model with
% better predictive power. This indicated that the differential equation
% does not describe the behaviour of the data well and perhaps a more
% complex equation would provide a better description of the data.

%% Plot the fitted curve with its confidence and prediction intervals
%
% This Figure illustrates the cancer cases per 100,000 population
% indicated by the circles. The fitted curve produced by Data2LD 
% with $\hat{\rho} = 0.99$ (solid line), the approximated 95% pointwise 
% confidence interval for the curve (dashed line) and the approximated 
% 95% pointwise prediction interval for the curve (grey region).
%

figure(4)
hp = patch([yCell{1}(:,1); yCell{1}(end:-1:1,1);yCell{1}(1)],...
    [Results_cell{4}(:,4); Results_cell{4}(end:-1:1,5); Results_cell{4}(1)],[0.85 0.85 0.85]);
hold on;
plot(yCell{1}(:,1),[Results_cell{4}(:,2),Results_cell{4}(:,3)],'k--')
plot(yCell{1}(:,1),Results_cell{4}(:,1),'k-')
plot(yCell{1}(:,1), yCell{1}(:,2), 'ko');
xlim([min(yCell{1}(:,1))-0.1,max(yCell{1}(:,1))+0.1])
ylim([min(Results_cell{4}(:,4))-0.1,max(Results_cell{4}(:,5))+0.1])
xlabel('Age')
ylabel('Cancer cases per 100,000 population')
legend('Prediction Interval','Confidence Interval','Confidence Interval','Fitted Curve','Data','Location','northwest')


%% Plot the estimated parameters of the ODE and their approximated CI 

% This Figure illustrates the estimated age-varying rate of growth of 
% the cancer cases per 100,000 population produced by Data2LD 
% with $\hat{\rho} = 0.99$ (solid line) and the approximated 95% pointwise 
% confidence interval for the curve (dashed line).
%

Wbasismat = eval_basis(yCell{1}(:,1),basisobjC);
StdErrW   = sqrt(diag(Wbasismat*Results_cell{2}(:,3:end)*Wbasismat'));

beta_fd   = fd(Results_cell{2}(1:getnbasis(basisobjC),1),basisobjC);
beta_0    = eval_fd(yCell{1}(:,1), beta_fd);

figure(5)
stairs(yCell{1}(:,1),beta_0,'k-');
hold on;
%stairs(yCell{1}(:,1), beta_0+1.96.*StdErrW  ,'--k')
%stairs(yCell{1}(:,1), beta_0-1.96.*StdErrW  ,'--k')
%plot([knots; knots], [min(beta_0-1.96.*StdErrW); max(beta_0+1.96.*StdErrW)], '--r')
%xlim([min(yCell{1}(:,1))-0.1,max(yCell{1}(:,1))+0.1])
%ylim([min(beta_0-1.96.*StdErrW)-0.01,max(beta_0+1.96.*StdErrW)+0.01])
xlabel('Age')
ylabel('$\hat{\beta}(t)$','Interpreter','latex')
 
 
 %% Check fit to the ODE
fit    = Results_cell{4}(:,1);
Dfit   = eval_fd(yCell{1}(:,1), getfd(Results_cell{1}),1);
RHS    = beta_0.*fit;

figure()
plot(yCell{1}(:,1),Dfit,'*k')
hold on;
plot(yCell{1}(:,1),RHS,'k--')
legend('Left-Hand side of the ODE','Right-Hand side of the ODE','Location','best')
xlabel('Age')
ylabel('Velocity of the cancer cases per 100,000 population')

%The difference between the LHS and RHS of the ODE is small indicating
%that the ODE is satified.

%% Check the residuals 

%% Check the normality
res = yCell{1}(:,2)-Results_cell{4}(:,1);
figure()
subplot(2,1,1)
h = normplot(res);
h(1).LineWidth = 2;
h(2).LineWidth = 2;
h(3).LineWidth = 2;
subplot(2,1,2)
hist(res)

% Our qqplot looks indicates that the residuals are not normal distribution, as it’s
% not generally a straight line. The histogram does not have the features
% of a normal distribution (its not symmetric and bell-shaped)

%% Check the fitted vrs residuals plot

figure()
plot(Results_cell{4}(:,1),yCell{1}(:,2)-Results_cell{4}(:,1),'ko')

% Should see random scatter in a constant band around zero that is no
% patterns, and no funnels. We can see clear non-linear pattern here
% indicating that the model does not adquetly capture the complexity of the
% data 
