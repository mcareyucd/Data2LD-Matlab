%% Demonstration of Data2LD: ﻿Incidence of cancer vrs age for American females in 2001
%
% Most people know that their chances of getting cancer increase as they age. 
% In fact, by looking at data complied by the National Cancer Institute you 
% can readily see that incidence of cancer increases dramatically between 
% the ages of 35 and 80.
%
% A first-order ordinary differential equation 
%
% $$\frac{\textrm{d}y}{\textrm{d}x} = \beta y$$
%
% has the solution
%
% $$y(x) = C \exp^{\beta x},$$ 
%
% where $C$ is any constant. This equation has an infinitely many 
% solutions depending on the coefficient $C$.
%
% However, this solution exhibits:
%
% * exponential decay when $\beta<0$ and 
% * exponential growth when $\beta>0$
%
% Could we use this simple first-order ordinary-differential equation to decribe the  
% increase in the cancer incidences with advancing age?
%
% This page provides a detailed description of the MATLAB calculations
% necessary to run Data2LD code for estimating the parameters and the
% solution of a linear differential equation to answer this question.

%% Set-up
clear
clear functions
close all

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
% $$Dx(t) =  \beta_{0} x(t)$$
%
% The parameter $\beta_{0}$ in the linear ODE conveys the exponentail growth 
% (if positive) or decay (if negative) in the 
% cancer cases per 100,000 population as the individuals increase in age.
%
% Our aim is to estimate the cancer cases  $x(t)$
% and the parameter $\theta=\beta_{0}$
% from the data in Figure (1).

%% Set up basis functions for X
rng       = [yCell{1}(1,1),yCell{1}(end,1)];
knots     = linspace(rng(1),rng(end),17);
norder    = 4;
nbasis    = length(knots)+ (norder - 2);
Xbasisobj = create_bspline_basis(rng,nbasis,norder,knots);

%Create a cell for the basis for each variable
XbasisCell = cell(1);
XbasisCell{1} = Xbasisobj;

%Plot the basis
figure(2)
plot(Xbasisobj)

%%
% For the expansion of $x(t)$ we used order four B-spline functions which by their nature have quadratic 1st derivatives.
% We used 17 equally spaced knots between the first age and last age. 

%% Run Data2LD

%Intial parameter estimates
intial_theta = {0.1};

%Which parameters would we like to estimate. In this case all parameters.
estimate = 1;

%Run Data2LD
[Res,Results_cell] = D2LD_ODE(yCell,XbasisCell,intial_theta,...
    estimate);

% Confidence intervals for theta
confidence_intervals(Results_cell);

%%
% The GCV indicates a model with a lower $\rho$ provides a model with
% better predictive power. This indicates that the differential equation
% does not describe the behaviour of the data well and perhaps a more
% complex equation would provide a better description of the data.

%% Plot the fitted curve with its confidence and prediction intervals
%
% This Figure illustrates the cancer cases per 100,000 population
% indicated by the circles.
% The fitted curve produced by Data2LD with $\hat{\rho} = 0.99$ (solid
% line), the approximated 95% pointwise confidence interval for the curve
% (dashed line) and the approximated 95% pointwise prediction interval for
% the curve (grey region).
%
% It seems the differential equation is too simply to capture the behaviour
% of the cancer cases maybe we should try letting the parameter $\beta$
% vary over age. This implies that the rate of growth/decline of cancer
% incidences changes with advancing age. See Cancer Cases Time Varying for the analysis with a 
% time-varying parameter.

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

%% Plot the estimated parameters of the ODE and their approximated CI with respect to rho
%
% The values of the parameter and its approximated 95% confidence interval
% over values of $\rho$ converging to 1.0. Note the reduction in the
% width of the approximated 95% confidence intervals as $\rho$  increase and the overlap
% of the final estimates of the ODE parameter, which illustrates that the parameter value
% has stabilised for high values of $\rho$.
%
% The optimal estimate of the parameter $\beta$ is $0.037$ with a 95%
% confidence interval of 
Res_Mat         = table2array(Res);
df              = Results_cell{3,1}(3);
t_v             = tinv(0.975,length(yCell{1}(:,1))-df);
indexOfInterest = (Res_Mat(:,1) > 0.95);

figure(5)
plot(Res_Mat(indexOfInterest,1),Res_Mat(indexOfInterest,2),'-*k')
hold on;
errorbar(Res_Mat(indexOfInterest,1),Res_Mat(indexOfInterest,2),...
    t_v.*Res_Mat(indexOfInterest,5),'ko')
xlabel('$\rho$','interpreter','latex')
ylabel('$\hat{\beta_{0}}$','interpreter','latex')
xlim([min(Res_Mat(indexOfInterest,1)),max(Res_Mat(indexOfInterest,1))])

%% Check fit to the ODE

Dfit   = eval_fd(yCell{1}(:,1), getfd(Results_cell{1}),1);
RHS    = Results_cell{2}(1).*Results_cell{4}(:,1);

figure()
plot(yCell{1}(:,1),Dfit,'*k')
hold on;
plot(yCell{1}(:,1),RHS,'k--')
legend('Left-Hand side of the ODE','Right-Hand side of the ODE','Location','best')
xlabel('Age')
ylabel('Velocity of the cancer cases per 100,000 population')

%The difference between the LHS and RHS of the ODE is negligible indicating
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
