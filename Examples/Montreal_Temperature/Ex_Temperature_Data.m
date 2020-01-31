%% Demonstration of Data2LD: Montreal Temperature Data
%
% Much of the solar radiation reaching the earth is absorbed by the 
% land and sea and then gradually released back into the atmosphere.
%
% Energy balance models, attributable to (Budyko, M.I. 1969), model the rate of 
% change in the average temperature by balancing the absorbed solar radiation 
% and the emitted terrestrial radiation.
%
% The absorbed solar radiation is approximated by 
% $\alpha\cos(2\pi(t+192)/365),$  where the translation $(t+192)$ places 
% the minimum of the cosine at winter solstice. 
%
% The emitted terrestrial radiation can be denoted by $-\beta(t) x(t),$ where $\beta(t)$ is a 
% time-varying parameter. 
%
% We model the 5-day averages of the temperature in Montreal with the 
% first order differential equation:
% 
% $$\frac{\textrm{d}y}{\textrm{d}t} = \beta(t) y  +
% \alpha \cos\left(\frac{2\pi(t+192)}{365}\right)$$
%
% Could we use this first-order ordinary-differential equation to decribe the  
% temperature in Montreal?
%
% This page provides a detailed description of the MATLAB calculations
% necessary to run Data2LD code for estimating the parameters and the
% solution of a linear differential equation to answer this question.

%%  Load the data
%Set up
clear
clear functions
close all

%  load the data files
fid    = fopen('dailtemp.dat','rt');
tempav = fscanf(fid,'%f');
tempav = reshape(tempav, [365,35]);

daytime   = (1:365)'-0.5;
dayrange  = [0,365];
dayperiod = 365;

%  set up centers of 73 5-day blocks
daytime73 = linspace(2.5,362.5,73)';

%  define 8-character names for stations
place = [ ...
    'Arvida     '; 'Bagottville'; 'Calgary    '; 'Charlottown'; ...
    'Churchill  '; 'Dawson     '; 'Edmonton   '; 'Fredericton'; ...
    'Halifax    '; 'Inuvik     '; 'Iqaluit    '; 'Kamloops   '; ...
    'London     '; 'Montreal   '; 'Ottawa     '; 'Pr. Albert '; ...
    'Pr. George '; 'Pr. Rupert '; 'Quebec     '; 'Regina     '; ...
    'Resolute   '; 'Scheffervll'; 'Sherbrooke '; 'St. Johns  '; ...
    'Sydney     '; 'The Pas    '; 'Thunder Bay'; 'Toronto    '; ...
    'Uranium Cty'; 'Vancouver  '; 'Victoria   '; 'Whitehorse '; ...
    'Winnipeg   '; 'Yarmouth   '; 'Yellowknife'];

%  set up indices that order the stations from east to west to north
geogindex = ...
    [24,  9, 25, 34,  4,  8, 22,  1,  2, 19, 23, 14, 15, 28, 13, ...
     27, 33, 26,  5, 20, 16, 29,  7,  3, 12, 30, 31, 17, 18, 32, ...
     6, 35, 11, 10, 21];
place  = place(geogindex,:);
tempav = tempav(:,geogindex);

%Specify the month
monthletter = ['(Jan)'; '(Feb)'; '(Mar)'; '(Apr)'; '(May)'; '(Jun)'; '(Jul)'; '(Aug)'; '(Sep)'; '(Oct)'; '(Nov)'; '(Dec)'];

%For Montreal caluate 5 day average temp
st=12; %  Montreal
tempav73 = zeros(73,1);
m2 = 0;
for i=1:73
    m1 = m2 + 1;
    m2 = m2 + 5;
    tempavi = mean(tempav(m1:m2,st));
    tempav73(i,:) = tempavi;
end

% indices for the winter-centered data
winterind73      = [ 37:73, 1:36];
wintermonthindex = [7:12,1:6];
tempav73         = tempav73(winterind73,:);
monthletter      = monthletter(wintermonthindex,:);

% Create a cell containing the data
yCell = cell(1);
yCell{1}(:,1)  = daytime73;
yCell{1}(:,2)  = tempav73;

% Plot the data
figure(1)
phdl=plot(yCell{1}(:,1), yCell{1}(:,2), 'ko');
set(phdl, 'LineWidth', 2)
xlabel('Time (days)')
ylabel('Temperature')

%%
% Each data point in Figure (1) illustrates the 5-day average of the daily temperature in
% Montreal averaged over 1960 to 1994. 

%% The first order linear differential equation
%
% The first order linear differential equation, as used in this example is
%
% $$Dx(t) =  \beta(t) x(t) + \alpha \cos\left(\frac{2\pi(t+192)}{365}$$
%
% The parameter $\beta(t)$ represents the emitted terrestrial
% radiation. We represent this with a linear combination of fourier basis functions
% $﻿\beta(t)  \approx \sum_j   b_{j} \psi_{j}(t) $ 
%
% Our aim is to estimate the temperature $x(t)$
% and the parameter $\theta= [b_{1},\ldots,b_{J},\alpha]$
% from the data in Figure (1).

%% Set up basis functions for X
norder    = 5;
knots     = 0:5:365;
nbreaks   = length(knots);
nbasis    = norder + nbreaks - 2;
rng       = [0,365];
basis     = create_bspline_basis(rng,nbasis,norder,knots);

%%
% For the expansion of $x(t)$ we used order four B-spline functions which by their nature have cubic 1st derivatives.
% We used 74 equally spaced knots between the first day of the year and last day of the year. 

%Create a cell for the basis for each variable
XbasisCell = cell(1);
XbasisCell{1} = basis;

%Plot the basis
figure(2)
plot(basis)

%% The forcing function

%%
% The explanatory/forcing function is $\cos\left(\frac{2\pi(t+192)}{365}\right)$
uvec      = -cos((2*pi/365)*(daytime73+192));
Ucosbasis = create_fourier_basis(rng, 3);
Ucosfd    = smooth_basis(daytime73, uvec, Ucosbasis);

%Plot the explanatory/forcing function
figure(2)
plot(Ucosfd)

% The explanatory/forcing function is $1$
basisobjC = create_constant_basis(rng);
confd     = fd(ones(getnbasis(basisobjC),1),basisobjC);

%Create a cell for the explanatory/forcing function for each variable
FCell = cell(1,2);
FCell{1,1} = confd;
FCell{1,2} = Ucosfd;

%% The basis function for beta

%%
% For the expansion of $\beta(t)$ we used order five fourier basis functions.
nWbasis   = 7;
Wbasisobj = create_fourier_basis(rng, nWbasis);
forfd     = fd(zeros(getnbasis(Wbasisobj),1),Wbasisobj);

%Plot the basis
figure(2)
plot(Wbasisobj)
 
%% Run Data2LD

%Intial parameter estimates
intial_theta = {forfd,1,1};

%Which parameters would we like to estimate. In this case all parameters.
estimate = [ones(getnbasis(Wbasisobj),1)',1,1];

%Run Data2LD
[Res,Results_cell] = D2LD_ODE(yCell,XbasisCell,intial_theta,...
    estimate,FCell);

%%
% The GCV indicates a model with a lower $\rho$ provides a model with
% better predictive power. This indicated that the differential equation
% does not describe the behaviour of the data well and perhaps a more
% complex equation would provide a better description of the data.

%% Plot the fitted curve with its confidence and prediction intervals
%
% This Figure illustrates the Montreal temperature indicated by the circles.
% The fitted curve produced by Data2LD with $\hat{\rho} = 0.99$ (solid
% line), the approximated 95% pointwise confidence interval for the curve
% (dashed line) and the approximated 95% pointwise prediction interval for
% the curve (grey region).
%
% It seems the differential equation is too simply to capture the behaviour
% of the cancer cases maybe we should try letting the parameter $\beta$
% vary over age. This implies that the rate of growth/decline of cancer
% incidences changes with advancing age. See () for the analysis with a 
% time-varying parameter.

figure(4)
hp = patch([yCell{1}(:,1); yCell{1}(end:-1:1,1);yCell{1}(1)],...
    [Results_cell{4,1}(:,4); Results_cell{4,1}(end:-1:1,5); Results_cell{4,1}(1)],[0.85 0.85 0.85]);
hold on;
plot(yCell{1}(:,1),[Results_cell{4,1}(:,2),Results_cell{4,1}(:,3)],'k--')
plot(yCell{1}(:,1),Results_cell{4,1}(:,1),'k-')
plot(yCell{1}(:,1), yCell{1}(:,2), 'ko');
xlim([min(yCell{1}(:,1))-0.1,max(yCell{1}(:,1))+0.1])
ylim([min(Results_cell{4,1}(:,4))-0.1,max(Results_cell{4,1}(:,5))+0.1])
xlabel('Time (days)')
ylabel('X-direction')
legend('Prediction Interval','Confidence Interval','Confidence Interval','Fitted Curve','Data','Location','northwest')

%% Plot the estimate of beta with its confidence intervals

Wbasismat = eval_basis(daytime73,Wbasisobj);
indW      = 3:(2+getnbasis(Wbasisobj));
StdErrW   = sqrt(diag(Wbasismat*Results_cell{2,1}(1:getnbasis(Wbasisobj),indW)*Wbasismat'));
beta_fd   = fd(Results_cell{2,1}(1:getnbasis(Wbasisobj),1),Wbasisobj);
beta_0    = eval_fd(daytime73, beta_fd);

mtnind = linspace(yCell{1}(1,1),yCell{1}(end,1),12);

h = figure(1);
plot(daytime73,beta_0,'-k')
hold on;
plot(daytime73, beta_0+1.96.*StdErrW  ,'--k')
plot(daytime73, beta_0-1.96.*StdErrW  ,'--k')
xlim([min(daytime73),max(daytime73)])
line([[5,mtnind(2:end)]', [5,mtnind(2:end)]'], [min(beta_0-1.96.*StdErrW)-0.002 min(beta_0-1.96.*StdErrW)]-0.0015, 'Color', 'k')
text([5,mtnind(2:end)], (min(beta_0-1.96.*StdErrW)-0.0015) * ones(size(mtnind)), monthletter, 'Rotation', 90, 'VerticalAlignment', 'middle')
ylim([min(beta_0-1.96.*StdErrW)-0.002,max(beta_0+1.96.*StdErrW)+0.001])
xlabel('Days (Month)')
ylabel('$\hat{\beta}(t)$','Interpreter','latex')
legend('Estimated Function','95% pointwise confidence interval','Location','best')

