%% Demonstration of Data2LD: Head Impact Analysis
%
% 
% This page provides a detailed description of the MATLAB calculations
% necessary to run Data2LD code for estimating the parameters and the 
% solution of linear differential equations. 

%% Set-up
clear 
clear functions

%%  Load the data
load motorcycledata.txt

motot = motorcycledata(:,2);  %  time in milliseconds
motoy = motorcycledata(:,3);  %  deformation in 0.1 millimeters

% Adjust the data for baseline and plot

impact  = 14;  %  impact time
baseind = motot < impact;
basey   = mean(motoy(baseind));

% Remove baseline, change time, and convert to centimeters
motoy   = (basey - motoy)/100.0;  
n       = length(motoy);

% Create a cell containing the data
yCell = cell(1);
yCell{1}(:,1)  = motot;
yCell{1}(:,2)  = motoy;

% Plot the data
figure(1)
phdl=plot(yCell{1}(:,1), yCell{1}(:,2), 'ko', [0,60], [0,0], 'k:');
set(phdl, 'LineWidth', 2)
axis([min(yCell{1}(:,1))-0.01,max(yCell{1}(:,1))+0.01,-1.0,1.5])
xlabel('\fontsize{16} Time (milliseconds)')
ylabel('\fontsize{16} Acceleration (cm/msec^2)')

%% 
% Figure (1) illustrates $133$ observations of head acceleration 
% measured 14 milliseconds before and 42.6 milliseconds 
% after a blow to the cranium of a cadaver. The dashed line represents 
% the impulse function which denotes the blow to the cranium that lasted 1 
% millisecond. The experiment, a simulated motor-cycle crash, is described
% in detail by Schmidt et al (1981). 

%% The second order linear differential equation
% 
% The second order linear differential equation, as used in this example is
% 
% $$D^2x(t) =  \beta_{0} x(t) + \beta_{1}  Dx(t) + \alpha  u(t)$$
%
% The three parameters $\beta_{0}$, $\beta_{1} $ and $\alpha$ in 
% the linear ODE convey the period of the oscillation, 
% the change in its amplitude, as $t \rightarrow \infty$ the oscillations 
% decay exponentially to zero, and the size of the impact from the unit 
% impulse respectively. 
%
% Our aim is to estimate the acceleration $x(t)$ 
% and the parameters $\theta=[\beta_{0},\beta_{1},\alpha]$ 
% from the data in Figure (1).


%% Set up basis functions for X
motorng   = [yCell{1}(1,1),yCell{1}(end,1)];
knots     = [motorng(1),14,14,14,15,15,linspace(15,motorng(end),11)];
norder    = 6;
nbasis    = length(knots)+ (norder - 2);
Xbasisobj = create_bspline_basis(motorng,nbasis,norder,knots);

%Create a cell for the basis for each variable
XbasisCell = cell(1);
XbasisCell{1} = Xbasisobj;

%Plot the basis
figure()
plot(Xbasisobj)

%% 
% For the expansion of $x(t)$ we used order five B-spline functions, which by their nature have 
% discontinuous third derivatives if all knots are singletons.  To achieve curvature discontinuity at the 
% impact point and at that point plus one, we placed three knots at these locations.  We put no knots 
% between the first observation and the impact point, where the data indicate a flat trajectory, and eleven 
% equally spaced knots between the impact point plus one and the final time of observation.

%% Create a functional data object for the external/forcing function

%The known external function is a unit impulse 
Ubasis = create_bspline_basis(motorng, 3, 1, [motorng(1),14,15,motorng(2)]);
Ufd = fd([0;1;0],Ubasis);

%Plot the external function
figure()
plot(Ufd)

%Create a cell for the external functions
EFCell = cell(1);
EFCell{1} = Ufd;

%%
% 
% Three order one B-splines over the knots [0, 14, 15, 56] with 
% coefficient vector [0,1,0] represented the unit pulse function $u(t)$
% 

%% Run Data2LD

%Intial parameter estimates
intial_theta = {0.01,0.01,0.01};

%Which parameters would we like to estimate. In this case all parameters.
estimate = [1,1,1];

%Run Data2LD
[Res,Results_cell] = D2LD_ODE(yCell,XbasisCell,intial_theta,...
    estimate,EFCell);

%Confidence intervals for theta
confidence_intervals(Results_cell);

%% Plot the fitted curve, the confidence and prediction intervals

figure()
hp = patch([yCell{1}(:,1); yCell{1}(end:-1:1,1);yCell{1}(1)],...
    [Results_cell{4}(:,4); Results_cell{4}(end:-1:1,5); Results_cell{4}(1)],[0.85 0.85 0.85]);
hold on;
plot(yCell{1}(:,1),[Results_cell{4}(:,2),Results_cell{4}(:,3)],'k--')
plot(yCell{1}(:,1),Results_cell{4}(:,1),'k-')
plot(yCell{1}(:,1), yCell{1}(:,2), 'ko');
xlim([min(yCell{1}(:,1))-0.1,max(yCell{1}(:,1))+0.1])
ylim([min(Results_cell{4}(:,4))-0.1,max(Results_cell{4}(:,5))+0.1])
xlabel('\fontsize{16} Time (milliseconds)')
ylabel('\fontsize{16} Acceleration (cm/msec^2)')
legend('Prediction Interval','Confidence Interval','Confidence Interval','Fitted Curve','Data')
set(findall(gca, 'Type','text'), 'FontSize', 20)


%%
% 
%  This Figure illustrates the accelerometer readings of the brain tissue before and after a blow to the cranium are indicated by the circles. 
%  The fitted curve produced by Data2LD with $\hat{\rho} = 0.99$ (solid
%  line), the approximated 95\% pointwise confidence interval for the curve
%  (dashed line) and the approximated 95\% pointwise prediction interval for the curve (grey region)
% 

%% Plot the estimated parameters of the ODE and their approximated CI with respect to rho
Res_Mat  = table2array(Res);
for i=1:3
subplot(3,1,i)
indexOfInterest = (Res_Mat(:,1) > 0.95);
plot(Res_Mat(indexOfInterest,1),Res_Mat(indexOfInterest,i+1),'-*k')
hold on;
errorbar(Res_Mat(indexOfInterest,1),Res_Mat(indexOfInterest,i+1),...
    t_v.*Res_Mat(indexOfInterest,i+4),'ko')
xlabel('$\rho$','interpreter','latex')
xlim([min(Res_Mat(indexOfInterest,1)),max(Res_Mat(indexOfInterest,1))])
end
subplot(3,1,1)
ylabel('Stiffness $\hat{\beta_{0}}$','interpreter','latex')
subplot(3,1,2)
ylabel('Damping $\hat{\beta_{1}}$','interpreter','latex')
subplot(3,1,3)
ylabel('External $\hat{\alpha}$','interpreter','latex')

%%
% 
%  The values of the three parameters and their approximated 95\% confidence intervals for
%  the head impact data over values of $\rho$ converging to 1.0. Note the reduction in the
%  width of the approximated 95\% confidence intervals as $\rho$  increase and the overlap
%  of the final estimates of the ODE parameters, which illustrates that the parameter values
%  have stabilised for high values of $\rho$.

%% Check the residuals 

%% Check the normality of the residuals
res = yCell{1}(:,2)-Results_cell{4}(:,1);
figure()
subplot(2,1,1)
h = normplot(res);
h(1).LineWidth = 2;
h(2).LineWidth = 2;
h(3).LineWidth = 2;
subplot(2,1,2)
hist(res)

% Our qqplot looks pretty good and indicates a normal distribution, as it’s
% generally a straight line.

%% Check the fitted vrs residuals plot

figure()
plot(Results_cell{4}(:,1),res,'ko')

% Should see random scatter in a constant band around zero that is no
% patterns, and no funnels.