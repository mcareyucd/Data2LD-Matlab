%% Demonstration of Data2LD: Handwriting Data
%
%﻿Generating handwritten characters involves translating a stream of 
% symbols which represent the characters on a cognitive level into a 
% stream of muscle activation commands. Assuming that each of the 46 
% events required an activation command and that the dependence 
% between the acceleration of the pen movements in the $X,$ $Y$ and $Z$ 
% direction and the activation commands can be represented by the 
% unknown external functions $\alpha_{X}(t)$, $\alpha_{Y}(t)$ and $\alpha_{Z}(t)$. 
% The handwritten characters can be represented by a system of three linear ODEs,
%	
% $D^2X(t) =  -\beta_{0} X(t) + \alpha_{X}(t)$
% $D^2Y(t) =  -\beta_{1} Y(t) + \alpha_{Y}(t)$
% $D^2Z(t) =  -\beta_{2} Z(t) + \alpha_{Z}(t)$
%
% where the parameters $\beta_{0}, \beta_{1}$ and $\beta_{2}$ represent the 
% natural frequency of the pen movements along the X, Y and Z dimensions 
% and $\alpha_{X}(t)$, $\alpha_{Y}(t)$ and $\alpha_{Z}(t)$ represent the 
% contribution of the activation commands.
%
% Could we use this first-order dynamical system to decribe the  
% handwritten characters?
%
% This page provides a detailed description of the MATLAB calculations
% necessary to run Data2LD code for estimating the parameters and the
% solution of the dynamical system to answer this question.

%% Set-up
clear 
clear functions

%%  Load the data
load chinastat.txt

% create an array object containing observations in centimetres
XYZ = reshape(chinastat, [2401,50,3])./10;

% set up the time values
n        = 2401;  %  number of observations points at 200 herz
sec      = linspace(0,6,n)';
centisec = sec*100;  %  at better choice for using splines
records  = 1:50;

% set up functions to be smoothed
X = XYZ(:,:,1);
Y = XYZ(:,:,2);
Z = XYZ(:,:,3);

% Create a cell containing the data
record=1;
yCell = cell(3,1);
yCell{1}(:,1)  = centisec;
yCell{1}(:,2)  = XYZ(:,record,1);
yCell{2}(:,1)  = centisec;
yCell{2}(:,2)  = XYZ(:,record,2);
yCell{3}(:,1)  = centisec;
yCell{3}(:,2)  = XYZ(:,record,3);

% plot XY coordinates below a threshold

% there are three types of events:
% --  strokes
% --  cusps
% --  lifts

% Here we position a number at each of 45 equally spaced times
% located within 45 equal-sized intervals

nevent    = 46;
eventtime = linspace(0,600,2*nevent+1);
eventind  = 2:2:2*nevent;
thresh    = 0.27;

% plot each record
figure(1)
% indices of pen positions with Z below thresh
inddn  = find(XYZ(:,record,3) <= thresh);
% X anfd Y coordinates at the 45 equally spaced times
Xeventctr = interp1(centisec,X(:,record),eventtime(eventind));
Yeventctr = interp1(centisec,Y(:,record),eventtime(eventind));
% plot the positions for which Z is below thresh
plot(X(inddn,record), Y(inddn,record), 'k.')
% now plot numbers at equally spaced time points
hold on
for k=1:nevent
    text(Xeventctr(k),Yeventctr(k),['\fontsize{13} ',num2str(k)])
end
hold off
axis([-14,20,-8,7])
xlabel('\fontsize{16} X (cm)')
ylabel('\fontsize{16} Y (cm)')
title(['\fontsize{16} Record ',num2str(record)])

figure(2)
plot3(yCell{1}(:,2), yCell{2}(:,2), yCell{3}(:,2), 'k-', ...
      yCell{1}(:,2), yCell{2}(:,2), yCell{3}(:,2), 'ko')
xlabel('\fontsize{16} X (cm)')
ylabel('\fontsize{16} Y (cm)')
zlabel('\fontsize{16} Z (cm)')

%%
% Each data point in Figure (1) illustrates the X and Y positions from 
% $0$ to $6$ seconds observed at 400 herz while writing ``statistics'' 
% in Mandarin script. Gaps correspond to vertical positions greater than 
% 0.47 cm. The indices of the 46 intervals are positioned at their centres.

%% The first order dynamical system
%
% The first order dynamical system, as used in this example is
%	
% $D^2X(t) =  -\beta_{0} X(t) + \alpha_{X}(t)$
% $D^2Y(t) =  -\beta_{1} Y(t) + \alpha_{Y}(t)$
% $D^2Z(t) =  -\beta_{2} Z(t) + \alpha_{Z}(t)$
%
% The﻿the parameters $\beta_{0}, \beta_{1}$ and $\beta_{2}$ represent 
% the natural frequency of the pen movements along the X, Y and Z 
% dimensions and $\alpha_{X}(t)$, $\alpha_{Y}(t)$ and $\alpha_{Z}(t)$ 
% represent the contribution of the activation commands. Approximate 
% the unknown external functions $\alpha_{X}(t), \alpha_{Y}(t),$ and 
% $\alpha_{Z}(t)$  by basis function expansions, that is, 
% $\alpha_{X}(t) = \sum_{\ell=1}^{L_{1}} a_{\ell 1} \xi_{\ell  1}(t)$, 
% $\alpha_{Y}(t) = \sum_{\ell=1}^{L_{2}} a_{\ell 2} \xi_{\ell  2}(t)$ 
% and $\alpha_{Z}(t) = \sum_{\ell=1}^{L_{3}} a_{\ell 3} \xi_{\ell  3}(t),$ 
% where $\xi_{\ell 1}(t)$, $\xi_{\ell  2}(t)$ and $\xi_{\ell  3}(t)$ are 
% known basis functions and $a_{\ell 1},a_{\ell 2}$ and $a_{\ell 3}$ are 
% the corresponding coefficients.
%
% Our aim is to estimate the yield curve $X(t),Y(t)$ and $Z(t)$
% and the parameter $\theta= [\beta_{0},\beta_{1},\beta_{2},a_{1 1},\ldots,
% a_{L_{1} 1},a_{1 2},\ldots,a_{L_{1} 2},a_{1 3},\ldots,a_{L_{1} 3}]$
% from the data in Figure (1).

%%  First define each of the terms in the equations

%% Set up basis functions for X
rng       = [yCell{1}(1,1),yCell{1}(end,1)];
knots     = linspace(yCell{1}(1,1),yCell{1}(end,1),196);
norder    = 6;
nbasis    = length(knots)+ (norder - 2);
Xbasisobj = create_bspline_basis(rng,nbasis,norder,knots);

%%
% For the expansion of $X,Y,Z$ we used order six B-spline functions which by their nature have quintic 1st derivatives.
% We used 196 equally spaced knots between the first observation and last observation.  32 basis functions per second, 
% 4 per 8 herz cycle. 

%Create a cell for the basis for each variable
XbasisCell = cell(3,1);
XbasisCell{1,1} = Xbasisobj;
XbasisCell{2,1} = Xbasisobj;
XbasisCell{3,1} = Xbasisobj;

%Plot the basis
figure(2)
plot(Xbasisobj)

%% Create a functional data object for the external/forcing function

%Constant basis for the forcing function u
basisobjC = create_constant_basis(rng);
confd     = fd(ones(getnbasis(basisobjC),1),basisobjC);

%B-spline basis for paramterts associated with the forcing function
%$\alpha$
nAorder   = 1;
nAbasis   = 46;
nAknots   = nAbasis + 1;
Aknots    = linspace(rng(1),rng(2),nAknots);
Abasisobj = create_bspline_basis(rng,nAbasis,nAorder,Aknots);

%%
% 
% The coefficients of the external functions are﻿step functions, with knots at 46 equally spaced 
% locations over the domain [$0,600$] for each coordinate
% 
% The external functions can the be represented by constant functions over
% the domain [$0,600$].

%Create a cell for the external functions u
FCell = cell(3,1);
FCell{1,1} = confd;
FCell{2,1} = confd;
FCell{3,1} = confd;

%% Run Data2LD

%Intial parameter estimates
intial_theta = cell(3,3);
%parameters associated with y
intial_theta{1,1} = 0.01;
intial_theta{2,1} = 0.01;
intial_theta{3,1} = 0.01;
%parameters associated with Dy
intial_theta{1,2} = {};
intial_theta{2,2} = {};
intial_theta{3,2} = {};
%parameters associated with non-homogenous terms of each ODE
intial_theta{1,3} = fd(-0.01.*ones(nAbasis,1),Abasisobj);
intial_theta{2,3} = fd(-0.01.*ones(nAbasis,1),Abasisobj);
intial_theta{3,3} = fd(-0.01.*ones(nAbasis,1),Abasisobj);

%Which parameters would we like to estimate. In this case 
% all but the terms associated with Dy parameters.
estimate = [ones(3,1),zeros(3,1),ones(3,nAbasis)];

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

figure(3)
plot(yCell{1}(:,2), yCell{2}(:,2),'ko')
xlabel('\fontsize{16} X (cm)')
ylabel('\fontsize{16} Y (cm)')
hold on;

figure(4)
hp = patch([yCell{1}(:,1); yCell{1}(end:-1:1,1);yCell{1}(1)],...
    [Results_cell{4}(:,4); Results_cell{4}(end:-1:1,5); Results_cell{4}(1)],[0.85 0.85 0.85]);
hold on;
plot(yCell{1}(:,1),[Results_cell{4}(:,2),Results_cell{4}(:,3)],'k--')
plot(yCell{1}(:,1),Results_cell{4}(:,1),'k-')
plot(yCell{1}(:,1), yCell{1}(:,2), 'ko');
xlim([min(yCell{1}(:,1))-0.1,max(yCell{1}(:,1))+0.1])
ylim([min(Results_cell{4}(:,4))-0.1,max(Results_cell{4}(:,5))+0.1])
xlabel('Time (days)')
ylabel('Temperature')
legend('Prediction Interval','Confidence Interval','Confidence Interval','Fitted Curve','Data','Location','northwest')



