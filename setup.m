%Add all subfolder to the path
folder=fileparts(which('D2LD_ODE.m'));
addpath(genpath(folder));

% remember the original working directory
pwdir = pwd;

% add FDAM package to the path
wherefile=fileparts(which('fdaM.zip'));
funname = [wherefile,'/fdaM.zip'];
unzip(funname)
funname = [wherefile,'/fdaM'];
addpath(funname)

%Compile the C code for innerloop
wherefile=fileparts(which('inner_loop.c'));
funname = [wherefile,'/inner_loop.c'];
[mexdir, mexname] = fileparts(funname);

try
% try to compile the mex file on the fly
warning('trying to compile MEX file from %s', funname);
cd(mexdir);
mex(funname);
cd(pwdir);
success = true;

catch
% compilation failed
disp(lasterr);
error('could not locate MEX file for %s', mexname);
cd(pwdir);
success = false;
end
