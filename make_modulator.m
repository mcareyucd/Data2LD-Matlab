function modStr = make_modulator(funobj, parvec, estimate, name)
if nargin < 4, name = 'beta'; end
if nargin < 3, estimate = 1;  end
if ~strcmp(class(funobj),'function_handle') && ...
   ~isa_basis(funobj) && ...
   ~isa_fd(funobj)    && ...
   ~isa_fdPar(funobj)
   error(['FUNOBJ is not a function handle, basis, fd object or '...
          ' a basis object.'])
end
modStr.fun      = funobj;
modStr.parvec   = parvec;
modStr.estimate = estimate;
modStr.name     = name;