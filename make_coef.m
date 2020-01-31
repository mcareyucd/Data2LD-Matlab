function coefStr = make_coef(funobj, parvec, estimate, coeftype)
if nargin < 4, coeftype = 'beta'; end
if nargin < 3, estimate = 1;  end
if ~strcmp(class(funobj),'function_handle') && ...
   ~isa_basis(funobj) && ...
   ~isa_fd(funobj)    && ...
   ~isa_fdPar(funobj)
   error(['FUNOBJ is not a function handle, basis, fd object or '...
          ' a basis object.'])
end
coefStr.fun      = funobj;
coefStr.parvec   = parvec;
coefStr.estimate = estimate;
coefStr.coeftype = coeftype;