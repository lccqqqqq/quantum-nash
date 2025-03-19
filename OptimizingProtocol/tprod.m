function A = tprod(varargin)
A = varargin{1};
for k = 2:nargin
     A = kron(A, varargin{k});
end
end
