function [q,dif,iter] = ergodicdist(Q,tol,maxit,alg)
% ErgodicDist - Ergodic distribution of discrete Markov Chains
%
% q=ergodicdist(Q,tol,maxiter,alg)
%
% Input: Q -   right stochastic matrix (ROWS sum to one): the prob.
%              of reaching state j from state i is P(i,j). ErgodicDist runs
%              faster if Q is sparse.
%        tol   - Numerical tolerance (default = 10^(-8) )
%        maxit - Maximum no. of iterations (only if alg = 1).
%        alg   - choice of method: 1: iterative, 2: direct. (default = 1)
%
% Output: q - (nx1) vector containing the ergodic distribution
%
% Author: Marco Maffezzoli. Ver. 1.0.1, 11/2012.
% Modified by Alessandro Di Nola on 07/2022.
% aledinola13@gmail.com

if isempty(tol)
    tol = 1e-8;
end
if isempty(maxit)
    maxit = 10000;
end
if isempty(alg)
    alg = 1;
end
if ~issparse(Q)
    warning('Stochastic matrix Q is not sparse: ergodicdist runs slow')
end

h=size(Q,1);
switch alg
    case 1
        q=zeros(1,h);
        q(1,1)=1;
        %q = ones(1,h)/h;
        dif=1;
        iter = 0;
        while dif>tol && iter<=maxit
            iter = iter+1;
            z=q*Q;
            dif=max(abs(z-q));
            q=z;
        end
        q=q';
    case 2
        q=(speye(h)-Q'+ones(h))\ones(h,1);
        iter = 1;
        dif  = 0;
    otherwise
        error('Please set alg = 1 or 2!')
end

end %end function