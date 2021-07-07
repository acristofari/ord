% -------------------------------------------------------------------------
%
% This file is part of ORD, which is a derivative-free solver for
% optimization problems of the following form:
%
%                                 min f(x)
%                            s.t. x in conv{a_1,...,a_m}
%
% where f(x) is a (black-box) continuously differentiable function and
% conv{a_1,...,a_m} is the convex hull of some given vectors a_1,...,a_m,
% called atoms.
%
% -------------------------------------------------------------------------
%
% Reference paper:
%
% A. Cristofari, F. Rinaldi (2021). A Derivative-Free Method for Structured
% Optimization Problems. SIAM Journal on Optimization, 31(2), 1079-1107.
%
% -------------------------------------------------------------------------
%
% Authors:
% Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
% Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
%
% Last update of this file:
% July 7th, 2021
%
% Licensing:
% This file is part of ORD.
% ORD is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% ORD is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with ORD. If not, see <http://www.gnu.org/licenses/>.
%
% Copyright 2021 Andrea Cristofari, Francesco Rinaldi.
%
% -------------------------------------------------------------------------

clear all, clc;
rng(1);

% In this file, we first show how to call ORD to solve a user-defined problem.
% Then, we also show how to call DF-SIMPLEX as standalone.

% (1) Define an objective function
%--------------------------------------------------------------------------
% Extended Rosenbrock function
n = 50; % dimension
c = 1e1; % parameter of the objective function
obj = @(x) sum(c*(x(2:2:n)-x(1:2:n-1).^2).^2+(1e0-x(1:2:n-1)).^2); % objective function
%--------------------------------------------------------------------------

% (2) Define an n-by-m matrix of atoms, where each column is an n-dimensional atom
%--------------------------------------------------------------------------
% atoms randomly generated
m = 20*n; % number of atoms
A = 1e1*rand(n,m); % matrix of atoms
%--------------------------------------------------------------------------

% (3) Choose an atom to use as starting point
%--------------------------------------------------------------------------
i0 = randi(m); % index of the atom to use as starting point
%--------------------------------------------------------------------------

% (4) call ORD
[x_ord,y_ord,f_ord,n_f_ord,it_ord,t_elap_ord,flag_ord] = ORD(obj,A,i0);

%--------------------------------------------------------------------------
% *** EXAMPLE OF HOW TO CHANGE ORD PARAMETERS ***
% (see the file 'syntax_ord.txt' to know which parameters can be changed
% and their default values)
%
% Instead of calling ORD by the above instruction, do the following:
%
% - create a structure having as field names the names of the parameters
%   to be changed and assign them new values, e.g.,
%
%     opts.verbosity = false;
%
% - pass the structure to ORD as fourth input argument, e.g.,
%
%     [x,y,f,n_f,it,t_elap,flag] = ORD(obj,A,i0,opts);
%--------------------------------------------------------------------------

% write statistics to the screen
fprintf(['\n********************** FINAL RESULTS **********************' ...
         '\nAlgorithm: ORD' ...
         '\nf =  %-.4e'   ...
         '\nobjective function evaluations = %-i' ...
         '\niterations = %-i' ...
         '\ncpu time (s) = %-.2e' ...
         '\nflag = %-i' ...
         '\n***********************************************************\n\n'], ...
         f_ord,n_f_ord,it_ord,t_elap_ord,flag_ord);


% To call DF-SIMPLEX as standalone, follow the same steps as above but call
% DF-SIMPLEX instead of ORD (note that, in DF-SIMPLEX, the starting point
% is passed as a vector of coefficients that express a convex combination
% of the atoms, see the file 'syntax_df_simplex.txt' for further details).
% Namely,

[x_dfsimplex,y_dfsimplex,f_dfsimplex,n_f_dfsimplex,it_dfsimplex,t_elap_dfsimplex,flag_dfsimplex] = ...
    DF_SIMPLEX(obj,A,double(1:m==i0)');
fprintf(['\n********************** FINAL RESULTS **********************' ...
         '\nAlgorithm: DF-SIMPLEX' ...
         '\nf =  %-.4e'   ...
         '\nobjective function evaluations = %-i' ...
         '\niterations = %-i' ...
         '\ncpu time (s) = %-.2e' ...
         '\nflag = %-i' ...
         '\n***********************************************************\n\n'], ...
         f_dfsimplex,n_f_dfsimplex,it_dfsimplex,t_elap_dfsimplex,flag_dfsimplex);

% DF-SIMPLEX parameters can be changed as explained above for ORD, that is,
% by passing as fourth input argument a structure with the parameters to be
% changed (see the file 'syntax_df_simplex.txt' to know which parameters
% can be changed and their default values)