% -------------------------------------------------------------------------
% 
% This is the main file to run ORD and DF-SIMPLEX.
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
% Andrea Cristofari (e-mail: andrea.cristofari@uniroma2.it)
% Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
%
% Last update of this file:
% November 19, 2026
% 
% -------------------------------------------------------------------------

clear all, clc;

rng(1);

% In this file, it is first shown how to call ORD to solve a user-defined problem
% Then, it is shown how to call DF-SIMPLEX as standalone.

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

% (3) Choose an atom to be used as starting point
%--------------------------------------------------------------------------
i0 = randi(m); % index of the atom to be used as starting point
%--------------------------------------------------------------------------

% (4) call ORD
[x_ord,f_ord,ord_info] = ord(obj,A,i0);

%--------------------------------------------------------------------------
% *** EXAMPLE OF HOW TO CHANGE ORD PARAMETERS ***
% (see the file 'usage_ord.txt' to know which parameters can be changed
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
%     [x_ord,f_ord,ord_info] = ord(obj,A,i0,opts);
%--------------------------------------------------------------------------

% write statistics to the screen
fprintf(['\n********************** FINAL RESULTS **********************' ...
         '\nAlgorithm: ORD' ...
         '\nf =  %-.4e'   ...
         '\nobjective function evaluations = %-i' ...
         '\niterations = %-i' ...
         '\nflag = %-i' ...
         '\n***********************************************************\n\n'], ...
         f_ord,ord_info.n_f,ord_info.it,ord_info.flag);


% To call DF-SIMPLEX as standalone, follow the same steps as above but call
% DF-SIMPLEX instead of ORD (note that, in DF-SIMPLEX, the starting point
% is passed as a vector of coefficients expressing a convex combination of
% the atoms, see the file 'usage_df_simplex.txt' for further details).
% Namely,

[x_df_simplex,f_df_simplex,df_simplex_info] = df_simplex(obj,A,double(1:m==i0)');
fprintf(['\n********************** FINAL RESULTS **********************' ...
         '\nAlgorithm: DF-SIMPLEX' ...
         '\nf =  %-.4e'   ...
         '\nobjective function evaluations = %-i' ...
         '\niterations = %-i' ...
         '\nflag = %-i' ...
         '\n***********************************************************\n\n'], ...
         f_df_simplex,df_simplex_info.n_f,df_simplex_info.it,df_simplex_info.flag);

% DF-SIMPLEX parameters can be changed as explained above for ORD, that is,
% by passing as fourth input argument a structure with the parameters to be
% changed (see the file 'usage_df_simplex.txt' to know which parameters
% can be changed and their default values)