# Structured Derivative-Free Optimization & Black-Box Adversarial Attacks

_Optimize, Refine & Drop_ (ORD) is a derivative-free solver for
optimization problems of the following form:

         min f(x)
    s.t. x in conv{a_1,...,a_m}

where _f(x)_ is a (black-box) continuously differentiable function and
_conv{a\_1,...,a\_m}_ is the convex hull of some given vectors _a\_1,...,a\_m,_
called _atoms_.

ORD uses an inner approximation approach that, at each iteration, approximately minimizes _f(x)_
by the DF-SIMPLEX algorithm, with a growing precision, over the convex hull of a suitably chosen subset of atoms,
using proper rules to add and remove atoms.

ORD can be also used to solve black-box adversarial attacks (see below).

## Reference paper

A. Cristofari, F. Rinaldi (2021). _A derivative-free method for structured
optimization problems._ SIAM Journal on Optimization, to appear. Pre-print
available at <https://arxiv.org/abs/2005.05224>.

## Authors

* Andrea Cristofari (e-mail: [andrea.cristofari@unipd.it](mailto:andrea.cristofari@unipd.it))
* Francesco Rinaldi (e-mail: [rinaldi@math.unipd.it](mailto:rinaldi@math.unipd.it))

## Licensing

ORD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
ORD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with ORD. If not, see <http://www.gnu.org/licenses/>.

Copyright 2021 Andrea Cristofari, Francesco Rinaldi.

## How to use ORD

1. This directory should contain the following files:

    * `COPYING.txt`,
    * `syntax_ord.txt`,
    * `syntax_df_simplex.txt`,
    * `README.md`,
    * `ORD.m`,
    * `main.m`,
    * `DF_SIMPLEX.m`.

2. See the file `syntax_ord.txt` to know how to call ORD, change algorithm parameters and get output values.

3. See the file `main.m` for an example of how to call ORD to solve a user-defined problem
  (this file also contains an example of how to call DF-SIMPLEX as standalone, see below).

## When using ORD and when using DF-SIMPLEX

DF-SIMPLEX is the algorithm used at each iteration of ORD to solve the reduced problems,
employing sparse directions that contain positive generators of the tangent cone at the current iterate
and a specific line search.
So, DF-SIMPLEX can even be used as standalone to solve the same class of problems as ORD.
To call DF-SIMPLEX as standalone, see the file `syntax_df_simplex.txt` and the example in the file `main.m`.

Since ORD and DF-SIMPLEX follow different approaches, in general terms we can say that
ORD is preferable when the optimal solutions can be expressed as the convex combination
of just a few atoms. For example, this is the case when the number of atoms is much larger than
the number of variables, or when the problem has a a particular structure that induces sparsity,
such as l1-norm constraint.

## An application: black-box adversarial attacks

As an example of application of ORD, consider black-box adversarial attacks.
In this problems, the goal is to perturb the inputs of a given classifier to generate samples that lead to
misclassification.
In particular, the _maximum allowable l1-norm attack_ can be formulated as

         min f(x0 + x)
    s.t. ||x||_1 <= eps

where _f_ is a suitably chosen attack loss function, _x0_ is a correctly classified sample,
_||x||\_1_ is the _l_1-norm of _x_ and _eps_ is a  positive parameter.
Note that _f_ is a black-box function when the internal configuration of the classifier is unknown.
Moreover, the feasible set can be expressed as the convex combination of the vertices of the _l_1-ball.
So, ORD can be used to solve this class of problems.

When considering image classification, usually there is an additional constraint of the form
_0 <= x0 + x <= 1_, which can be removed by applying a proper variable transformation.
Further details can be found in the reference paper (Cristofari, Rinaldi, 2021) and in the references therein.