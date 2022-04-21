# Structured Derivative-Free Optimization & Black-Box Adversarial Attacks

_Optimize, Refine & Drop_ (ORD) is a derivative-free solver for
optimization problems of the following form:

<img src="https://latex.codecogs.com/svg.image?\min&space;f(x)&space;\\\text{s.t.&space;}&space;&space;x&space;\in&space;\text{conv}&space;\{a_1,\ldots,a_m\}">

where f(x) is a _black-box function_  (assumed to be continuously differentiable)
and conv{a<sub>1</sub>,...,a<sub>m</sub>} is the convex hull of some given vectors a<sub>1</sub>,...,a<sub>m</sub>
called _atoms_.

ORD uses an inner approximation approach that, at each iteration, approximately minimizes f(x)
by the DF-SIMPLEX algorithm, with a growing precision, over the convex hull of a suitably chosen subset of atoms,
using proper rules to add and remove atoms.

ORD can also be used for _black-box adversarial attacks_ (see below).

## Reference paper

[A. Cristofari, F. Rinaldi (2021). _A Derivative-Free Method for Structured Optimization Problems._
SIAM Journal on Optimization, 31(2), 1079-1107](https://epubs.siam.org/doi/abs/10.1137/20M1337417).

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

Copyright 2021-2022 Andrea Cristofari, Francesco Rinaldi.

## How to use ORD

1. This directory should contain the following files:

    * `COPYING.txt`,
    * `df_simplex.m`,
    * `main.m`,
    * `ord.m`,
    * `README.md`,
    * `usage_ord.txt`,
    * `usage_df_simplex.txt`.

2. See the file `usage_ord.txt` to know how to call ORD in Matlab, change
   algorithm parameters and get output values.

3. See the file `main.m` for an example.
   To run the example, just call `main.m` in Matlab.

   File `main.m` also contains an example of how to call DF-SIMPLEX as standalone, see below.

## When using ORD and when using DF-SIMPLEX

DF-SIMPLEX is the algorithm used at each iteration of ORD to solve the reduced problems,
employing sparse directions that contain positive generators of the tangent cone at the current iterate
and a specific line search.
So, DF-SIMPLEX can even be used as standalone to solve the same class of problems as ORD.
To call DF-SIMPLEX as standalone, see the file `usage_df_simplex.txt` and the example in the file `main.m`.

Since ORD and DF-SIMPLEX follow different approaches, in general terms we can say that
ORD is preferable when the optimal solutions can be expressed as the convex combination
of just a few atoms. For example, this is the case when the number of atoms is much larger than
the number of variables, or when the problem has a a particular structure that induces sparsity,
such as l1-norm constraint.

## An application: black-box adversarial attacks

As an example of application of ORD, consider black-box adversarial attacks.
In these problems, the goal is to perturb the inputs of a given classifier to generate samples that lead to
misclassification.
In particular, the _maximum allowable l1-norm attack_ can be formulated as

<img src="https://latex.codecogs.com/svg.image?\min&space;f(x_0&plus;x)&space;\\\text{s.t.&space;}&space;||x||_1&space;\le&space;\varepsilon">

where f is a suitably chosen attack loss function, x<sub>0</sub> is a vector representing a correctly classified sample,
||x||<sub>1</sub> is the &ell;<sub>1</sub>-norm of x and &epsilon; is a positive parameter.
Note that f is a black-box function when the internal configuration of the classifier is unknown.
Since the feasible set can be expressed as the convex combination of the vertices of the &ell;<sub>1</sub>-ball,
ORD can be used to solve this class of problems.

When considering image classification, usually there is an additional constraint of the form
0 &le; x<sub>0</sub> + x &le; 1, which can be removed by applying a proper variable transformation.
Further details can be found in the reference paper [(Cristofari, Rinaldi, 2021)](https://epubs.siam.org/doi/abs/10.1137/20M1337417) and in the references therein.