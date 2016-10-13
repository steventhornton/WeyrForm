# Weyr Canonical Form
This repository contains a [Maple](http://www.maplesoft.com/products/maple/) function for computing the [Weyr canonical form](https://wikipedia.org/wiki/Weyr_canonical_form) of a square matrix.

The function `WeyrForm` expects a square matrix as input. It will return the Weyr form of the input by default. The optional output argument can be used to specify what should be returned:
- `output = W` will return the Weyr form
- `output = Q` will return an invertible matrix `Q` such that `MatrixInverse(Q).A.Q = W`
- `output = [W, Q]` will return both the Weyr form and invertible matrix `Q`
- `output = [Q, W]` will return the invertible matrix followed by the Weyr form

The `Example.mw` Maple worksheet includes 4 examples.

A poster explaining how this implementation works can be found [here](https://s3.amazonaws.com/stevenethornton.github/WeyrForm.pdf).

## Example
```
read("WeyrForm.mpl");
with(LinearAlgebra):

A := Matrix([[-1,  0,  0,  2,  1,  0,  0],
             [-1,  0,  0,  1,  1,  0,  0],
             [ 6, -2, -2,  4,  2, -2,  4],
             [-2,  1,  1,  0,  0,  1, -1],
             [ 3, -1, -1,  2,  1, -1,  2],
             [-5,  2,  2, -5, -3,  2, -4],
             [ 2,  0,  0, -4, -2,  0,  0]]);

Q, W := WeyrForm(A, output = [Q, W]);
```

## Relevant References
- O'Meara, K. et al (2011). Advanced Topics in Linear Algebra: Weaving Matrix Problems through the Weyr Form. Oxford University Press, USA.
- Helene Shapiro. The Weyr Characteristic. The American Mathematical Monthly, 106(10):919–929, 1999.
- Eduard Weyr. Répartition des matrices en espèces et formation de toutes les espèces. CR Acad. Sci. Paris, 100:966–969, 1885.

## License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License  along with this program.  If not, see http://www.gnu.org/licenses/.
