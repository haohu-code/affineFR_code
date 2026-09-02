# Third-party software and benchmark data

The MIT License in `LICENSE` applies only to the original Affine FR code in
this package. Third-party software and benchmark data remain subject to their
respective copyright notices, licenses, and terms of use.

## YALMIP

This package includes YALMIP version 20230622 in `YALMIP-master/`. YALMIP is
copyright Johan Löfberg and is distributed under the terms in
`YALMIP-master/license.txt`; it is not covered by this package's MIT License.
The YALMIP license also requires YALMIP to be referenced in published work.

Suggested reference:

> J. Löfberg, "YALMIP: A Toolbox for Modeling and Optimization in MATLAB,"
> Proceedings of the 2004 IEEE International Symposium on Computer Aided
> Control Systems Design, pp. 284--289, 2004.

Project website: <https://yalmip.github.io/>

## MIPLIB 2017

The MIPLIB 2017 instance files required by the MIPLIB experiments are not
distributed with this package. Users must obtain them from the official
MIPLIB website and place them in `MIPLIB/data/` as described in
`MIPLIB/README.md`. The downloaded files are not covered by this package's MIT
License: <https://miplib.zib.de/>.

Reference:

> A. Gleixner et al., "MIPLIB 2017: Data-Driven Compilation of the 6th
> Mixed-Integer Programming Library," Mathematical Programming Computation,
> 2021. <https://doi.org/10.1007/s12532-020-00194-3>

## SATLIB

The files in `SAT/data/` include benchmark instances obtained from SATLIB.
They are not covered by this package's MIT License. SATLIB asks researchers
who use the library to cite the following publication:

> H. H. Hoos and T. Stützle, "SATLIB: An Online Resource for Research on
> SAT," in SAT 2000, I. P. Gent, H. van Maaren, and T. Walsh, eds.,
> pp. 283--292, IOS Press, 2000.

SATLIB website: <https://www.cs.ubc.ca/~hoos/SATLIB/>

## External solvers

Gurobi and Mosek are runtime dependencies but are not included in this
package. Each must be installed separately and is governed by its provider's
license and terms.

- Gurobi: <https://www.gurobi.com/>
- Mosek: <https://www.mosek.com/>
