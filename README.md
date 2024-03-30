# approxEllipsotopeContainment

Code generating the figures and data for the paper "Approximability of the Containment Problem for Zonotopes and Ellipsotopes".

This code requires versions of Mosek, Yalmip, CORA, and AROC that were available at the time of writing this readme (this repository will NOT be updated in the future, as these functionalities will be integrated into CORA and AROC).

Also, note that even with those libraries installed, you will get an error due to AROC; this is because CORA 2024 changed the mptPolytope class to the (better) polytope class, therefore one needs to make some arrangements in the AROC code: whenever a mptPolytope object is created, one needs to replace 'mptPolytope' by 'polytope'. This will be corrected in the next AROC version.