This hash uses a tree structure to efficiently hash unit vectors that fall in a region defined by the paper :

https://www.semanticscholar.org/paper/A-PARTITION-OF-THE-UNIT-SPHERE-INTO-REGIONS-OF-AREA-Leopardi/9cb561cd1f0430d5556c1c128bbe6baaeda7edca

This algorithm provides an implementation of partitioning the unit sphere into regions, and the hash mechanism hashes unit vectors within the same region
into the same bucket. The hash of a unit vector is of the form "[a1,b1][a2,b2]...[an,bn]", where ai and bi are the left and right angular bounds of dimension i.
Note that these bounds are angular, because we are using n-dimensional spherical coordinates.
