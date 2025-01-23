% this demo is a complete test of FR, and a `complete' test of SDP up to
% some RAM constraints

% compare affine FR, partial FR (D) and partial FR (DD)
maxsize = 10000;
solveSDP = false;
BIPtable = BIP_test(maxsize,solveSDP);

% apply affine FR, partial FR (D) and partial FR (DD)
% and if there is a positive reduction, then solve Shor's relaxation
% with and without reduction for comparison.
maxsize = 350;
solveSDP = true;
BIPtable = BIP_test(maxsize,solveSDP);
