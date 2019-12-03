function Iwind = makeWindFun(SZX,SZY)
% This is just a helper function to make a random scalar function. This is
% called twice to generate a random wind field.
%
% Copyright (c) 2012, MathWorks, Inc. 
%

windFineness = 0.1;
if ~exist('SZX','var')
SZX = 50;
SZY = 50;
end

N = 50; % Various parameters used in generating a random "smooth" matrix
NL = 40; 
NP = 500;
rx = randn(NL,N);
rx = interpft(rx,NP);
ry = randn(NL,N);
ry = interpft(ry,NP);
I = (rx*ry');

[xgi,ygi] = meshgrid(linspace(1,2 + 498*windFineness,SZX+1),linspace(1,2 + 498*windFineness,SZY+1));
Iwind = 10*interp2(1:500,1:500,I,xgi,ygi);

