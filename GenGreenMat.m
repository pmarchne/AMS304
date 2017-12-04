function [GreenMat] = GenGreenMat(NodesCoor,k)
%Generate the matrix associated to the Green Kernel
N = length(NodesCoor);

% Extract x,y and z coordinates
XCoor = meshgrid(NodesCoor(:,1));
YCoor = meshgrid(NodesCoor(:,2));
ZCoor = meshgrid(NodesCoor(:,3));

% compute the interactions between nodes for the 3 spatial coordinates x,y and z
Sx = XCoor' - XCoor; 
Sy = YCoor' - YCoor; 
Sz = ZCoor' - ZCoor;

% the final matrix should have zero diagonal terms
S = sqrt(Sx.^2 + Sy.^2 + Sz.^2) ;%- eye(N);
GreenMat = exp(1i*k*S)./S ;% + eye(N)*exp(-1i * k);

end

