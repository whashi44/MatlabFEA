function [N,Dxi,Deta,Neta] = f0_shape(ngp,nl,xi,eta)
% This function calculate the shape function and derivative of shape function evaluated at the gaussian points
% This shape function will only work for deg of 1

%% Preallocation (row:nl, column:ngp^2)
N     = zeros(nl,ngp^2); % shape function
Dxi   = zeros(nl,ngp^2); % derivative of shape function with respect to xi
Deta  = zeros(nl,ngp^2); % derivative of shape function with respect to eta

%% Assign shape function
% Initilize
Nxi     = 0.5*[1-xi(:) 1+xi(:)];    % shape function for xi direction in master element, [N1xi N2xi]
Neta    = 0.5*[1-eta(:) 1+eta(:)];  % shape function for eta direction in master element, [N1eat N2eta]
D       = 0.5*[-1 1];               % derivative of the shape fucntion in both xi and eta, which are the same so we only need 1 row [D1 D2]
SFM     = [1 2 2 1; 1 1 2 2]';      % SFM = Shape Function Matrix to assign shape function

%% Calculate and evaluate shape functions
% Using the SFM, we can calculate 2D shape function with specific
% combination 
N(:,:)      = Nxi(:,SFM(:,1))'.*Neta(:,SFM(:,2))';  % shape unction in 2D 
Dxi(:,:)    = D(SFM(:,1))'.*Neta(:,SFM(:,2))';      % derivative of shape function in 2D in xi direction 
Deta(:,:)   = Nxi(:,SFM(:,1))'.*D(SFM(:,2))';       % derivative of shape function in 2D in eta direction 

end