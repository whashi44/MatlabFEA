function [x, y, xg, yg] = f0_mapXY(ne,ngp,nl,CM,ENM,N)
%{ 
this function grab coordiantes for each node for each element and map it for
gaussian integration. 
Input: ne, ngp, nl, CM, ENM, N
Outpu: x, y, xg, yg 
%} 

%% Preallocation
x   = zeros(nl,ne^2);         % element node x coordinates
y   = zeros(nl,ne^2);         % element node y coordinates
xg  = zeros(ne^2,ngp^2);      % gaussian mapped x
yg  = zeros(ne^2,ngp^2);      % gaussian mapped y

%% Assembly
%{
Assemblying coordinate for node 1-4 for each element
This looks confusing but since we know that number of local node is always 4,
 it works and this way we can speed up the code even more! 
%} 
x(:,:) =[CM(ENM(1,:),1) CM(ENM(2,:),1) CM(ENM(3,:),1) CM(ENM(4,:),1)]';
y(:,:) =[CM(ENM(1,:),2) CM(ENM(2,:),2) CM(ENM(3,:),2) CM(ENM(4,:),2)]';
xg(:,:) = x(:,:)'*N(:,:); %(nl by ne^2)'*(nl by ngp^2)
yg(:,:) = y(:,:)'*N(:,:); %(nl by ne^2)'*(nl by ngp^2)
 end
