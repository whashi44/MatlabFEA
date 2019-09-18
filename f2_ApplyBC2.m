function F = f2_ApplyBC2 (ne,CM,BM,F,ngp,w,yg,Neta,t)

% This function apply Drichlet and Neuman Boundary Condition to only
% Forcing term 
% Boundary Matrix column: 1:bottom, 2:right, 3:top, 4:left 
%% Neuman Boundary Condition 
% gN = -xsin(x-1+t)         (:,1) top face 
% gN = -ysin(-y+t)          (0,:) left face 
xyg1 = yg(1:ne:ne^2,1:ngp); %mapped coordinates for 1D  

% Computing top face  
F(BM(1:end-1,3))          = F(BM(1:end-1,3)) -((-xyg1(:,:).*(sin(xyg1(:,:)-1+t)))*(w(1,:)'.*Neta(1:ngp,1))*0.5*1/ne);  
F(BM(2:end,3))            = F(BM(2:end,3))   -((-xyg1(:,:).*(sin(xyg1(:,:)-1+t)))*(w(1,:)'.*Neta(1:ngp,2))*0.5*1/ne);                      
% Computing left face 
F(BM(1:end-1,4))          = F(BM(1:end-1,4)) -((-xyg1(:,:).*(sin(-xyg1(:,:)+t))) *(w(1,:)'.*Neta(1:ngp,1))*0.5*1/ne);  
F(BM(2:end,4))            = F(BM(2:end,4))   -((-xyg1(:,:).*(sin(-xyg1(:,:)+t))) *(w(1,:)'.*Neta(1:ngp,2))*0.5*1/ne); 

%% Drichlet BC 
% gD = cos(x+t)             (:,0) bottom face
% gD = cos(1-y+t)           (1,:) right face 
F(BM(:,1),1)           = cos(CM(BM(:,1),1)+t); %bottom node forcing term 
F(BM(:,2),1)           = cos(1-CM(BM(:,2),2)+t); %top node forcing term 

end
