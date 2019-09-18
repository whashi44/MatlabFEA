function [K,F,M] = f1_BuildKFM(ne,ngp,nl,nn,xg,yg,N,Dxi,Deta,w12,ENM,J,t)
% This function build K and F matrix
%% Preallocatoin
F = zeros(nn,1);            % preallocation global force matrix
k = zeros(2,2);             % preallocation k tensor

%% Evaluation and Assembly using proper sparce matrix
% Preallocation
counter = 1;
R = zeros(nl^2*ne^2,1);   % Row vector
C = zeros(nl^2*ne^2,1);   % Column vector
Kval = zeros(nl^2*ne^2,1);% stiffness value
MVal = zeros(nl^2*ne^2,1);% mass value
D1 = zeros(nl,2,ngp^2);   % derivative of shape function 1
D2 = zeros(2,nl,ngp^2);   % derivative of shape function 2

% Prepare derivative of shape function
for ig = 1:ngp^2
    D1(:,:,ig) = [Dxi(:,ig) Deta(:,ig)];
    D2(:,:,ig) = [Dxi(:,ig)';Deta(:,ig)'];
end

%% Assembly
for ie = 1:ne^2
    ke = zeros(nl);            % preallocate local stiffness matrix
    for ig = 1:ngp^2
        k(1,1) = yg(ie,ig);    % k tensor y value
        k(2,2) = xg(ie,ig);    % k tensor x value
        ke(:,:)   = ke(:,:)   + w12(ig)*(D1(:,:,ig) * k *  D2(:,:,ig)); % local stiffness matrix
        % ke(il,jl) = ke(il,jl) +
    end
    
    for il = 1:nl
        for jl = 1:nl
            R(counter) = ENM(il,ie);  %store row information
            C(counter) = ENM(jl,ie);  %store column information
            Kval(counter) = ke(il,jl);%store ke value
            MVal(counter) =  N(jl,:).*N(il,:)*w12(:,1)*J;  %store me value
            counter = counter+1;
        end
    end
    
    F(ENM(:,ie),1)          = F(ENM(:,ie),1)+N(:,:)*(((xg(ie,:)+yg(ie,:)).*cos(xg(ie,:)-yg(ie,:)+t))'.*(w12(:,1)))*J;
    %                                4x1    + 4x9  *(       1x9+1x9      .*         1x9-1x9        )'.*9x1       )*1
end

K = sparse(R,C,Kval,nn,nn); %create sparse matrix for optimization of performance
M = sparse(R,C,MVal,nn,nn); %sparse(row, column, value, matrix dim, matrix dim)

end
