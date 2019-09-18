function F = f1_BuildF(ne,nn,xygp,xygm,N,Jw12,ENM,it)
% This function build F matrix for the time loop
%% Preallocatoin
F = zeros(nn,1);            % preallocation global force matrix
%% Evaluation and Assembly
for ie = 1:ne^2
  % F(ENM(:,ie),1)  = F(ENM(:,ie),1)+N(:,:)*(((xg(ie,:)+yg(ie,:)).*cos(xg(ie,:)-yg(ie,:)+it)-sin(xg(ie,:)-yg(ie,:)+it))'.*w12(:,1)*J);
    F(ENM(:,ie),1)  = F(ENM(:,ie),1)+N(:,:)*(((xygp(ie,:))       .*cos(xygm(ie,:)+it)       -sin(xygm(ie,:)+it))'       .*Jw12(:,1));  
end

end
