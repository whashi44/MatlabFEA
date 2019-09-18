function L2 = f4_CalcError(ne,N,xg,yg,w12,ENM,T,J,t)
% This function calculates L2 error 
L2 = 0;
for ie = 1:ne^2
    T_a(:,1) = T(ENM(:,ie),1); %approximate solution     
    T_t(:,1) = cos(xg(ie,:)-yg(ie,:)+t); %theoretical solution 
    L2 = L2 + J*((T_t-N(:,:)'*T_a(:,1)).^2)'*w12(:,1);
    %          *( 9x1-(4x9)' *4x1         )'*(9x1    .*9x1    )  
end

L2 = sqrt(L2);

end
