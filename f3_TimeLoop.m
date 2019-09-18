function [T, T_graph] = f3_TimeLoop(ne,ngp,nn,xg,yg,N,w12,w,Neta,CM,BM,ENM,J,t,K,F,M,dt,alpha,tend)
T_graph = zeros(nn,size(t,2));
T = K\F; %evaluate T
T_graph(:,1) = T;
%% Time for loop
i = 2;

xygp = xg(:,:)+yg(:,:);
xygm = xg(:,:)-yg(:,:);
Jw12 = w12(:,1)*J;

switch alpha %switch to speed up process for alpha = 1 
    case 1
        for it = dt:dt:tend
            
            F = f1_BuildF(ne,nn,xygp,xygm,N,Jw12,ENM,it);
            F = f2_ApplyBC2 (ne,CM,BM,F,ngp,w,yg,Neta,it);
            T = (M+dt*K)\(dt*F+M*T); %alpha = 1 so omit
            % Graphical purpose
            T_graph(:,i) = T;   %store T value for each time frame for movie purpose
            i = i+1;            %increment i
        end
        
    otherwise
        for it = dt:dt:tend
            V = M\(F-K*T);
            T_hat = T+(1-alpha)*dt*V;
            F = f1_BuildF(ne,nn,xg,yg,N,w1,w2,ENM,J,it);
            F = f2_ApplyBC2 (ne,CM,BM,F,ngp,w,yg,Neta,it);
            T = (M+alpha*dt*K)\(alpha*dt*F+M*T_hat);
            
            % Graphical purpose
            T_graph(:,i) = T;   %store T value for each time frame for movie purpose
            i = i+1;            %increment i
        end
        
end

end
