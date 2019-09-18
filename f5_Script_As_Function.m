function [L2,nn] =f5_Script_As_Function(ne,ngp,timestep,alpha)

%unchangeable variables
deg = 1;        % degree of the shape function
nl  = 4;        % number of local nodes, quadrilateral = 4
h   = [0 1];    % spatial domain
tdomain = [0 1];
%% Derived variables
nn  = ((ne*deg)+1)^2; % number of total nodes
dt  = 1/timestep;        % step size
tend = tdomain(2);
t = 0:dt:tend;


% Calculated using functions
[xi,eta,w12,w]        = f0_gauss(ngp);                %abscissae and weight
[N,Dxi,Deta,Neta]       = f0_shape(ngp,nl,xi,eta) ;     %shape function and derivative of shape function
[CM,ENM,BM]             = f0_ID(ne,nl,h,nn);            %Coordinate Matrix (CM), Elemental Node Matrix (ENM) and Boundary Matrix (BM)
[x, y, xg, yg]          = f0_mapXY(ne,ngp,nl,CM,ENM,N); %local x y coordinates, mapped xg yg coordinates

%Jacobian is element and gaussian point independent, so we only need one, check read me for the detail of this line
J   = det([x(:,1)'*Dxi(:,1) x(:,1)'*Deta(:,1); y(:,1)'*Dxi(:,1) y(:,1)'*Deta(:,1)]);

%% Buildkfm
[K,F,M] = f1_BuildKFM(ne,ngp,nl,nn,xg,yg,N,Dxi,Deta,w12,ENM,J,t(1));

%% Apply Drichlet and Neuman Boundary Condition
[K,F] = f2_ApplyBC (ne,CM,BM,K,F,ngp,w,yg,Neta,t(1));

%% Calculate Temperature
[T, T_graph] = f3_TimeLoop(ne,ngp,nn,xg,yg,N,w12,w,Neta,CM,BM,ENM,J,t,K,F,M,dt,alpha,tend);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1 = toc;     % end measuring
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Plot Results
% 

%%calculate error
L2 = f4_CalcError(ne,N,xg,yg,w12,ENM,T,J,tend);
fprintf(' The L2 error for %d x %d x %d (ne by ne by timestep)is: %d \n',ne,ne,timestep,L2);

end 
