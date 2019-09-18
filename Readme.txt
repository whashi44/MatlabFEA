%%%%%%%%%%%%%%%%%%%%%%
This program approximates the solution to the time-dependent second-order
partial differential heat equation using Strong form, Weak form, divergent
theorem, discretization, gallerkin's approximation, and Gaussian quadrature. 



%%%%%%%%%%%%%%%%%%%%%%%
% Variable definition %
%%%%%%%%%%%%%%%%%%%%%%%
ne    = number of elements 
ngp   = number of gaussian points 
deg   = degree of shape function 
nl    = number of local nodes 
h     = domain of the problems
nn    = total number of nodes 
J     = Jacobian 
L2    = L2 error 
ie    = element index 
il    = i node index
jl    = j node index 
ig    = gaussian index 
xi    = (ngp^2 x 1) 	gaussian point in x coordinate for master element		
ea    = (ngp^2 x 1) 	gaussian point in y coordinate for master element		
w1    = (ngp^2 x 1) 	weight in x direction for master element 	  		
w2    = (ngp^2 x 1) 	weight in y direction for master element 	   
N     = (nl x ngp^2) 	shape function 					  		
Dxi   = (nl x ngp^2) 	derivative of shape function with respect to xi	  
Deta  = (nl x ngp^2) 	derivative of shape function with respect to eta  
CM    = (nn x 2)        Coordinate Matrix that stores global coordinate  		
ENM   = (nl x ne^2) 	Elemental Node Matrix for assemlying local to global matrix 
BM    = (ne+1 x face)	Boundary Matrix to assign boundary condition 
x     = (nl x ne^2)     local x coordinates for each node for each element
y     = (nl x ne^2) 	local y coordinates for each node for each element 
xg    = (ne^2 x ngp^2)  x coordinates mapped for gaussian integration 
yg    = (ne^2 x ngp^2)  y coordinates mapped for gaussian integration 
K     = (nn x nn) 	global stiffness matrix 
F     = (nn x 1)	global force matrix 
M     = (nn x nn) 	global mass matrix 
T     = (nn x 1)	temeprature matrix 
ke    = (nl x nl) 	local stiffness matrix 
fe    = (nl x 1)	local force matrix 
me    = (nl x 1)        local mass matrix 
timestep = number of division user want for the time steps 
alpha = constant used to decide which time discretization method user wants
dt = time steps 

------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%
%   Functions  %
%%%%%%%%%%%%%%%%

f0_gauss.m
Description: 
Grab the abscissae and weight from 1D gaussian table and transform it into 2D 

Input: 
number of gaussian points (ngp)

Output: 
abscissae x direction (xi)
abscissae y directoin (eta) 
weight x direction (w1)
weight y direction (w2) 

Required function:  
None 

Comment: 


------------------------------------------------------------------------------------------
f0_shape.m
Description: 
calculate the shape function matrices 

Input: 
number of gaussian points (ngp), number of local nodes (nl) 

Output: 
shape function matrix (N)(nl by ngp^2)
derivative of shape function with respect to xi (Dxi)(nl by ngp^2) 
derivative of shape function with respect to eta (Deta)(nl by ngp^2) 

comment: 
Only work for degree 1 

shape function explanation 
Nxi = 0.5*[1-xi(:)  1+xi(:)];  %shape function for xi direction in master element, [N1xi N2xi]
Neta = 0.5*[1-eta(:) 1+eta(:)];%shape function for eta direction in master element, [N1eat N2eta]
D   = 0.5*[-1 1];              %derivative of the shape fucntion in both xi and eta, which are the same [D1 D2]
SFM = [1 2 2 1; 1 1 2 2]';     %shape function matrix to assign shape function

N(il,:)   = Nxi(:,SFM(il,1)).*Neta(:,SFM(il,2)); 
Dxi(il,:) = D(SFM(il,1)).*Neta(:,SFM(il,2));
Deta(il,:) = Nxi(:,SFM(il,1)).*D(SFM(il,2));

ABOVE CODE IS SAME AS BELOW BUT FASTER --------

% Anonymous function  
N1 =@(x) 0.5*(1-x);
N2 =@(x) 0.5*(1+x); 
D1 =-0.5;
D2 =0.5;

for i = 1:ngp^2     %for every gaussian point 
    % shape function 
    N(1,i) =N1(xi(i))*N1(eta(i)) ;
    N(2,i) =N2(xi(i))*N1(eta(i)) ;
    N(3,i) =N2(xi(i))*N2(eta(i)) ;
    N(4,i) =N1(xi(i))*N2(eta(i)) ;
    
    % derivative of shape function with respect to xi 
    Dxi(1,i) = D1*N1(eta(i));  
    Dxi(2,i) = D2*N1(eta(i));
    Dxi(3,i) = D2*N2(eta(i));
    Dxi(4,i) = D1*N2(eta(i));
    
    % derivative of shape function with respect to eta
    Deta(1,i) = N1(xi(i))*D1;
    Deta(2,i) = N2(xi(i))*D1;
    Deta(3,i) = N2(xi(i))*D2;
    Deta(4,i) = N1(xi(i))*D2; 
end 

end
------------------------------------------------------------------------------------------
f0_ID.m
Description: 
calculate 3 matrices that help us identify points and boundaries 

Comment: 
ENM = zeros(nl,ne^2+ne-1); %create an array 
j = 1; 
ENM(1,:) = (j:(ne+1)^2-(ne+2));
ENM(2,:) = (j+1:(ne+1)^2-(ne+1));
ENM(3,:) = (j+ne+2:(ne+1)^2);
ENM(4,:) = (j+ne+1:(ne+1)^2-1);  
ENM(:,ne+1:ne+1:end)=[];  %remove uneeded columns 

ABOVE CODE IS SAME AS BELOW BUT FASTER --------

%% Elemental Node Matrix 
ENM = zeros(nl,ne^2);
x = pts(1);
y = pts(1);
for ie =1:ne^2
    for jl = 1:nl        
        ENM(jl,ie) = find(abs(CM(:,1)-x)<1e-5 & abs(CM(:,2)-y)<1e-5); %finding the desired value to assign coordinates 
        if jl == 1
            x = x+(1/ne); %go right
        elseif jl == 2
            y = y+(1/ne); %go up
        elseif jl == 3
            x = x-(1/ne); %go left
        elseif jl == 4
            y=y-(1/ne);  %go down
            x=x+(1/ne);  %go right
        end
    end    
    if mod(ie,ne) ==0 %if the current
        x = pts(1); %go to the left most node
        y = y+(1/ne); %go up       
    end
end


------------------------------------------------------------------------------------------
f0_mapXY
Decription:
Allocate local coodinate (x,y) for each element and map (x,y) for gaussian integration (xg,yg) 

Input:

Output:


Comment: 
%Line clarification 


x(:,:) =[CM(ENM(1,:),1) CM(ENM(2,:),1) CM(ENM(3,:),1) CM(ENM(4,:),1)]';
y(:,:) =[CM(ENM(1,:),2) CM(ENM(2,:),2) CM(ENM(3,:),2) CM(ENM(4,:),2)]';
xg(:,:) = x(:,:)'*N(:,:); %(nl by ne^2)'*(nl by ngp^2)
yg(:,:) = y(:,:)'*N(:,:); %(nl by ne^2)'*(nl by ngp^2)


ABOVE CODE IS SAME AS BELOW BUT FASTER --------


x(:,:) =[CM(ENM(1,:),1)';CM(ENM(2,:),1)';CM(ENM(3,:),1)';CM(ENM(4,:),1)'];
y(:,:) =[CM(ENM(1,:),2)';CM(ENM(2,:),2)';CM(ENM(3,:),2)';CM(ENM(4,:),2)'];
xg(:,:) = x(:,:)'*N(:,:); %(nl by ne^2)'*(nl by ngp^2)
yg(:,:) = y(:,:)'*N(:,:); %(nl by ne^2)'*(nl by ngp^2)

ABOVE CODE IS SAME AS BELOW BUT FASTER --------

x(1,:) = CM(ENM(1,:),1);   % 1st node, second 1 indicate x coordinate
x(2,:) = CM(ENM(2,:),1);   
x(3,:) = CM(ENM(3,:),1); 
x(4,:) = CM(ENM(4,:),1); 
y(1,:) = CM(ENM(1,:),2);  % 1st node, 2 indicate y coordinates 
y(2,:) = CM(ENM(2,:),2); 
y(3,:) = CM(ENM(3,:),2); 
y(4,:) = CM(ENM(4,:),2); 
xg(:,:) = x(:,:)'*N(:,:); %(nl by ne^2)'*(nl by ngp^2)
yg(:,:) = y(:,:)'*N(:,:); %(nl by ne^2)'*(nl by ngp^2)

ABOVE CODE IS SAME AS BELOW BUT FASTER --------

for ie=1:ne^2 %for each element
      x(:,ie)   = CM(ENM(:,ie),1);  % x(nl,ie) = CM(ENM(nl,ie),1) 1 = x coordinates
     y(:,ie)   = CM(ENM(:,ie),2);  % y(nl,ie) = CM(ENM(nl,ie),2) 2 = y coordinates
     xg(ie,:)  = N(:,:)'*x(:,ie);  % xg(ie,npg^2) = N(nl,ngp^2)'*x(nl,ie)
      % ne^2 by ngp^2     nl*ne^2' * nl by ngp^  
      yg(ie,:)  = N(:,:)'*y(:,ie);  % xg(ie,npg^2) = N(nl,ngp^2)'*x(nl,ie)
end % ie

for ie=1:ne^2
    x(:,ie) = CM(ENM(:,ie),1); 
    y(:,ie) = CM(ENM(:,ie),2); 
    xg (ie,:) = N'*x(:,ie);    
    yg (ie,:) = N'*y(:,ie);   
end

ABOVE CODE IS SAME AS BELOW BUT FASTER --------
 for ie =1:ne^2 
     for il = 1:nl 
          x(il,ie) = CM(ENM(il,ie),1); 
          y(il,ie) = CM(ENM(il,ie),2); 
      end
end 
 for ie =1:ne^2 
     for ig = 1:ngp^2
           xg(ie,ig) = N(:,ig)'*x(:,ie);
           yg(ie,ig) = N(:,ig)'*y(:,ie);
     end 
end 

 but instead of doing this loop, we can simply loop over the number of element 
 while performing matrix calculation. 
 x(:,ie): (row:nl, column:element number) x coordinates for all the local nodes for element ie 
 y(:,ie): same as x but y coordinates 
 ENM(:,ie): (row:local nodes, column:element number) elemental node matrix for all the local nodes for element ie 
 CM(ENM(:,ie),1): (row:global coordinate, column: coordinate of x) for all global coordinate x that is in this local nodes for element ie 
 CM(ENM(:,ie),2): (row:global coordinate, column: coordinate of y) for all global coordinate y that is in this local nodes for element ie 







------------------------------------------------------------------------------------------
f1_BuildKFM
Decription:
1. calculate local stiffness, force and mass matrix 
2. assign them to global stifness, force and mass matrix
3. create sparse matrix of K & M for performance optimization 

Comment: 
tensor explanation from buildKFM
k = [xg(ie,ig) 0;0 yg(ie,ig)]; 

k tensor is an array, but creating array for every loop is computationally expensive, 
so I simply substitute value for only first and last cell. 

k(2:3) = 0; %this is already 0 due to preallocation, so we don't need it. 




counter = 0;
R = zeros(nl^2*ne^2,1);
C = zeros(nl^2*ne^2,1);
Kval = zeros(nl^2*ne^2,1);
MVal = zeros(nl^2*ne^2,1);


for ie = 1:ne^2
    ke = zeros(nl);            % preallocate local stiffness matrix
    me = zeros(nl);
    
     for ig = 1:ngp^2
        k(1,1) = yg(ie,ig);    % k tensor y value
        k(2,2) = xg(ie,ig);  % k tensor x value
        ke(:,:) = ke(:,:) + w1(ig)*w2(ig)*([Dxi(:,ig) Deta(:,ig)] * k *  [Dxi(:,ig) Deta(:,ig)]'); % local stiffness matrix
     end %ig

    for il = 1:nl
        for jl = 1:nl
            counter = counter+1;       
            me(il,jl) = me(il,jl) + N(jl,:).*N(il,:)*(w1(:,1).*w2(:,1))*J;
            R(counter) = ENM(il,ie); %store row information 
            C(counter) = ENM(jl,ie); %store column information 
            Kval(counter) = ke(il,jl);  %store ke value          
            MVal(counter) = me(il,jl);  %store me value 
        end
    end
    % F(ENM(il,ie),1)=F(ENM(il,ie),1)+fe(il);
    F(ENM(:,ie),1)          = F(ENM(:,ie),1)+((((xg(ie,:)+yg(ie,:)).*cos(xg(ie,:)-yg(ie,:)+t)).*w1(:,1)'.*w2(:,1)')*N(:,:)')'*J;
    %                                             4x1    +  ((   1x9  +   1x9  ) .*         1x9         )   .*    9x1'.*    9x1' *(4x9 )')'1x1    )
    
    
    
end
 K = sparse(R,C,Kval,nn,nn); %create sparse matrix for optimization of performance
 M = sparse(R,C,MVal,nn,nn); %create sparse matrix for optimization of performance

ABOVE CODE IS SAME AS BELOW BUT FASTER --------


for ie = 1:ne^2
    ke = zeros(nl);            % preallocate local stiffness matrix
    me = zeros(nl);

    for ig = 1:ngp^2
        k(1,1) = yg(ie,ig);    % k tensor y value
        k(2,2) = xg(ie,ig);  % k tensor x value
        ke(:,:) = ke(:,:) + w1(ig)*w2(ig)*([Dxi(:,ig) Deta(:,ig)] * k *  [Dxi(:,ig) Deta(:,ig)]'); % local stiffness matrix
     end %ig

    for il = 1:nl
        for jl = 1:nl
            me(il,jl) = me(il,jl) + N(jl,:).*N(il,:)*(w1(:,1).*w2(:,1))*J;
         end
    end

    %me(:,:) = N(:,:)*(w1(:,1)*w2(:,1)')*N(:,:)'*J; %didn't work 
    %me(:,:) = N(:,:)*(w1(:,1).*w2(:,1))*(N(:,:)*(w1(:,1).*w2(:,1)))'*J; %didn't work 
    % Global Assembly
    % K(ENM(il,ie),ENM(jl,ie))=K(ENM(il,ie),ENM(jl,ie))+ke(il,jl);
    K(ENM(:,ie),ENM(:,ie))  = K(ENM(:,ie),ENM(:,ie))+ke(:,:);

    M(ENM(:,ie),ENM(:,ie))  = M(ENM(:,ie),ENM(:,ie))+me(:,:);
    % F(ENM(il,ie),1)=F(ENM(il,ie),1)+fe(il);
    F(ENM(:,ie),1)          = F(ENM(:,ie),1)+((((xg(ie,:)+yg(ie,:)).*cos(xg(ie,:)-yg(ie,:)+t)).*w1(:,1)'.*w2(:,1)')*N(:,:)')'*J;
    %                                             4x1    +  ((   1x9  +   1x9  ) .*         1x9         )   .*    9x1'.*    9x1' *(4x9 )')'1x1    )



end
 K = sparse(K); %create sparse matrix for optimization of performance
 M = sparse(M); %create sparse matrix for optimization of performance

 
ABOVE CODE IS SAME AS BELOW BUT FASTER --------

%% Preallocatoin
K = zeros(nn,nn);           % preallocation global stiffness matrix
M = zeros(nn,nn);
F = zeros(nn,1);            % preallocation global force matrix
k = zeros(2,2);             % preallocation k tensor
%% Evaluation and Assembly
for ie = 1:ne^2
    ke = zeros(nl);            % preallocate local stiffness matrix
    me = zeros(nl);
    
    for ig = 1:ngp^2
        k(1,1) = yg(ie,ig);    % k tensor y value
        k(2,2) = xg(ie,ig);  % k tensor x value
        %ke(:,:) = ke(:,:) + w1(ig)*w2(ig)*(Di' * k *  Dj); % local stiffness matrix
        ke(:,:) = ke(:,:) + w1(ig)*w2(ig)*([Dxi(:,ig) Deta(:,ig)] * k *  [Dxi(:,ig) Deta(:,ig)]'); % local stiffness matrix
        %ke(:,:) = ke(:,:) + w1(ig)*w2(ig)*(D * k *  D'); 
        %  nlxnl      nlxnl     % ngp^2x1 * ngp^2x1  *(2xnl'*2x2*2xnl)
        %me(:,:)= me(:,:) + J*w1(ig)*w2(ig)*N(:,ig)*N(:,ig)'; %somehow, for loop in ngp is slower than for loop in il and jl 
        % 4x4    4x4       1*   1  *   1  * 4x1   *  (4x1)'
    end %ig
    %
    for il = 1:nl
        for jl = 1:nl
            me(il,jl) = me(il,jl) + N(jl,:).*N(il,:)*(w1(:,1).*w2(:,1))*J;
        end
    end
    
    % Global Assembly
    % K(ENM(il,ie),ENM(jl,ie))=K(ENM(il,ie),ENM(jl,ie))+ke(il,jl);
    K(ENM(:,ie),ENM(:,ie))  = K(ENM(:,ie),ENM(:,ie))+ke(:,:);
    M(ENM(:,ie),ENM(:,ie))  = M(ENM(:,ie),ENM(:,ie))+me(:,:);
    % F(ENM(il,ie),1)=F(ENM(il,ie),1)+fe(il);
    F(ENM(:,ie),1)          = F(ENM(:,ie),1)+((((xg(ie,:)+yg(ie,:)).*cos(xg(ie,:)-yg(ie,:)+t)).*w1(:,1)'.*w2(:,1)')*N(:,:)')'*J;
    %                                             4x1    +  ((   1x9  +   1x9  ) .*         1x9         )   .*    9x1'.*    9x1' *(4x9 )')'1x1    )

    
    
end
K = sparse(K); %create sparse matrix for optimization of performance
M = sparse(M); %create sparse matrix for optimization of performance



ABOVE CODE IS SAME AS BELOW BUT FASTER --------

%% Preallocatoin
K = zeros(nn,nn);           % preallocation global stiffness matrix
M = zeros(nn,nn);
F = zeros(nn,1);            % preallocation global force matrix
k = zeros(2,2);             % preallocation k tensor
Di = zeros(2,nl);            % preallocation for Derivative of shape for i
Dj = zeros(2,nl);            % preallocation for Derivative of shape for j
%% Evaluation and Assembly
for ie = 1:ne^2
    ke = zeros(nl);            % preallocate local stiffness matrix
    me = zeros(nl);
    fe = zeros(nl,1);          % preallocate local force matrix
    
    for ig = 1:ngp^2
        k(1,1) = yg(ie,ig);    % k tensor y value
        k(2,2) = xg(ie,ig);  % k tensor x value
        Di(1,:) = Dxi(:,ig);   % (1 by nl) = (nl by 1)
        Di(2,:) = Deta(:,ig);
        Dj(1,:) = Dxi(:,ig);
        Dj(2,:) = Deta(:,ig);
        ke(:,:) = ke(:,:) + w1(ig)*w2(ig)*(Di' * k *  Dj); % local stiffness matrix

    end %ig
    %
    for il = 1:nl
        for jl = 1:nl
            me(il,jl) = me(il,jl) + N(jl,:).*N(il,:)*(w1(:,1).*w2(:,1))*J;
        end
    end
    
    fe(:,1) = fe(:,1) + ((((xg(ie,:)+yg(ie,:)).*cos(xg(ie,:)-yg(ie,:)+t)).*w1(:,1)'.*w2(:,1)')*N(:,:)')'*J;
    %  4x1     4x1    +  ((   1x9  +   1x9  ) .*         1x9         )   .*    9x1'.*    9x1' *(4x9 )')'1x1    )
    
    % Global Assembly
    % K(ENM(il,ie),ENM(jl,ie))=K(ENM(il,ie),ENM(jl,ie))+ke(il,jl);
    K(ENM(:,ie),ENM(:,ie))  = K(ENM(:,ie),ENM(:,ie))+ke(:,:);
    M(ENM(:,ie),ENM(:,ie))  = M(ENM(:,ie),ENM(:,ie))+me(:,:);
    % F(ENM(il,ie),1)=F(ENM(il,ie),1)+fe(il);
    F(ENM(:,ie),1)          = F(ENM(:,ie),1)+fe(:,1);
    
    
    
end
K = sparse(K); %create sparse matrix for optimization of performance
M = sparse(M); %create sparse matrix for optimization of performance
end


------------------------------------------------------------------------------------------
f1_BuildF
Decription:
1. calculate local force matrix 
2. assign them to global force matrix

Comment: 

for ie = 1:ne^2
    F(ENM(:,ie),1)  = F(ENM(:,ie),1)+((((xg(ie,:)+yg(ie,:)).*cos(xg(ie,:)-yg(ie,:)+t)-sin(xg(ie,:)-yg(ie,:)+t)).*w1(:,1)'.*w2(:,1)')*N(:,:)')'*J;
end


ABOVE CODE IS SAME AS BELOW BUT FASTER --------

for ie = 1:ne^2
    fe = zeros(nl,1);
    fe(:,1) = fe(:,1) + ((((xg(ie,:)+yg(ie,:)).*cos(xg(ie,:)-yg(ie,:)+t)-sin(xg(ie,:)-yg(ie,:)+t)).*w1(:,1)'.*w2(:,1)')*N(:,:)')'*J;    
    
    % Global Assembly
    % F(ENM(il,ie),1)=F(ENM(il,ie),1)+fe(il);
    F(ENM(:,ie),1)          = F(ENM(:,ie),1)+fe(:,1);        
end

------------------------------------------------------------------------------------------
f2_ApplyBC
Decription:
Apply Neuman (Natural) and Drichlet (Essential) Boundary Condition 

Comment: 

Drichlet Boundary Condition 

%% Drichlet BC 
% gD = cos(x) for (:,0)=bottom
% gD = cos(x-1) for (:,1) = top
K(BM(:,1),:)           = 0; %bottom node zeroing 
K(BM(:,1),BM(:,1)) = K(BM(:,1),BM(:,1))+eye(size(BM,1),size(BM,1)); %Modifying only the diagonal term by adding the identity matrix of the size of boundary matrix row by boundary matrix row
F(BM(:,1),1)           = cos(CM(BM(:,1),1)); %bottom node forcing term 
K(BM(:,3),:)           = 0; %top node zeroing 
K(BM(:,3),BM(:,3)) = K(BM(:,3),BM(:,3))+eye(size(BM,3),size(BM,3)); %Modifying only the diagonal term by adding the identity matrix of the size of boundary matrix row by boundary matrix row
F(BM(:,3),1)           = cos(CM(BM(:,3),1)-1); %top node forcing term 


ABOVE CODE IS SAME AS BELOW BUT FASTER --------

for ie = 1:ne+1    
    % Bottom node 1 
    K(BM(ie,1),:)          = 0;
    K(BM(ie,1),BM(ie,1))    = 1; %Bottom node 1
    F(BM(ie,1),1)           = cos(CM(BM(ie,1),1));
    
    % Bottom node 2 
    K(BM(ie+1,1),:)         = 0;
    K(BM(ie+1,1),BM(ie+1,1)) = 1; %Bottom node 2 
    F(BM(ie+1,1),1)          = cos(CM(BM(ie+1,1),1));
    
    % Top node 1 
    K(BM(ie,3),:)            = 0;
    K(BM(ie,3),BM(ie,3))     = 1; %Top node 1 
    F(BM(ie,3),1)            = cos(CM(BM(ie,3),1)-1);
    
    % Top node 2 
    K(BM(ie+1,3),:)          = 0;
    K(BM(ie+1,3),BM(ie+1,3)) = 1; %Top node 2 
    F(BM(ie+1,3),1)          = cos(CM(BM(ie+1,3),1)-1);
end


---------------




yg1= yg(1:ne:ne^2,1:ngp);%mapped y function for 1D 
 F(BM(1:end-1,2))          = F(BM(1:end-1,2))-((sin(1-yg1(:,:)))*(w(1,:)'.*Neta(1:ngp,1)))*0.5*1/ne;  
 F(BM(2:end,2))            = F(BM(2:end,2))  -((sin(1-yg1(:,:)))*(w(1,:)'.*Neta(1:ngp,2)))*0.5*1/ne; 

ABOVE CODE IS SAME AS BELOW BUT FASTER --------

iyg = 1;
for ie = 1:ne 
     F(BM(ie,2))    = F(BM(ie,2))-0.5*1/ne*(sin(1-yg(iyg,1:ngp))*(w(1,:)'.*Neta(1:ngp,1))); 
                                %     J *     1x3             *  (1x3)'.*3x1   
     F(BM(ie+1,2))  = F(BM(ie+1,2))-0.5*1/ne*(sin(1-yg(iyg,1:ngp))*(w(1,:)'.*Neta(1:ngp,2)));
                                  %   J   *     1x3             *  (1x3)'.*3x1   
iyg=iyg+ne;  %counter for yg, because I reused the 2D mapped y value, so I have too many and I don't need all of them. 
end 


------------------------------------------------------------------------------------------
f2_ApplyBC2
Decription:
Apply Neuman (Natural) and Drichlet (Essential) Boundary Condition only to force matrix 

Comment: 

yg1= yg(1:ne:ne^2,1:ngp);%mapped y function for 1D 
 F(BM(1:end-1,2))          = F(BM(1:end-1,2))-((sin(1-yg1(:,:)))*(w(1,:)'.*Neta(1:ngp,1)))*0.5*1/ne;  
 F(BM(2:end,2))            = F(BM(2:end,2))  -((sin(1-yg1(:,:)))*(w(1,:)'.*Neta(1:ngp,2)))*0.5*1/ne; 

ABOVE CODE IS SAME AS BELOW BUT FASTER --------

iyg = 1;
for ie = 1:ne 
     F(BM(ie,2))    = F(BM(ie,2))-0.5*1/ne*(sin(1-yg(iyg,1:ngp))*(w(1,:)'.*Neta(1:ngp,1))); 
                                %     J *     1x3             *  (1x3)'.*3x1   
     F(BM(ie+1,2))  = F(BM(ie+1,2))-0.5*1/ne*(sin(1-yg(iyg,1:ngp))*(w(1,:)'.*Neta(1:ngp,2)));
                                  %   J   *     1x3             *  (1x3)'.*3x1   
iyg=iyg+ne;  %counter for yg, because I reused the 2D mapped y value, so I have too many and I don't need all of them. 
end 



------------------------------------------------------------------------------------------
f3_PlotResult
Decription:
Plot both approximate and theoretical solution to 3D graph using patch and view function 

Comment: 



------------------------------------------------------------------------------------------
f3_AnimateResult
Decription:
Animate the solution using get frame function

Comment: 



------------------------------------------------------------------------------------------
CalcError
Decription:
Calculate the L2 error 

Comment: 




------------------------------------------------------------------------------------------
Script_As_Function 
Description:
take number of element and calculate the entire matlab script, useful for calculating L2 for decreasing number of element


------------------------------------------------------------------------------------------
makeTable
Description:
This code run the main code with increasing number of element to plot the error graph and create table of number of element and L2 error 


------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Line clarification                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Explanation of jacobian calculation 
J   = det([x(:,1)'*Dxi(:,1) x(:,1)'*Deta(:,1); y(:,1)'*Dxi(:,1) y(:,1)'*Deta(:,1)]);

Since the det of the value do not change based on the x and y coordinates, 
the jacobian becomes constant for any element for any gaussian points, so we only need 1 
The above line is doing the same thing as function below 

function J = f0_jacob (ne,ngp,x,y,Dxi,Deta) 
delx_xi = x'*Dxi(:,ig) ;
delx_eta = x'*Deta(:,ig); 
dely_xi = y'*Dxi(:,ig);
dely_eta = y'*Deta(:,ig);
J_matrix = [delx_xi delx_eta; dely_xi dely_eta]
J = det(J_matrix);

If Jacobian changes, then below code will be helpful to store Jacobian matrix 
J = zeros(ne^2,ngp^2); 
for ie = 1:ne^2
    for ig = 1:ngp^2 
        J(ie,ig) = det([x(:,ie)'*Dxi(:,ig) x(:,ie)'*Deta(:,ig); y(:,ie)'*Dxi(:,ig) y(:,ie)'*Deta(:,ig)]);
    end
end 
end 

    