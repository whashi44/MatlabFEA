function [CM,ENM,BM] = f0_ID(ne,nl,h,nn)
%% Coordinate Matrix
%{
This is a global matrix, the counting system work 1, 2, 3 for each
element going from left to right, and then if you reach the right edge
then you will start from the left and repeat. 
%}
CM  = zeros(nn,2);  % Coordinate Matrix column 1 = x, column 2 = y.
pts = 0:1/ne:1;     % Array of 1D coordinates 
cj  = 1;            % counter for for loop to assign coordinates correctly 

for ie = 1:ne+1 
    for jl = 1:ne+1
        CM(cj,1) = pts(jl); % store x coordinate
        CM(cj,2) = pts(ie); % store y coordinate
        cj       = cj+1;    % increment counter 
    end
end

%% Elemental Node Matrix 
%{
There are only 4 nodes per element, and if you look at the ENM matrix, the
patter is you start with a certain number, and then you increment by 1 for
each element, and for each ne+1 column, you increment by 2. That means, if
you create a matrix first with starting number and increment of 1, until
specific ending number, and then remove the column for each ne+1, then you
get ENM matrix. 
 %}
j           = 1;                        % for readability 
ENM         = zeros(nl,ne^2+ne-1);      %create an array 
ENM(1,:)    = (j:(ne+1)^2-(ne+2));      %node 1 start from 1, end with that value  
ENM(2,:)    = (j+1:(ne+1)^2-(ne+1));    %node 2 start from 2, 
ENM(3,:)    = (j+ne+2:(ne+1)^2); 
ENM(4,:)    = (j+ne+1:(ne+1)^2-1);  
ENM(:,ne+1:ne+1:end)=[];                %remove each ne+1 column starting from ne+1 
    
%% Boundary Matrix
%{
Assign boundary matrix, row represents number of total nodes in 1D 
and column represents 1:bottom, 2:right, 3:top, 4:left 
h is domain, [0 1] so h(1) = 0, h(2) = 1
%}
BM      = zeros(ne+1,4);       % Preallocatin
BM(:,1) = find(CM(:,2)==h(1)); % (:,0), bottom
BM(:,2) = find(CM(:,1)==h(2)); % (1,:), right
BM(:,3) = find(CM(:,2)==h(2)); % (:,1), top
BM(:,4) = find(CM(:,1)==h(1)); % (0,:), left

end

