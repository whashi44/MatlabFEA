function f4_PlotResult(ne,nl,x,y,ENM,T,tend,timestep) 
% This function plot the approximate and theoretical solutoins 

%% creating points 
T_approx    = zeros(nl,ne^2);  %preallocate approximate solution 
T_theor     = zeros(nl,ne^2);  %preallocate theoretical solution 
for ie = 1:ne^2 %for every element in 2D 
    T_approx(:,ie)  = T(ENM(:,ie));         %grab the approximate solution from T table using    
    T_theor(:,ie)   = cos(x(:,ie)-y(:,ie)+tend); %evaluate the theoretical solution 
end 

%% graphing 
figure('Name', 'Analytical vs. Theoretical Solution');  
hold on % for creating multiple graphs in same figure 
%view(2);% for viewing the patch graph in 3D 
view(3);
%patch(x,y,T_approx,T_approx, 'LineStyle', 'none');
patch(x,y,T_approx,T_approx,'FaceColor', 'none'); %graphing the approximate solution 
patch(x,y,T_theor,T_theor, 'LineStyle', 'none');   %graphing the theoretical solution 
hold off  

xlabel('x coordinates');
ylabel('y coordinates');
zlabel('Temperature');
%legend('Analytical','Theoretical','Location','best'); 
title(['Grid size: ' num2str(ne) 'x' num2str(ne) '  Total number of element: ' num2str(ne^2) ' Timestep: ' num2str(timestep)]);
colorbar;
%movegui('east');
fprintf('\n Created %d by %d by %d mesh grid. The total number of element is %d \n',ne,ne,timestep,ne^2);

end