function t3 = f5_makeTableTemporal(dtM,ngp,ne,alpha)
% This code run the main code with increasing number of element to plot the
% error graph and create table of number of element and L2 error 
fprintf('....Generating L2 error for timestep %d to %d with 2^n increment for %d number of elements, \n %d number of gaussian points and alpha = %d \n\n',dtM(1), dtM(end),ne,ngp,alpha);
L2M = zeros(1,size(dtM,2));
nnM = zeros(1,size(dtM,2));
cycleM = zeros(1,size(dtM,2));
neM = zeros(1,size(dtM,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic; % begin measuring time
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate L2 error for increasing step size 
for i = 1:size(neM,2)
[L2,nn] = f5_Script_As_Function(ne,ngp,dtM(i),alpha);
L2M(i) = L2; 
nnM(i) = nn;
neM(i) = ne;
cycleM(i) = i;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t3 = toc; %end measuring time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
log_L2 = log(L2M); 
log_h = log(dtM); 

%% Plot log log graph 
figure('Name', 'loglog plot of step size vs. L2 error');  
plot(log_h,log_L2);
legend({'L2'});
xlabel('log stepsize'); 
ylabel('log error'); 
title(['time steps: ' num2str(dtM) ', number of element: ' num2str(ne)]);
xlim auto;
ylim auto;

%% calculate convergence rate 
coefL2 = polyfit(log_h,log_L2,1);
slopeL2 = -coefL2(1);

cycle = cycleM';
time_step = (1./dtM)';
mesh_size =neM'; 
cell_number =(neM.^2)'; 
dof_number=nnM'; 
L2_error = L2M'; 

T = table(cycle,time_step,mesh_size,cell_number,dof_number,L2_error); 
disp(T); 
fprintf('\n The convergence rate is %d \n',slopeL2);
%movegui('center');
fprintf(' Elapsed time for %d to %d  timestep for %d element is:             %d seconds\n',dtM(1),dtM(end),ne,t3);

end 