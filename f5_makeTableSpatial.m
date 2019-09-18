function f5_makeTableSpatial(neM,ngp,timestep,alpha)
% This code run the main code with increasing number of element to plot the
% error graph and create table of number of element and L2 error 
fprintf('....Generating L2 error for number of element from %d to %d with 2^n increment for  %d timestep, \n %d number of gaussian points, and alpha = %d \n\n',neM(1), neM(end),1/timestep,ngp,alpha);

nnM = zeros(size(neM,2),1);
tM = zeros(size(neM,2),1);
cycleM = zeros(size(neM,2),1);
L2M = zeros(1,size(neM,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic; % begin measuring time
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate L2 error for increasing step size 
for i = 1:size(neM,2)
[L2,nn] = f5_Script_As_Function(neM(i),ngp,timestep,alpha);
L2M(i) = L2; 
nnM(i) = nn;
tM(i) = 1/timestep; 
cycleM(i) = i;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t2 = toc; %end measuring time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
log_L2 = log(L2M); 
log_h = log(neM); 

%% Plot log log graph 
figure('Name', 'loglog plot of step size vs. L2 error');  
plot(log_h,log_L2);
legend({'L2'});
xlabel('log stepsize'); 
ylabel('log error'); 
title(['timestep: ' num2str(timestep) ', number of elements: ' num2str(neM)]);
xlim auto;
ylim auto;

%% calculate convergence rate 
coefL2 = polyfit(log_h,log_L2,1);
slopeL2 = -coefL2(1);

cycle = cycleM;
time_step = tM;
mesh_size =neM'; 
cell_number =(neM.^2)';
dof_number=nnM; 
L2_error = L2M'; 

T = table(cycle,time_step,mesh_size,cell_number,dof_number,L2_error); 
disp(T); 
fprintf('\n The convergence rate is %d \n',slopeL2);
%movegui('center');
fprintf(' Elapsed time for %d by %d to %d by %d element for %d timestep is:             %d seconds\n',neM(1),neM(1),neM(end),neM(end),timestep,t2);

end 