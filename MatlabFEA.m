%{
 Please refer to the readMe file for the detailed explanation of this code
 Matlab matrix is always (row,column), use space or "," for row and ";" for column
 If you want to see the L2 error with increasing element number, please run makeTable.m file.
 Functions include:
 f0_gauss:          Calculate gaussian abscissae and weight
 f0_ID:             Calculate Coordinate Matrix(CM), Elemental Node Matrix(ENM), Boundary Matrix(BM)
 f0_mapXY:          Calculate x and y coordinate for each node for each element and map it for gaussian integration
 f0_shape:          Calculate shape function and derivative of shape function with respect to xi and eta
 f1_BuildF:         Evaluate force matrix only during time loop
 f1_BuildKFM:       Evaluate local stiffness (ke) and force matrix (fe) and mass matrix (me),
 assembly to create global stiffness (K),force (F) and mass matrix
 f2_ApplyBC:        Apply Neuman (Natural) and Drichlet (Essential) Boundary Condition to K and F
 f2_ApplyBC2:       Apply Boundary Condition to only F during time loop
 f3_TimeLoop:       Loop through time to calculate the final Temperature
 f4_PlotResult:     Plot the result using patch function
 f4_AnimateResult:  Animate the result using patch and getframe function
 f4_CalcError:      Calculate L2 error
 f5_Script_As_Function:Take ne,ngp,timestep as input and run the entire script. Useful for makeTable
 f5_makeTableSpacial: Calculate L2 for increasing element size, generate slope of that log L2 vs. log step size, and create table
 f5_makeTableTemporal: Calculate L2 for increasing time step size, generate slope for log L2 vs. log step size, and create table
%}
clear; clc;     % Initilize script, clear is beter than clear all in terms of performance
format compact; % Reduce line spacing for dense command window
format short;   % Show less decimals
close all;      % close all the figures
%% User input variables
% user input variables for asking user to input values
% ne        = input('Number of element in 1D? \n(ex.If you like 64 by 64 grid, type in "64"): ');
% timestep  = input('How many time steps? \n(ex. If you like to split the domain by 10, type in "10"): ');
% alpha     = input('which time discretization method? \n (forward Euler (explicit) = 0 Crank-Nicolson (Midpoint) = 0.5 backward Euler (Implicit) = 1) : ');
% ngp       = input('Number of gaussian point in 1D? \(ex. If you want 3 by 3 gaussian point, type in "3", up to 6): ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;        % begin measuring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% direct input variables
ne          = 64;    % number of element in 1D (ex. if you want 64 by 64, type in "ne = 64")
timestep    = 10;    % number of step in time (ex. 40 means you divide the time domain by 40 step sizes)
alpha       = 1;     % Indication of which time discretization method to use forward Euler (explicit) = 0; backward Euler (Implicit) = 1;, Crank-Nicolson (Midpoint) = 0.5
ngp         = 3;     % number of gaussian point in 1D (ex. if you want 3 by 3 = 9 gaussian points, type in "ngp =3")
%% makeTable variable
% neM   = [8 16 32 64];   %number of element matrix: used for the makeTableSpatial function
% dtM   = [10 20 40 80];  %number of time step matrix: used for the makeTableTemporal function
%% Constants
deg     = 1;        % degree of the shape function
nl      = 4;        % number of local nodes, quadrilateral = 4
h       = [0 1];    % spatial domain
tdomain = [0 1];    % temporal domain
%% Derived variables
nn      = ((ne*deg)+1)^2; % number of total nodes
dt      = 1/timestep;     % step size
tend    = tdomain(2);
t       = 0:dt:tend;
%% Problem Initialization
% Calculated using functions
[xi,eta,w12,w]          = f0_gauss(ngp);                %Obain abscissae and weight
[N,Dxi,Deta,Neta]       = f0_shape(ngp,nl,xi,eta) ;     %Obtain shape function and derivative of shape function
[CM,ENM,BM]             = f0_ID(ne,nl,h,nn);            %Obain Coordinate Matrix (CM), Elemental Node Matrix (ENM) and Boundary Matrix (BM)
[x, y, xg, yg]          = f0_mapXY(ne,ngp,nl,CM,ENM,N); %Obtain local x y coordinates, mapped xg yg coordinates

%Jacobian is element and gaussian point independent, so we only need one, check read me for the detail of this line
J       = det([x(:,1)'*Dxi(:,1) x(:,1)'*Deta(:,1); y(:,1)'*Dxi(:,1) y(:,1)'*Deta(:,1)]);
%% Build K,F,M Matrix
[K,F,M] = f1_BuildKFM(ne,ngp,nl,nn,xg,yg,N,Dxi,Deta,w12,ENM,J,t(1));
%% Apply Drichlet (First) and Neuman (Second) Boundary Condition
[K,F]   = f2_ApplyBC (ne,CM,BM,K,F,ngp,w,yg,Neta,t(1));
%% Calculate Temperature
[T, T_graph] = f3_TimeLoop(ne,ngp,nn,xg,yg,N,w12,w,Neta,CM,BM,ENM,J,t,K,F,M,dt,alpha,tend);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1 = toc;     % end measuring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Displaying performance
fprintf('\n Elapsed time for %d x %d x %d (ne by ne by timestep)is: %d seconds\n', ne,ne,timestep,t1);
%% Plot Results
f4_PlotResult(ne,nl,x,y,ENM,T,tend,timestep);
%% Animate Result
%gif = f4_AnimateResult(ne,nl,x,y,ENM,T_graph,tend,t);
%% Calculate error
L2 = f4_CalcError(ne,N,xg,yg,w12,ENM,T,J,tend);
fprintf('\n The L2 error for %d x %d x %d (ne by ne by timestep)is: %d \n',ne,ne,timestep,L2);
%% Output results into table format

tablei = input('\n\n\nWould you like to make table of spatial, temporal or both? \n Type "1" for spatial  \n Type "2" for temporal \n Type "3" for both \n Type "4" to exit  \n Type here:');
switch tablei
    case 1
neM = input('\nInput element number with increasing order \n PLEASE PUT SQUARE BRACKETS [], "[4 8 16 32]"  NOT "4 8 16 32":');
timestep = input('\nHow many time steps? \n(ex. If you like to split the domain by 10, type in "10"): ');
f5_makeTableSpatial(neM,ngp,timestep,alpha);
    case 2
dtM = input('\nInput timesteps with increasing order \n  PLEASE PUT SQUARE BRACKETS [], "[10 20 40 80]"  NOT "10 20 40 80": ');
ne = input('\nNumber of element in 1D? \n(ex.If you like 64 by 64 grid, type in "64"): ');
f5_makeTableTemporal(dtM,ngp,ne,alpha);
    case 3
neM = input('\nInput element numbers with increasing order \n Ex. PLEASE PUT SQUARE BRACKETS [], "[4 8 16 32]"  NOT "4 8 16 32": ');
timestep = input('\nHow many time steps? \n(ex. If you like to split the domain by 10, type in "10"): ');
f5_makeTableSpatial(neM,ngp,timestep,alpha);
dtM = input('\nInput timesteps with increasing order \n PLEASE PUT SQUARE BRACKETS [], "[10 20 40 80]"  NOT "10 20 40 80": ');
ne = input('\nNumber of element in 1D? \n(ex.If you like 64 by 64 grid, type in "64"): ');
f5_makeTableTemporal(dtM,ngp,ne,alpha);
    case 4
        disp('OK! See you then!');
    otherwise
        disp('Not valid number');
end



