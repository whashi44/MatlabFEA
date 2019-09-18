function [xi,eta,w12,w] = f0_gauss(ngp)
%This function grab the gauss table from 1D and transform into 2D
% Input: ngp (number of gauss point)
% Output: xi,eta (abscissae) and w12 (weight in 2D), and w(weight in 1D)

%% 1D gaussina quadarture
switch ngp %switch cases by number of gaussian points
    case 1  %if number of gauss point is 1
        gp  = 0; %abscissae is 0
        w   = 2; %weight is 2
    case 2   % if number of gauss point is 2
        gp  = [-1 +1]/sqrt(3);
        w   = [1 1];
    case 3
        gp  = [-1 0 +1]*sqrt(3/5);
        w   = [5 8 5]/9;
    case 4
        gp  = [-sqrt(3/7+2/7*sqrt(6/5)) +sqrt(3/7+2/7*sqrt(6/5)) -sqrt(3/7-2/7*sqrt(6/5)) +sqrt(3/7-2/7*sqrt(6/5))];
        w   = [18-sqrt(30) 18-sqrt(30) 18+sqrt(30) 18+sqrt(30)]/36;
    case 5
        gp  = [-sqrt(5+2*sqrt(10/7)) sqrt(5+2*sqrt(10/7)) 0 -sqrt(5-2*sqrt(10/7)) sqrt(5-2*sqrt(10/7)) ]/3;
        w   = [322-13*sqrt(70) 322-13*sqrt(70) 128/225*900 322+13*sqrt(70) 322+13*sqrt(70) ]/900;
    case 6
        gp  = [-0.9324695142 0.9324695142 -0.6612093865 0.6612093865 -0.2386191861 0.2386191861];
        w   = [0.1713244924 0.1713244924 0.3607615730 0.3607615730 0.4679139346 0.4679139346];
end

%% 2D gaussian quadrature
%{
We want eta and xi interchangeably
For example for ngp = 3, use same xi for first 3, and use 3 different eta for 3. 
Then you will get the xi and eta combination that we want 
%} 
% Preallocation
xi  = zeros(ngp^2,1);
eta = zeros(ngp^2,1);
w1  = zeros(ngp^2,1);
w2  = zeros(ngp^2,1);
cj  = 1;

%eta, xi and weight allocation
for ig =1:ngp
    for jg = 1:ngp
        xi(cj)  = gp(ig);  
        eta(cj) = gp(jg);
        w1(cj)  = w(ig);
        w2(cj)  = w(jg);
        cj      = cj+1;
    end %jg
end %ig
w12 = w1(:,1).*w2(:,1);

end

