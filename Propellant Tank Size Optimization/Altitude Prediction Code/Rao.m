%% Rao Function
% Rao Parabolic Approximation Method for sizing nozzles
function [xx, coeffs, R_x] = Rao(R_t, epsilon, Lf_ratio,theta_E,theta_N,theta_C)
% returns coefficients for parabola, and the radius vector
% note: using fixed exit and entrance angles, need to get numeric
% data for plots so these can be chosen automatically


global n_rao

n = n_rao;
% *************** Don't change anything in this section ******************
R_dt = 0.382*R_t;
R_ut = 1.5*R_t;
R_cu = R_ut + R_t;
R_cd = R_dt + R_t;

% geometry (see notes for figure)
L_n = Lf_ratio*(R_t*(sqrt(epsilon)-1)+...
    R_dt*(sec(15*pi/180)-1))/tan(15*pi/180);

% solving for parabola coeffs
x_N = tan(theta_N)*R_dt/(sqrt(1+tan(theta_N))^2);
R_N = R_cd - sqrt(R_dt^2-x_N^2);
c = (tan(theta_N)-tan(theta_E))/(2*(x_N-L_n));
b = tan(theta_N)-2*c*x_N;
a = R_N-x_N*(b+c*x_N);

coeffs = [a;b;c];

x_co = tan(theta_C)*R_ut/sqrt(1+tan(theta_C)^2);

xx = linspace(-1.5*x_co,L_n,n);
for i = 1:n   
    if(xx(i) >= x_N)
        R_x(i) = a + b*xx(i) + c*xx(i)^2;        
    elseif (xx(i) < x_N && xx(i) >= 0)
        R_x(i) = R_cd-sqrt(R_dt^2-xx(i)^2);       
    else
        R_x(i) = R_cu-sqrt(R_ut^2-xx(i)^2);
    end
end


