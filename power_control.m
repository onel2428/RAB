function p_opt = power_control(d_i,theta_i, M, p_T)
%Per-beam power control mechanism for the RAB scheme proposed in [REF].

%% Input parameters:
% d_i    -> Distance between the i-th EH device and the PB
% theta_i-> azimuth angle relative to the boresight of the transmitting
%           ULA at which the i-th EH device is
% M      -> Number of transmit antennas of the PB
% p_T    -> total transmit power budget of the PB

%% Output parameters:
% p_opt  -> 2xN power allocation vector. First/second row corresponds to
%           the exact/approximate solution


   
%% Main code

% This 3 level loop computes [REF, eq.14]
pii =@(aa,jj,ll) cos(mod(jj-ll,2)*pi+(ll-jj)*pi*sin(aa));
for j2 = 1:M
    Pi = ones(length(d_i),1);
    % this 2 level loop computes [REF, eq.11] for each input angle
    for j=1:M
        for l=j+1:M
            Pi = Pi + (2/M)*pii(theta_i+j2*pi/M,j,l);
        end
    end
    a_i(:,j2)=d_i.^(-2).*Pi/M; %iteratively creates a matrix with all a_ij
                               %entries using [REF, eq.16]. Here the path
                               %loss exponent is set to 2.
end

%% Exact solution

% LP solution of problem P in [REF, eq.17]. 
% This requires CVX solver (http://cvxr.com/cvx/)
% Alternatively one could directly use MatLab LP solver
cvx_begin quiet
    variable p(M) 
    variable xi
    minimize( - xi )    
    for i=1:length(d_i)           
        xi-a_i(i,:)*p <= 0;
    end
    ones(1,M)*p-p_T <= 0;    
    p >= 0;
cvx_end
p_opt(1,:) = p; % Exact optimal power allocation

%% Approximate solution

p = [];
% Next code computes the approximate power allocation as described in
% [REF], specifically around [REF, eqs.18,19]
theta_i(theta_i>pi/2 & theta_i< 3*pi/2 )= theta_i(theta_i>pi/2 & theta_i< 3*pi/2 ) - pi;
theta_i(theta_i >= 3*pi/2 )= theta_i(theta_i>= 3*pi/2 ) - 2*pi;
for j=1:M
   a1 = theta_i<=pi/2-j*pi/M+pi/(2*M);
   a2 = theta_i>=pi/2-j*pi/M-pi/(2*M);
   a = a1.*a2;
   if sum(a)<1
       p(j) = 0;
   else
       p(j)=max(d_i(a==1))^2;
   end
end

% Approximate power allocation [REF, eq.19]
p = p_T*p/sum(p);
p_opt(2,:) = p; 

%% References:
%[REF]    - O. L. A. López, H. Alves, S. Montejo-Sánchez, R. D. Souza and 
%           M. Latva-aho, "CSI-free Rotary Antenna Beamforming for Massive
%           RF Wireless Energy Transfer," in IEEE Internet of Things 
%           Journal, doi: 10.1109/JIOT.2021.3107222.