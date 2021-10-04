function p_opt = power_control_SAR(d_i,theta_i, M, p_T, S, delta, beam)
%Per-beam power control mechanism for the RAB scheme proposed in [REF].

%% Input parameters:
% d_i    -> Distance between the i-th EH device and the PB
% theta_i-> azimuth angle relative to the boresight of the transmitting
%           ULA at which the i-th EH device is
% M      -> Number of transmit antennas of the PB
% p_T    -> total transmit power budget of the PB
% S      -> MxM SAR matrix as defined after [REF, eq.20]
% delta  -> SAR contraint as defined before [REF, eq.21]
% beam   -> vector with the indeces of the beams that are SAR-constrained

%% Output parameters:
% p_opt  -> 1xN optimum power allocation vector

%% Main code
% Computation of [REF, eq.21] (z=1)
s_z = 0;
for i=1:M
    for j=1:M
        if mod(i+j,2)==0
            s_z = s_z + S(i,j);
        else
            s_z = s_z - S(i,j);
        end
    end
end

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
    a_i(:,j2)=d_i.^(-2).*Pi/M;  %iteratively creates a matrix with all a_ij
                                %entries using [REF, eq.16]. Here the path
                                %loss exponent is set to 2.
end

% LP solution of problem P in [REF, eq.17] + SAR constraints [REF, eq.22]. 
% This requires CVX solver (http://cvxr.com/cvx/)
% Alternatively one could directly use MatLab LP solver
cvx_begin quiet
    variable p(M) 
    variable xi
    minimize( - xi )    
    for i=1:length(d_i)           
        xi - a_i(i,:)*p <= 0;
    end
    ones(1,M)*p - p_T <= 0; 
    p(beam)*s_z/M <= delta;  %SAR constraint per beam [REF, eq.22]
    p>=0;
cvx_end

p_opt = p;  %Exact optimal power allocation

%% References:
%[REF]    - O. L. A. López, H. Alves, S. Montejo-Sánchez, R. D. Souza and 
%           M. Latva-aho, "CSI-free Rotary Antenna Beamforming for Massive
%           RF Wireless Energy Transfer," in IEEE Internet of Things 
%           Journal, doi: 10.1109/JIOT.2021.3107222.


