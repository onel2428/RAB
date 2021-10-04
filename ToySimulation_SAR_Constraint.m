%This script allows generating a figure similar to [REF1, Fig.13]
%With the output data here, one can also generate a figure similar to [REF1, Fig.14]

%% System parameters
M = 4;         %Number of PB antennas
p_T = 3;       %Number of transmit antennas at the PB
N = 32;        %Number of EH devices
samples = 250; %Number of Monte Carlo samples (number of different deployment realizations)
beam = [1 2];  %Beams 1 and 2 under SAR constraints
% SAR matrix [REF1, eq.23]:
S = [1.6 -1.2*1i -0.42 0; 1.2*1i 1.6 -1.2*1i -0.42; -0.42 1.2*1i 1.6 -1.2*1i; 0 -0.42 1.2*1i 1.6];

%% Main section of the code script
DELTA = 0.25:0.25:2;
Full_CSI = zeros(1,length(DELTA));
RAB_OptPC = zeros(1,length(DELTA));
for z=1:length(DELTA)  %results as a function of the SAR constraint
    
    delta = DELTA(z);
    h = [];
    h_rot = [];
    full_csi = zeros(1,samples);
    rab_optPC = zeros (N,samples);
    
    for j=1:samples   
        [delta j]
        rand('seed',j)
        
        % Devices are randomly and uniformly distributed in a 10 m-radius
        % area around the PB
        radius = 10*sqrt(rand(N,1));
        theta = 2*pi*rand(N,1);
        beta = 1e-4*radius.^(-2);    %path loss model
        
        %LOS channel generation for each device
        h = zeros(N,M);
        h_rot = zeros(N,M,M);
        for i=1:N
            phi = -(0:(M-1))*pi*sin(theta(i));  %[REF1, eq.2]     
            h(i,:) = exp(1i*pi/4)*exp(1i*phi);  %[REF1, eq.1 with varphi_0=pi/4]
                
            %For the rotary mechanism, the LOS channel changes for each beam
            %direction m
            for m=1:M
                phi_rot = -(0:(M-1))*pi*sin(theta(i)+pi*m/M)+mod(0:(M-1),2)*pi;
                h_rot(i,:,m) = exp(1i*pi/4)*exp(1i*phi_rot);                 
            end
        end
             
        
        %Full-CSI (SDP-based) [REF2]
        [~, E] = CSIBeamf_SDP_SAR(M,N,h,beta,S,delta);
        full_csi(j) = p_T*E;
       
        % Optimal power control [REF1, solution to P eq.17]
        popt = power_control_SAR(radius, theta, M, p_T, S, delta, beam);

        %RAB performance for each EH device (under optimal power control)
        for i=1:N                      
            E=beta(i).*abs(sum(squeeze(h_rot(i,:,:))./sqrt(M))).^2;                        
            rab_optPC(i,j) = sum(popt(1,:).*E);            
        end
    end
    
    % Results in dBm
    Full_CSI(z) = 10*log10(nanmean(full_csi)*1e3);   
    RAB_OptPC(z) = 10*log10(nanmean(min(rab_optPC,[],1))*1e3);     
   
end

%% Plot of the results
figure
set(gcf, 'Units', 'centimeters'); 
axesFontSize = 14;
legendFontSize = 12;
afFigurePosition = [2 7 14 10]; 
set(gcf, 'Position', afFigurePosition,'PaperSize',[14 10],'PaperPositionMode','auto'); % [left bottom width height], setting printing properties 

h1=plot(DELTA,Full_CSI,'--pr','LineWidth',1.5,'MarkerSize',7,'MarkerFaceColor','w'); hold on
h2=plot(DELTA,RAB_OptPC,'--x','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','Color',[255,140,0]/255); hold on

grid on
box on

set(gca,'FontSize',12)
xlabel('SAR constraint, $\delta$ [W/Kg]','Interpreter', 'latex','fontsize',axesFontSize)
ylabel('average worst-case RF energy (dBm)','Interpreter', 'latex','fontsize',axesFontSize)

hl=legend([h1 h2],'full-CSI', '$\mathrm{RAB}$ (optimum power control - (17))' ); 
set(hl,'interpreter','latex','FontSize',legendFontSize);

xlim([min(DELTA) max(DELTA)])
xticks(DELTA)
yticks(-32:2:-20)
ylim([-32 -20])

%% References:
%[REF1]   - O. L. A. López, H. Alves, S. Montejo-Sánchez, R. D. Souza and 
%           M. Latva-aho, "CSI-free Rotary Antenna Beamforming for Massive
%           RF Wireless Energy Transfer," in IEEE Internet of Things 
%           Journal, doi: 10.1109/JIOT.2021.3107222.
%[REF2]   - O. L. A. López, F. A. Monteiro, H. Alves, R. Zhang and M. Latva-Aho,
%           "A Low-Complexity Beamforming Design for Multiuser Wireless Energy 
%           Transfer," in IEEE Wireless Communications Letters, vol. 10,
%           no. 1, pp. 58-62, Jan. 2021, doi: 10.1109/LWC.2020.3020576.