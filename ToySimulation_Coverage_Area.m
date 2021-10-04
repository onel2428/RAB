%This script allows generating a figure similar to [REF1, Fig.7]
%The performance curves corresponding to the CSI-free schemes are not
%displayed with this code. 

%% System parameters
M = 4;        %Number of PB antennas
R = 10;       %EH devices are deployed in a RxR m^2 area
p_T = 3;      %Total transmit power of the PB
% PB is at coordinate (0,0)
x0 = 0;
y0 = 0;

%% Main section of the code script
% (xi,yi) are the coordiantes of the point in the area where the RF energy 
% is measured via simulation
E_RAB = zeros(length(-R:0.2:R));
row = 1;
for xi = -R:0.1:R    
    col = 1;    
    for yi = -R:0.1:R                   
        d = sqrt((x0-xi).^2+(y0-yi).^2);            % distance of point (xi,yi) with respect to origin (0,0) where PB is        
        beta = max(d,1).^(-2)*10^(-4);            % path loss   
        theta_i = atan((y0-yi)./(x0-xi));           % azimuth angle relative to the boresight of the transmitting ULA at which the i-th EH device is 
               
        %Simulation of the RF energy at a point (xi,yi) under LOS.
        %The precoder is [REF1,eq.5]
        channel_plus_precod = zeros(1,M);
        for j=1:M
            channel_plus_precod(j,:) = exp(-1i*pi*[0:(M-1)]*sin(pi/2-theta_i+pi*(j-1)/M)+1i*mod(0:M-1,2)*pi);
        end          
        E_RAB(row,col) = p_T*mean(beta.*abs(sum((channel_plus_precod')^2)/M));
        col = col + 1;
    end
    row = row + 1;
end


% CDF data points generated with histcount MatLab function
edgesdB = (-28:-4)-30;         %power values in dB
edges = [0 10.^(edgesdB/10)];  %power values in linear scale
[values,~] = histcounts(E_RAB(:),edges, 'Normalization', 'cdf');
edges=10*log10(edges(2:end));

% CDF curve
figure
set(gcf, 'Units', 'centimeters'); 
axesFontSize = 14;
legendFontSize = 12;
afFigurePosition = [2 7 14 10]; 
set(gcf, 'Position', afFigurePosition,'PaperSize',[14 10],'PaperPositionMode','auto'); 

plot(edges+30,100*(1-values),'x-','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','Color',[255,140,0]/255);
xticks([-28:2:-4]);  % in dBm
ylim([0 100])

grid on
box on

set(gca,'FontSize',12)
xlabel('average RF energy (dBm)','Interpreter', 'latex','fontsize',axesFontSize)
ylabel('minimum coverage area ($\%$)','Interpreter', 'latex','fontsize',axesFontSize)
xlim([-28 -4])


%% References:
%[REF1]   - O. L. A. López, H. Alves, S. Montejo-Sánchez, R. D. Souza and 
%           M. Latva-aho, "CSI-free Rotary Antenna Beamforming for Massive
%           RF Wireless Energy Transfer," in IEEE Internet of Things 
%           Journal, doi: 10.1109/JIOT.2021.3107222.
%[REF2]   - O. L. A. López, F. A. Monteiro, H. Alves, R. Zhang and M. Latva-Aho,
%           "A Low-Complexity Beamforming Design for Multiuser Wireless Energy 
%           Transfer," in IEEE Wireless Communications Letters, vol. 10,
%           no. 1, pp. 58-62, Jan. 2021, doi: 10.1109/LWC.2020.3020576.