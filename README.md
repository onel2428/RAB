# RAB
Rotary Antenna Beamforming for Wireless Energy Transfer

This compilation of MatLab scripts is used to generate the results illustrated in [REF1].
The scripts description is as follows.

full-CSI WET precoding:
  - CSIBeamf_SDP (solution of P2 eq.4 in [REF2])
  - CSIBeamf_SDP_SAR (solution of [REF2, P2 eq.4] + SAR constraints defined in [REF1, eq.20])

rotary energy beamforming (RAB):
  - power_control (optimal/suboptimal power allocation [REF1, eq.17/eq.19])
  - power_control_SAR (optimal power allocation with SAR constraints [REF1, eq.17] + SAR constraints [REF1, eq.22])

Scripts for figures generation:
  - ToySimulation_Coverage_Area.m (It allows generating a figure similar to [REF1, Fig.7]: Area coverage for different average RF energy requirements. 
  - ToySimulation_Num_Devices.m (It allows generating a figure similar to [REF1, Fig.8]: Average worst case RF energy at the set of IoT devices.
  - ToySimulation_SAR_Constraint.m (It allows generating a figure similar to [REF1, Fig.10]:   Average worst case RF energy at a set of 32 EH devices as a function of SAR constraint. 

References:

[REF1] - O. L. A. López, H. Alves, S. Montejo-Sánchez, R. D. Souza and M. Latva-aho, "CSI-free Rotary Antenna Beamforming for Massive RF Wireless Energy Transfer," in IEEE Internet of Things Journal, doi: 10.1109/JIOT.2021.3107222.

[REF2] - O. L. A. López, F. A. Monteiro, H. Alves, R. Zhang and M. Latva-aho, "A Low-Complexity Beamforming Design for Multiuser Wireless Energy Transfer," in IEEE Wireless Communications Letters, doi: 10.1109/LWC.2020.3020576.
