function [profile] = ComputeProteinAffinity(DNA_seq)

% function uses read in DNA_seq to calculate affinity for protein
% using PWM method
%
% call PWM method for base-pair level binding constants

global K0

% Calculate affinity using PWM given below
% Currently set up for CTCF using the CTCF PWM:
% Change the filename at the end for other transcription factors
 
[K_pwm,K_rev,~,~] = PWM_affinity(DNA_seq,K0(1),1.5,'CTCF_PWM.txt');

profile = K_pwm + K_rev;

end
