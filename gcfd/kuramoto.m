function [k] = kuramoto ( s1, s2, f1, f2 )
% Calculates "Kuramoto index", a measure for phase-phase synchronization of two signals
% (possibly of different but rationally related frequencies)
%
% as proposed in
%
% P. Tass, M. G. Rosenblum, J. Weule, J. Kurths, A. Pikovsky, J. Volkmann, A. Schnitzler, and H.-J. Freund, 
% Phys. Rev. Lett. 81, 3291, October 1998
% "Detection of n:m Phase Locking from Noisy Data: Application to Magnetoencephalography"
% 
% INPUT: 
%     s1,s2 - narrow-band signals of same length, Tx1
%     f1,f2 - integers proportional to their frequencies
%
% OUTPUT:
%     k - Kuramoto index
%
% References:
%
% Volk, D., Dubinin, I., Myasnikova, A., Gutkin, B., & Nikulin, V. V. (2018). 
% Generalized Cross-Frequency Decomposition: A Method for the Extraction of Neuronal Components Coupled at Different Frequencies. 
% Frontiers in neuroinformatics, 12.
%

c1 = hilbert(s1);
c2 = hilbert(s2);

phase_rel = f2*unwrap(angle(c1)) - f1*unwrap(angle(c2));

k = abs(mean(exp(1i.*phase_rel)));

end