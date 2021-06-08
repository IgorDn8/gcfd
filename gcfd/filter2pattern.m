function [pattern] = filter2pattern (filter, elec)
% Convert a spatial filter to a spatial pattern based on formula (19) from 
% Nikulin, Nolte, Curio "Cross-frequency decomposition: A novel technique for studying interactions between neuronal oscillations with different frequencies", Clinical Neurophysiology, 2012
% who, in turn, cite
% Parra LC, Spence CD, Gerson AD, Sajda P. "Recipes for the linear analysis of EEG" Neuroimage 2005;28:326?41
% 
% INPUT: 
%     filter - spatial filter (demixing coefficients), column vector of size N
%     elec   - signals for all electrodes, matrix of size TxN
%
% OUTPUT:
%     pattern - spatial pattern (mixing coefficients), column vector of size N
%
% References:
%
% Volk, D., Dubinin, I., Myasnikova, A., Gutkin, B., & Nikulin, V. V. (2018). 
% Generalized Cross-Frequency Decomposition: A Method for the Extraction of Neuronal Components Coupled at Different Frequencies. 
% Frontiers in neuroinformatics, 12.
%

V = elec * filter;
V = detrend(V,'constant');
%pattern = (elec' * V) / (V' * V);
pattern = (elec' * V) / sqrt(V' * V);

end