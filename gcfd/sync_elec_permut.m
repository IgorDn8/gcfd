function [  perm_error, perm_sync, mismatch ] = sync_elec_permut( filter, elec, sig_2 , Fr_1, Fr_2, cut_num, permute_num, lsqnonlin_param, filter_seeds_num) 
%   This function performs permutations based on Maris, E., Schoffelen, J. M., and Fries, P. (2007). 
%   "Nonparametric statistical testing of coherence differences" Journal of Neuroscience Methods 163, 161–175.
% 
% INPUT: 
%     filter - filter for target band
%     elec   - recording in target band
%     sig_2  - signal from reference band 
%     Fr_1  - integer proportional to target band 
%     Fr_2  - integer proportional to reference band
%     cut_num  - number of cut parts for permutations, for better performance make samples/cut_num integer  
%     permute_num  - number of permutations to be performed
%     lsqnonlin_param, filter_seeds_num  - parameters for non-linear optimization
%
% OUTPUT:
%     perm_error - p-value for the signal
%     perm_sync  - PLV values for all permuted signals
%     mismatch   - All differences between permuted an original filters
%
% References:
%
% Volk, D., Dubinin, I., Myasnikova, A., Gutkin, B., & Nikulin, V. V. (2018). 
% Generalized Cross-Frequency Decomposition: A Method for the Extraction of Neuronal Components Coupled at Different Frequencies. 
% Frontiers in neuroinformatics, 12.
%

% construct original signal
sig_1 = elec * filter;
sync = kuramoto ( sig_1, sig_2, Fr_1, Fr_2 );

% cut signals from target frequency band
cut_size = size(elec,1)/cut_num;
elec_split=cell(cut_num,1);
for i=1:cut_num
    elec_split{i} = elec(((i-1)*cut_size+1):(i*cut_size),:);
end

% permute signals
perm_sync=NaN(permute_num,1);
mismatch=NaN(permute_num,1);
for idx = 1:permute_num
%   disp(sprintf('Permutation %d of %d...', idx, permute_num));  
    perm_elec_1 = cell2mat(elec_split(randperm(cut_num)));
    
    elec_1 = detrend(perm_elec_1,'constant');
    m_1 = hilbert(elec_1);
    
    r = hilbert(sig_2);
    r = abs(r).*exp(1i.*(Fr_1*angle(r)));
    r=r(1:size(m_1,1),:);
    
    filter_1 = gcfd_optimize ( m_1, r, Fr_2, filter_seeds_num, lsqnonlin_param,0);
    perm_sig_1 = perm_elec_1 * filter_1;    
    perm_sync(idx) = kuramoto ( perm_sig_1, sig_2(1:size(m_1,1),:), Fr_1, Fr_2 ); 
    mismatch(idx) = vec_mismatch(filter_1,filter);
end
    
pd = fitdist(perm_sync,'Normal');
perm_error=1-cdf(pd,sync);

end

