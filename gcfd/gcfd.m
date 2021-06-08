function [W_P, W_Q, A_P, A_Q, gcfd_P, gcfd_Q, p_values] = gcfd(X, FrBase, FrP, FrQ, PairsNum, sampling_freq)

% GCFD - Generalized Cross-Frequency Decomposition
%
% [W_P, W_Q, A_P, A_Q, gcfd_P, gcfd_Q, p_values] = gcfd(X, FrBase, FrP, FrQ, PairsNum, sampling_freq)
%
% This function is for extraction of cross-frequency phase-to-phase synchronized components. 
% It reconstructs the time courses of synchronized neuronal components, their spatial filters and patterns. 
% The method extends the previous state of the art, Cross-Frequency Decomposition* (CFD)
% to the whole range of frequencies: the method works for any f1 and f2 whenever f1:f2 is a rational number. 
% 
% INPUT: 
%     X        -  a matrix of size TxM with T samples and M channels of data
%     FrBase   -  base frequency, Hz - to be multiplied by FrP, FrQ  
%     FrP      -  frequency ratio f1 - (reference band for warp) 
%     FrQ      -  frequency ratio f2 - (target band) 
%     PairsNum -  number of synchronized source pairs to be found 
%     sampling_freq - the sampling frequency (in Hz) of X (the data)
%                         
% OUTPUT:
%     W_P, W_Q -        the de-mixing matrix. Each column is a filter and the
%                       timecourse of found components is extracted with X * W
%     A_P, A_Q -        the patterns (mixing matrix) with the i'th column
%                       corrresponding to the i'th component
%     gcfd_P, gcfd_Q -  the bandpass filtered data projected onto found components
%                       This is two matrices of size TxM each with T samples and M number of synchronized pairs (PairsNum). 
%                       For example - first columns of gcfd_P and gcfd_Q correspond to the first pair of phase-to-phase synchronized components.
%     p_value  -        p-value for correspoding pair of components based
%                       on non-parametric permutation tests
% 
% References:
%
% Volk, D., Dubinin, I., Myasnikova, A., Gutkin, B., & Nikulin, V. V. (2018). 
% Generalized Cross-Frequency Decomposition: A Method for the Extraction of Neuronal Components Coupled at Different Frequencies. 
% Frontiers in neuroinformatics, 12.
%
% *
% Nikulin V. V., Nolte G., Curio G.
% Cross-frequency decomposition: A novel technique for studying interactions between neuronal oscillations with different frequencies
% Clinical Neurophysiology, 2012
%
% MIT License
%
% Copyright (c) [2018] [Denis Volk, Igor Dubinin, Alexandra Myasnikova, Boris Gutkin, Vadim Nikulin]
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


%% Default parameters

epochs = [1 size(X,1)];     % default epochs
fr_cutoff_low = 0.5;        % for preprocessing
fr_cutoff_high = 50.0;


ssd_primary_on = true;          % ssd decompostion in f1 band
ssd_secondary_on = true;        % ssd decompostion in f2 band for computation speed 
ssd_secondary_comp_num = 15;    % set to Inf to use all the ssd_q components

% parameters for non-linear optimization
lsqnonlin_param.Algorithm = 'levenberg-marquardt';
lsqnonlin_param.Display = 'off';
lsqnonlin_param.TolFun = 1e-8;
lsqnonlin_param.MaxFunEvals = 1000000;
lsqnonlin_param.MaxIter = 100;
filter_seeds_num=2;
output=true;

% parameters for permutation tests. For help read sync_elec_permut
stat_tests = false;   % to perform tests set true
permute_num = 1000;
cut_num = 100; 

%% Preprocess the raw data
[butter_preproc_b, butter_preproc_a] = butter(2,[fr_cutoff_low fr_cutoff_high]/(sampling_freq/2));
X = filtfilt (butter_preproc_b, butter_preproc_a, X);

%% Create filters for bands
%[butter_p_b,butter_p_a]=butter(2,[FrP*FrBase-1 FrP*FrBase+1]/(sampling_freq/2));    % FrP band-pass filter
%[butter_q_b,butter_q_a]=butter(2,[FrQ*FrBase-1 FrQ*FrBase+1]/(sampling_freq/2));    % FrQ band-pass filter

%% run SSD to recover FrP sources
if ssd_primary_on
    [W_P, pattern_p_ssd, ~, ~, X_ssd_p] = ssd(X, [FrP*FrBase-1 FrP*FrBase+1; FrP*FrBase-3 FrP*FrBase+3; FrP*FrBase-2 FrP*FrBase+2], sampling_freq, 2, epochs);
end

%% run SSD to recover FrQ sources
if ssd_secondary_on
    [W_Q, pattern_q_ssd, ~, ~, X_ssd_q] = ssd(X, [FrQ*FrBase-1 FrQ*FrBase+1; FrQ*FrBase-3 FrQ*FrBase+3; FrQ*FrBase-2 FrQ*FrBase+2], sampling_freq, 2, epochs);
end
        
%% genCFD DECOMPOSITION PART:
%% prepare the signals for P - Q search
gcfd_elec_q_num = NaN;
if ssd_secondary_on
    gcfd_elec_q_num = min([size(X_ssd_q,2) ssd_secondary_comp_num]);
    gcfd_elec_q = X_ssd_q(:,1:gcfd_elec_q_num);
else
    gcfd_elec_q_num = size(X,2);
    gcfd_elec_q = filtfilt(butter_q_b,butter_q_a,X);
end

filter_q = NaN(gcfd_elec_q_num,PairsNum);
gcfd_elec_q = detrend(gcfd_elec_q,'constant');
X_ssd_p = detrend(X_ssd_p,'constant');

%% complexify the signals (compute 'analytic signals')
m_Q = hilbert(gcfd_elec_q);
r_P = hilbert(X_ssd_p);

%% frequency warp the SSD component
r_P = abs(r_P).*exp(1i.*(FrQ*angle(r_P)));

%% Prepare the seeds for the nonlinear optimization
f0_default = ones(gcfd_elec_q_num,PairsNum)/gcfd_elec_q_num;

%% Non-linear optimization for each of FrP sources
for idx = 1:PairsNum
           fprintf('[%d/%d] Nonlinear least-squares optimization for SSD component %d of %d ...', idx, PairsNum, idx, PairsNum);           
           filter_q(:,idx) = gcfd_optimize ( m_Q, r_P(:,idx), FrP, filter_seeds_num, lsqnonlin_param,output);
end 

%% Filter to pattern for all FrQ sources
clear pattern_q_gcfd
pattern_q_gcfd = NaN(gcfd_elec_q_num,PairsNum); 
for idx = 1:PairsNum
    pattern_q_gcfd(1:gcfd_elec_q_num,idx) = filter2pattern (filter_q(:,idx), gcfd_elec_q);
end

%% Convert the pattern in SSD coordinates to the pattern in original coordinates
if ssd_secondary_on    
    pattern_q_gcfd (gcfd_elec_q_num+1:size(X_ssd_q,2),:) = 0;
    pattern_q = pattern_q_ssd * pattern_q_gcfd;
else
    pattern_q = pattern_q_gcfd;
end

%% Output patterns and componets
W_P = W_P(:,1:PairsNum);
W_Q = filter_q;

A_P = pattern_p_ssd(:,1:PairsNum);
A_Q = pattern_q;

gcfd_P = X_ssd_p(:,1:PairsNum);
gcfd_Q = gcfd_elec_q*W_Q; 

%% Permutation tests
p_values = NaN(1,PairsNum);
if stat_tests    
    for idx=1:PairsNum
        fprintf('[%d/%d] Permutation test %d of %d ...', idx, PairsNum, idx, PairsNum); 
        [p_values(1,idx), ~, ~] = sync_elec_permut( filter_q(:,idx), gcfd_elec_q, X_ssd_p(:,idx) , FrQ, FrP, cut_num, permute_num, lsqnonlin_param, filter_seeds_num);
    end
end

end