function [filter_q] = gcfd_optimize ( m, r, FrP, filter_seeds_num, lsqnonlin_param, output)
% Core function for Generalized Cross-Frequency Decomposition
%
% INPUT:
%     m - signals at frequency FrQ, complexified, matrix TxM
%     r - signal at frequency FrP, FrP:FrQ = p:q, complexified, freq warped up q times, matrix Tx1
%         r is supposed to be one "reference" component, known a-priori from other means
%     FrP - FrP, integer number
%
%     f0_default, randomize_filter_seeds, filter_seeds_num - settings for seeds for iterative nonlinear optimization
%     lsqnonlin_param - parameters for "lsqnonlin" function (MATLAB built-in)
% OUTPUT:
%     filter_q - spatial filters (demixing coefficients) for frequency FrQ, matrix MxN
%
% References:
%
% Volk, D., Dubinin, I., Myasnikova, A., Gutkin, B., & Nikulin, V. V. (2018). 
% Generalized Cross-Frequency Decomposition: A Method for the Extraction of Neuronal Components Coupled at Different Frequencies. 
% Frontiers in neuroinformatics, 12.
%

SensorsNum = size(m,2);
PairsNum = size(r,2);

assert(PairsNum == 1);

%% nonlinear least-squares optimization
%% define the functions with phase shifts for 'lsqnonlin'
fun = @(f) abs((m*f).^FrP - r(:)); 

%% assemle the options structure for the algorithm
oo_names = fieldnames(lsqnonlin_param);
oo_values = struct2cell(lsqnonlin_param);
ooo(1:2:2*length(oo_names)-1) = oo_names;
ooo(2:2:2*length(oo_names)) = oo_values;
options = optimoptions('lsqnonlin', ooo{:});

%% set up the seeds and call 'lsqnonlin'

if output
    fprintf('randomize_filter_seeds with %d seeds:\n', filter_seeds_num);
end

filter_tmp = NaN(SensorsNum, filter_seeds_num);
resnorm = NaN(filter_seeds_num,1);
   
for seed_idx = 1:filter_seeds_num
    if output
        fprintf('Random seeds %d/%d:\n', seed_idx, filter_seeds_num);
    end
    f0 = randn(SensorsNum,1);
    f0 = f0/sum(abs(f0(:)));    
    
    [filter_tmp(:,seed_idx),resnorm(seed_idx),~,~,~] = lsqnonlin(fun,f0,[],[],options);
end
    
[M,I] = min(resnorm(:));

if output
    if filter_seeds_num > 1
        fprintf('Residual norms for %d seed filters:\n', filter_seeds_num);
        disp(resnorm(:));
        fprintf('Resnorms std/mean = %.4f\n\n', std(resnorm(:)) / mean(resnorm(:)));
    end
end
    
filter_q = filter_tmp(:,I(1));
 
end
