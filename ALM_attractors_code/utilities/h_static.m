function [ h_out ] = h_static(y,tauD,tauf,U)
%Static transfer function including short-term synaptic plasticity
%   from : Tsodyks, Misha; Pawelzik, Klaus and Markram, Henry (1998). 
%   Neural Networks with Dynamic Synapses. Neural Computation. 10(4): 821-835. 
%   doi:10.1162/089976698300017502.doi:10.1162/089976698300017502
    
    h_out = (U + y.*U*tauf)./(1 + y.^2*U*tauD*tauf + y.*U*(tauD + tauf));
        
end

