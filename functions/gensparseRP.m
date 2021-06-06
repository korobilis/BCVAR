function [PHI_LARGE] = gensparseRP(p,m_l,KM)
%==========================================================================
% Function generates random projection matrix
%==========================================================================
% INPUTS:
%    RP_type:  Random projection type, 1: G+Dunson case; 
%                                      2: Achlioptas sparse RP matrix
%                                      3: RP matrix from N(0,1)
%    m_l: Dimension of compressed model
%    KM:  Dimension of original model
%==========================================================================
PHI_LARGE = zeros(m_l*p,KM);
for irep = 1:p
    % Very sparse RP matrix, based on Achlioptas
    s = sqrt(irep*(KM/p))/2;
    prob(1,1) = 1/(2*s); prob(2,1) = 1 - (1/s); prob(3,1) = 1/(2*s);
    cum_prob = cumsum(prob,1);
    PHI_val(1,1) = sqrt(s); PHI_val(2,1) = 0; PHI_val(3,1) = -sqrt(s);
                    
    PHI = rand(m_l,KM);
    GG = (PHI<=cum_prob(1,1));
    GG2 = (PHI<=cum_prob(2,1)) & (PHI>cum_prob(1,1));
    GG3 = (PHI>cum_prob(2,1));
    PHI(GG) = PHI_val(1); PHI(GG2) = PHI_val(2); PHI(GG3) = PHI_val(3);
    PHI_LARGE((irep-1)*m_l+1:m_l*irep,:) = PHI;
end
