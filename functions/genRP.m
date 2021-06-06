function [PHI] = genRP(RP_type,m_l,KM)
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

if RP_type == 1 
    % G+Dunson case, we need to put this inside the for loop below    
    psi = unifrnd(0.1,0.8,1,1);       
elseif RP_type == 2
    % Sparse RP matrix, based on Achlioptas
    prob(1,1) = 1/6; prob(2,1) = 2/3; prob(3,1) = 1/6;
    cum_prob = cumsum(prob,1);
    PHI_val(1,1) = sqrt(3); PHI_val(2,1) = 0; PHI_val(3,1) = -sqrt(3);
end

% -------| Generate PHI      
if RP_type == 1   % G+Dunson case
    prob(1,1) = psi^2; prob(2,1) = 2*(1-psi)*psi; prob(3,1) = (1-psi)^2;
    PHI_val(1,1) = -sqrt(1/psi); PHI_val(2,1) = 0; PHI_val(3,1) = sqrt(1/psi);
    cum_prob = cumsum(prob,1);
    PHI = rand(m_l,KM);
    GG = (PHI<=cum_prob(1,1));
    GG2 = (PHI<=cum_prob(2,1)) & (PHI>cum_prob(1,1));
    GG3 = (PHI>cum_prob(2,1));
    PHI(GG) = PHI_val(1); PHI(GG2) = PHI_val(2); PHI(GG3) = PHI_val(3);
elseif RP_type == 2    % Achlioptas (2003) Case I                     
    PHI = rand(m_l,KM);
    GG = (PHI<=cum_prob(1,1));
    GG2 = (PHI<=cum_prob(2,1)) & (PHI>cum_prob(1,1));
    GG3 = (PHI>cum_prob(2,1));
    PHI(GG) = PHI_val(1); PHI(GG2) = PHI_val(2); PHI(GG3) = PHI_val(3);
elseif RP_type == 3        
    PHI = rand(m_l,KM);
    GG = (PHI<=0.5);
    GG2 = (PHI>0.5);
    PHI(GG) = -1; PHI(GG2) = 1;   
elseif RP_type == 4
    PHI = unifrnd(-1,1,m_l,KM);
elseif RP_type == 5
    PHI = randi([-1,1],m_l,KM); 
end
if m_l > 1 && RP_type == 1 
    [PHI] = grams(PHI);    % Grahm-Shmidt orthonormalization
end
