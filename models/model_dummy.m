function val=model_dummy(x,y,keyf,keyxy)
 
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'f'; val = [];         % m_1 column vector
    case 'F'; val = [];         % m_2 column vector
    case 'G'; val = [];         % p column vector
    case 'G_tilde'; val = [];   % p_tilde column vector
    end    
else
    switch keyf
	case 'f'   
        switch keyxy
        case 'x' ; val = [];    % m_2 x n_1 matrix      
        case 'y' ; val = [];    % m_2 x n_2 matrix
        case 'xx'; val = [];    % m_2*n_1 x n_1 matrix
        case 'xy'; val = [];    % m_2*n_2 x n_1 matrix
        case 'yy'; val = [];    % m_2*n_2 x n_2 matrix
        end 
    case 'G'
        switch keyxy
        case 'x' ; val = [];    % p x n_1 matrix
        case 'y' ; val = [];    % p x n_2 matrix
        case 'xx'; val = [];    % p*n_1 x n_1 matrix
        case 'xy'; val = [];    % p*n_2 x n_1 matrix
        case 'yy'; val = [];    % p*n_2 x n_2 matrix
        end
    end
    
end