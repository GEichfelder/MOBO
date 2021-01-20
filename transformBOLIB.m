function val = transformBOLIB(x, y, keyf, keyxy)
%TRANSFORMBOLIB transform BOLIB Model to Model used in the Paper of
%G. Eichfelder.

global ModelBOLIB n_1 n_2 p p_tilde

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; val = feval(ModelBOLIB, y, x, 'F');
    case 'G'; val = feval(ModelBOLIB, y, x, 'g');    
    case 'f'; val = feval(ModelBOLIB, y, x, 'f');
    case 'G_tilde'; val = feval(ModelBOLIB, y, x, 'G');
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; val = feval(ModelBOLIB, y, x, 'F', 'y').'; 
        case 'y' ; val = feval(ModelBOLIB, y, x, 'F', 'x').';
        case 'xx'; val = feval(ModelBOLIB, y, x, 'F', 'yy');
        case 'xy'; val = feval(ModelBOLIB, y, x, 'F', 'xy').'; % easy because F is onedimensional
        case 'yy'; val = feval(ModelBOLIB, y, x, 'F', 'xx');    
        end  
    case 'G'
        switch keyxy
        case 'x' ; val = feval(ModelBOLIB, y, x, 'g', 'y'); 
        case 'y' ; val = feval(ModelBOLIB, y, x, 'g', 'x');   
        case 'xx'; val = feval(ModelBOLIB, y, x, 'g', 'yy');
        case 'xy'
            val_BOLIB = feval(ModelBOLIB, y, x, 'g', 'xy');
            val = zeros(p*n_2, n_1);
            for i = 1:p
                val((i-1)*n_2+1:i*n_2, :) = val_BOLIB((i-1)*n_1+1:i*n_1, :).';
            end
        case 'yy'; val = feval(ModelBOLIB, y, x, 'g', 'xx');    
        end          
	case 'f'   
        switch keyxy
        case 'x' ; val = feval(ModelBOLIB, y, x, 'f', 'y').';          
        case 'y' ; val = feval(ModelBOLIB, y, x, 'f', 'x').';   
        case 'xx'; val = feval(ModelBOLIB, y, x, 'f', 'yy');
        case 'xy'; val = feval(ModelBOLIB, y, x, 'f', 'xy').';  % easy because f is onedimensional
        case 'yy'; val = feval(ModelBOLIB, y, x, 'f', 'xx');    
        end           
	case 'G_tilde'   
        switch keyxy
        case 'x' ; val = feval(ModelBOLIB, y, x, 'G', 'y');             
        case 'y' ; val = feval(ModelBOLIB, y, x, 'G', 'x');
        case 'xx'; val = feval(ModelBOLIB, y, x, 'G', 'yy');
        case 'xy'
            val_BOLIB = feval(ModelBOLIB, y, x, 'G', 'xy');
            val = zeros(p_tilde*n_2, n_1);
            for i = 1:p_tilde
                val((i-1)*n_2+1:i*n_2, :) = val_BOLIB((i-1)*n_1+1:i*n_1, :).';
            end
        case 'yy'; val = feval(ModelBOLIB, y, x, 'G', 'xx');     
        end        
   end   
end