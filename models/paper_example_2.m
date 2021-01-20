function val=paper_example_2(x,y,keyf,keyxy)
% [dim_x dim_y dim_F dim_f dim_G dim_G_tilde] = [14 1 2 2 1 2]   
global V

A = zeros(6, 8);
A(1, 1) = 1;
A(2, 2) = 1;
A(3, 3) = 1;
sin_y = sin(y);
cos_y = cos(y);
A(6, 4) = -sin_y;
A(4, 5) = cos_y;
A(4, 6) = sin_y;
A(5, 7) = sin_y;
A(6, 7) = cos_y;
A(5, 8) = cos_y;
A(6, 8) = -sin_y;
b = [0; cos_y; sin_y; 1; 0; 0];

     
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'
        x_old = [0.1247; 0.1335; -0.0762; -0.1690; 0.2118; -0.0534; -0.1473; ...
         0.3170; -0.0185; -0.1800; 0.1700; -0.0718; 0.0058; 0.0985];
        val = [norm(x)^2; ...
               norm(x - x_old)^2];
    case 'G'; val = norm(A*V*x-b)^2 - 0.3;    
    case 'f'; val = [norm(A*V*x-b)^2; norm(x)^2];
    case 'G_tilde'; val = [y - pi; ...
                           0 - y];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; val = []; 
        case 'y' ; val = [];
        case 'xx'; val = [];
        case 'xy'; val = [];
        case 'yy'; val = [];    
        end  
    case 'G'
        switch keyxy
        case 'x' ; val = [2*(A*V*x-b).'*A*V];          
        case 'y' 
            dA = zeros(6, 8);
            dA(6, 4) = -cos_y;
            dA(4, 5) = -sin_y;
            dA(4, 6) = cos_y;
            dA(5, 7) = cos_y;
            dA(6, 7) = -sin_y;
            dA(5, 8) = -sin_y;
            dA(6, 8) = -cos_y;
            db = [0; -sin_y; cos_y; 0; 0; 0];
            val = [2*(A*V*x-b).'*(dA*V*x-db)];   
        case 'xx'; val = [2*(A*V).'*(A*V)];
        case 'xy'
            dA = zeros(6, 8);
            dA(6, 4) = -cos_y;
            dA(4, 5) = -sin_y;
            dA(4, 6) = cos_y;
            dA(5, 7) = cos_y;
            dA(6, 7) = -sin_y;
            dA(5, 8) = -sin_y;
            dA(6, 8) = -cos_y;
            db = [0; -sin_y; cos_y; 0; 0; 0];
            val = [2*((dA*V*x-db).'*A*V+(A*V*x-b).'*dA*V)];
        case 'yy'
            dA = zeros(6, 8);
            dA(6, 4) = -cos_y;
            dA(4, 5) = -sin_y;
            dA(4, 6) = cos_y;
            dA(5, 7) = cos_y;
            dA(6, 7) = -sin_y;
            dA(5, 8) = -sin_y;
            dA(6, 8) = -cos_y;
            db = [0; -sin_y; cos_y; 0; 0; 0];
            d2A = zeros(6, 8);
            d2A(6, 4) = sin_y;
            d2A(4, 5) = -cos_y;
            d2A(4, 6) = -sin_y;
            d2A(5, 7) = -sin_y;
            d2A(6, 7) = -cos_y;
            d2A(5, 8) = -cos_y;
            d2A(6, 8) = sin_y;
            d2b = [0; -cos_y; -sin_y; 0; 0; 0];
            val = [2*((dA*V*x-db).'*(dA*V*x-db)+ (A*V*x-b).'*(d2A*V*x-d2b))];        
        end          
	case 'f'   
        switch keyxy
        case 'x' ; val = [2*(A*V*x-b).'*A*V;
                          2*x.'];          
        case 'y' 
            dA = zeros(6, 8);
            dA(6, 4) = -cos_y;
            dA(4, 5) = -sin_y;
            dA(4, 6) = cos_y;
            dA(5, 7) = cos_y;
            dA(6, 7) = -sin_y;
            dA(5, 8) = -sin_y;
            dA(6, 8) = -cos_y;
            db = [0; -sin_y; cos_y; 0; 0; 0];
            val = [2*(A*V*x-b).'*(dA*V*x-db);
                   0];   
        case 'xx'; val = [2*(A*V).'*(A*V); ...
                          2*eye(14)];
        case 'xy'
            dA = zeros(6, 8);
            dA(6, 4) = -cos_y;
            dA(4, 5) = -sin_y;
            dA(4, 6) = cos_y;
            dA(5, 7) = cos_y;
            dA(6, 7) = -sin_y;
            dA(5, 8) = -sin_y;
            dA(6, 8) = -cos_y;
            db = [0; -sin_y; cos_y; 0; 0; 0];
            val = [2*((dA*V*x-db).'*A*V+(A*V*x-b).'*dA*V);...
                   zeros(1, 14)];
        case 'yy'
            dA = zeros(6, 8);
            dA(6, 4) = -cos_y;
            dA(4, 5) = -sin_y;
            dA(4, 6) = cos_y;
            dA(5, 7) = cos_y;
            dA(6, 7) = -sin_y;
            dA(5, 8) = -sin_y;
            dA(6, 8) = -cos_y;
            db = [0; -sin_y; cos_y; 0; 0; 0];
            d2A = zeros(6, 8);
            d2A(6, 4) = sin_y;
            d2A(4, 5) = -cos_y;
            d2A(4, 6) = -sin_y;
            d2A(5, 7) = -sin_y;
            d2A(6, 7) = -cos_y;
            d2A(5, 8) = -cos_y;
            d2A(6, 8) = sin_y;
            d2b = [0; -cos_y; -sin_y; 0; 0; 0];
            val = [2*((dA*V*x-db).'*(dA*V*x-db)+ (A*V*x-b).'*(d2A*V*x-d2b)); ...
                   0];    
        end           
	case 'G_tilde'   
        switch keyxy
        case 'x' ; val = [];             
        case 'y' ; val = [];
        case 'xx'; val = [];
        case 'xy'; val = [];
        case 'yy'; val = [];     
        end        
   end   
end

end




