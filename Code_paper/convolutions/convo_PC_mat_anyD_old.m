function Mb = convo_PC_mat_anyD_old(b,alpha_mat,output_format)

% Constructs the generalized convolution matrix associated to the 
% linearization coefficients in alpha_mat. 
%
% b can be a vector, or a higher dimensional tensor, and alpha_mat can be
% a matrix or a cell array of matrices. The length L_convo of the cell 
% array must be at most equal to the number L_b of dimensions of b. If
% L_convo = L_b, we consider for each dimension l the generalized
% convolution product whose coefficients are stored in alpha_mat{l}. If
% L_convo < L_b, we consider generalized convolution products only for the 
% dimensions starting from L_b-L_convo+1 in b. In that case, the shape of
% the output can be specified with a third optional argument.

if not(iscell(alpha_mat))
    alpha_mat = {alpha_mat};
    L_convo = 1;
else
    L_convo = length(alpha_mat);
end
    
if isvector(b) 
    L_b = 1;
    sz = [length(b), 1];
else
    L_b = ndims(b);
    sz = [size(b), 1];
end

shift = L_b-L_convo; 

Mb = transpose(tens2vect(b));

for l = L_b:-1:(1+shift)
    sz_block = prod(sz(l+1:end));
    sz_rest = prod(sz(1:l));
    sz_next = prod(sz(1:l-1));
    if exist('intval','file') && isintval(b(1))
        %%% the mat2cell/cell2mat conversion is slow with intvals, so we do
        %%% it "by hand" instead
        Mb_mid = mat2cell( mid(Mb), sz_block, sz_block*ones(1,sz_rest) );
        Mb_mid = cell2mat( reshape( Mb_mid, sz(l), sz_next ) );       
        Mb_rad = mat2cell( rad(Mb), sz_block, sz_block*ones(1,sz_rest) );
        Mb_rad = cell2mat( reshape( Mb_rad, sz(l), sz_next ) );
        Mb = midrad(Mb_mid,Mb_rad);
    else
        Mb = mat2cell( Mb, sz_block, sz_block*ones(1,sz_rest) );
        Mb = cell2mat( reshape( Mb, sz(l), sz_next ) );
    end
    Mb = my_kron( alpha_mat{l-shift}, speye(sz_block) ) * Mb;
    if exist('intval','file') && isintval(b(1))        
        %%% the mat2cell/cell2mat conversion is slow with intvals, so we do
        %%% it "by hand" instead
        Mb_mid = mat2cell( mid(Mb), sz_block*ones(1,sz(l)^2), sz_block*ones(1,sz_next) );
        Mb_mid = cell2mat( reshape( Mb_mid, sz(l), sz_rest ) );
        Mb_rad = mat2cell( rad(Mb), sz_block*ones(1,sz(l)^2), sz_block*ones(1,sz_next) );
        Mb_rad = cell2mat( reshape( Mb_rad, sz(l), sz_rest ) );
        Mb = midrad(Mb_mid,Mb_rad);
    else
        Mb = mat2cell( Mb, sz_block*ones(1,sz(l)^2), sz_block*ones(1,sz_next) );
        Mb = cell2mat( reshape( Mb, sz(l), sz_rest ) );
    end
end

if nargin == 3
    switch output_format
        case 'column'
            %Do nothing, the output is already arranged column-wise 
        case 'row'
            l = l-1;
            sz_block = prod(sz(l+1:end));
            sz_rest = prod(sz(1:l));
            if exist('intval','file') && isintval(b(1))  
                %%% the mat2cell/cell2mat conversion is slow with intvals, so we do
                %%% it "by hand" instead
                Mb_mid = mat2cell( mid(Mb), sz_block, sz_block*ones(1,sz_rest) );
                Mb_mid = cell2mat( transpose(Mb_mid) );   
                Mb_rad = mat2cell( rad(Mb), sz_block, sz_block*ones(1,sz_rest) );
                Mb_rad = cell2mat( transpose(Mb_rad) );  
                Mb = midrad(Mb_mid,Mb_rad); 
            else
                Mb = mat2cell( Mb, sz_block, sz_block*ones(1,sz_rest) );
                Mb = cell2mat( transpose(Mb) );
            end
        otherwise
            error('This output format is not valid')
    end
end

