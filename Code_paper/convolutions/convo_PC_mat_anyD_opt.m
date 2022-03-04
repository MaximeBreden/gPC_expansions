function Mb = convo_PC_mat_anyD_opt(b,alpha_mat,extra_dim,output_format)

% Constructs the generalized convolution matrix for b associated to the 
% linearization coefficients in alpha_mat. 
%
% b can be a vector, or a higher dimensional tensor, and alpha_mat can be
% a matrix or a cell array of matrices. If extra_dim > 0, then we do not
% consider generalized convolutions along the extra_dim first dimensions in
% b. For instance if b is a tensor containing Fourier x gPC coefficients,
% by putting extra_dim = 1 we produce generalized convolution matrices for
% each Fourier coefficient in b. In such a case, output_format selects in
% which way these matrices are concatenated.

if nargin < 4
    output_format = 'column';
    if nargin < 3
        extra_dim = 0;
    end
end

if not(iscell(alpha_mat))
    alpha_mat = {alpha_mat};
    dim_alpha = 1;
else
    dim_alpha = length(alpha_mat);
end

sz = [size(b), ones(1,1+dim_alpha+extra_dim-length(size(b)))];

Mb = transpose(tens2vect(b));

%% Specific (faster) code for the first dimension
l = dim_alpha+extra_dim;
sz_next = prod(sz(1:l-1));
Mb = alpha_mat{dim_alpha} * reshape(Mb, sz(l), sz_next);

if dim_alpha == 1
    % no extra dimensio needed, we reshape the output and exit here
    switch output_format
        case 'column'
            Mb = reshape(Mb,[sz(l),sz(l)*sz_next]);
        case 'row'
            Mb = reshape( permute( reshape(Mb,[sz(l),sz(l),sz_next]), [1,3,2]), [sz(l)*sz_next,sz(l)] );
        otherwise
            error('This format is not valid')
    end  
else
%% Next dimensions if needed
    if dim_alpha == 2 && extra_dim == 0
        % shortcut for pure 2D
        Mb = reshape( permute( reshape(Mb,[sz(l),sz(l),sz_next]), [1,3,2]), [sz(l)*sz_next,sz(l)] );
    else
        sz_block = prod(sz(l+1:end));
        sz_rest = prod(sz(1:l));
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
    for l = (dim_alpha-1:-1:1)+extra_dim
        sz_block = prod(sz(l+1:end));
        sz_rest = prod(sz(1:l));
        sz_next = prod(sz(1:l-1));
        if dim_alpha ~= 2 || extra_dim ~= 0
            % avoided for pure 2D
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
        end
        Mb = my_kron( alpha_mat{l-extra_dim}, eye(sz_block) ) * Mb;
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