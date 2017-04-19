% read_imagecube
%
% 8-class 6D synthetic image cube
function [X, Y] = read_imagecube()

    %% Import the data
    nrow = 128; % 128 rows in each image band
    ncol = 128; % 128 columns in each image band
    dim = 6; % 6 image bands
    datalen = nrow*ncol; % number of entries in the file

    %read data file
    fid  = fopen('~/pjh5/neural_nets/data/8class-caseI/8class-caseI-n_m0v5.viff','r'); 

    % Skip standard 1024-byte header
    [~,~] = fread(fid,1024);

    %       Reading into a datalen x 8 array for easy handling
    %       (making each image band into one continuous string of
    %       pixel values)
    %       This format is just like the 'iris' 
    x = zeros(datalen,dim); 
    for i = 1:dim
       [x(:,i),~] = fread(fid,datalen,'float32'); % 32 bit float
    end

    % extract by band for use
    full_data = zeros(nrow,ncol,6);
    for i = 1:nrow
        for k = 1:ncol
            full_data(i,k,1:6) = x((i-1) * ncol + k, 1:6);
        end
    end

    % Visualize Data
    %{
    for i = 1:dim
       figure;
       imagesc(full_data(:,:,i));
       caxis([0, 255]);
       colormap(gray);
       colorbar;
       title(['Band ', num2str(i)], 'fontsize', 16);
    end
    %}

    % generate the labels
    label                  = zeros(nrow,ncol);
    label(  1:64 ,  1:64 ) = 1;    %red
    label( 33:64 , 97:128) = 2;    %orange
    label(  1:64 , 65:96 ) = 3;    %green
    label(  1:32 , 97:128) = 4;    %yellow
    label( 65:128, 65:128) = 5;    %white
    label( 65:96 ,  1:32 ) = 6;    %blue
    label( 97:128,  1:32 ) = 7;    %purple
    label( 65:128, 33:64 ) = 8;    %gray

    % generate the color labels
                      %r    or  g   y   w   b    p gray
    label_color_map_r = [255;254;  0;252;255; 47;164;180];
    label_color_map_g = [  0;146;255;228;255;146;  0;180];
    label_color_map_b = [  0; 43;  0; 98;255;255;175;180];
    label_color_map   = [label_color_map_r label_color_map_g label_color_map_b]/255;
    % example of how to plot the labels for the whole image:
    % imagesc(label); colormap(label_color_map);
    
    
    % example of how to plot the labels of the region included:
    % to visualize the tr_tst_index you can use.
    % imagesc(tr_tst_index); colorbar;

    % convert from the data from matrix to vector form 
    Y = zeros(datalen, 1);
    X = zeros(datalen, dim);
    i = 1;
    for     row = 1:nrow
        for col = 1:ncol
            X(i, :) = full_data(row, col, :);
            Y(i) = label(row, col);
            i = i + 1;
        end
    end


end

