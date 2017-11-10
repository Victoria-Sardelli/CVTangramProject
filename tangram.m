function tangram
    % Analyzes tangram images to detect pieces based on color and shape
    %
    % Authors: Victoria Sardelli and Roxanne Meskell
    %
    % Note: Uses Dr. Kinsman's code to add salt and pepper noise to images

    % call function to analyze tangram images
    isolate_pieces();
    

    function im_dirty = flip_fraction_of_bits( im_clean, fractional_amt )
        %  Treating the image as a vector of values in the range [0,1], 
        %  flip "n_percent" of the bits.  Effectively this adds salt-and-pepper
        %
        %  Input image is a double, in the range [0,1].
        %
        % Author: Dr. Thomas Kinsman
        
        if ( fractional_amt <= 0 ) || ( fractional_amt >= 1 )
            error('fractional_amt of image is out of vaiable range for use.');
        end
    
        n_pix_per_plane = numel( im_clean(:,:,1) );
    
        n_to_flip       = round( fractional_amt * n_pix_per_plane );
        rnd_vec         = randperm(  n_pix_per_plane );
    
        im_dirty        = im_clean;
        v_to_flip       = rnd_vec( 1:n_to_flip );
    
        %
        %  To add grayscale noise, we need to modify the (red,green,blue) 
        % channels all the same way.
        %
        im_dirty(v_to_flip) = 1 - im_dirty(v_to_flip); % Modify the r channel.
        v_to_flip = v_to_flip + n_pix_per_plane; % Advance to the g channel.
    
        im_dirty(v_to_flip) = 1 - im_dirty(v_to_flip); % Modify the g channel
        v_to_flip = v_to_flip + n_pix_per_plane; % Advance to the b channel.
    
        im_dirty(v_to_flip) = 1 - im_dirty(v_to_flip); % Modify the b channel.
    end

    
    function [cluster_ctr, segmented_images] = kmeans_color_seg(im,k,numColors)
        % apply kmeans clustering algo to image to segment based on colors
       
        % convert image to L*a*b* color space
       im_lab = rgb2lab( im );
       imshow( im_lab );
       title( 'L*a*b* color space Image' );
       pause( 3 );
        
       % we want the color info from the a* b* space
       ab = double(im_lab(:,:,2:3));
       num_rows = size(ab,1);
       num_cols = size(ab,2);
       ab = reshape(ab, num_rows*num_cols, 2);
       
       % repeat clustering 3 times to avoid local minima
       [cluster_idx, cluster_ctr] = kmeans(ab,numColors,'distance',...
           'sqEuclidean','Replicates',3);
       
       % label every pixel with its cluster idx
       pixel_labels = reshape(cluster_idx, num_rows, num_cols);
       imshow(pixel_labels,[]), title('Image labeled by cluster index');
       pause(3);
       
       % using pixel_labels, separate objects by color
       segmented_images = cell(1,numColors);
       rgb_label = repmat(pixel_labels,[1 1 3]);
       
       for k = 1:numColors
           color = im;
           color(rgb_label ~= k) = 0;
           segmented_images{k} = color; 
       end
       
    end


    function isolate_pieces
        % isolates tangram pieces in image
       
        % read in original image
        im = imread('tangram2.png');
        imshow( im );
        title('Original Image');
        pause( 3 );
        
        % add salt-and-pepper noise
        im = flip_fraction_of_bits( im, 0.03 );
        imshow(im);
        title('Noisy Image');
        pause(3);
        
        % clean noisy image with median filter
        r = im(:,:,1);
        g = im(:,:,2);
        b = im(:,:,3);
        im(:,:,1) = medfilt2(r);
        im(:,:,2) = medfilt2(g);
        im(:,:,3) = medfilt2(b);
        imshow(im);
        title('Cleaned Image');
        pause(3);
        
       % cluster with kmeans
       [cluster_ctr, segmented_images] = kmeans_color_seg(im, 8, 8);
       
       % now that we have to figure out which cluster_idx has sky and 
       % not show that one
       mean_cluster_val = mean(cluster_ctr,2);
       [arr, idx] = sort(mean_cluster_val);
       im_tangram = zeros; % will hold final image
       
       % add code here to loop through segmented images
       % erode with strel cube size 2 
       % find square -> if found, keep track of which im it is
       % find triangles in remaining images
       % two images should remain -> parallelogram and background
       
       %for curr_idx = 1:8
            %if curr_idx ~= sky_num
                %im_tangram = im_tangram + segmented_images{curr_idx};
            %end
       %end
       
       imshow(im_tangram);
       title('Tangram Image');
       pause(3);
       
    end

end