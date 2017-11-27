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