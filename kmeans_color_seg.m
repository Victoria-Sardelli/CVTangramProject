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