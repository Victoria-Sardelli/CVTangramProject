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