function morphed_im = morph_tps(im_source, a1_x, ax_x, ay_x, w_x,...
                        a1_y, ax_y, ay_y, w_y, ctr_pts, sz)
    'morph_tps'
    morphed_im = zeros(sz(1), sz(2), 3);
    for x = 1:sz(2)
        for y = 1:sz(1)
            U = sqrt(((ctr_pts(:,1) - x) .* (ctr_pts(:,1) - x)) + ...
                     ((ctr_pts(:,2) - y) .* (ctr_pts(:,2) - y)))';
            U = spfun(@(r) (r.^2) .* log(r.^2), U);
            M = [1, x, y, U];
            f_x = min(sz(2), max(1, round(M*[a1_x;ax_x;ay_x;w_x])));
            f_y = min(sz(1), max(1, round(M*[a1_y;ax_y;ay_y;w_y])));
            morphed_im(y,x,:) = im_source(f_y, f_x,:);
        end
    end
end
