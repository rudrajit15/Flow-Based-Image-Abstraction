function [t_next1,imgret] = ETF2(image,win)
winby2 = floor(win/2);
inp = image;
sz = size(inp);
img = zeros(sz(1)+2*winby2,sz(2)+2*winby2);
sz1 = size(img);
i = 1:sz(1);
j = 1:sz(2);
img(winby2+i,winby2+j) = inp(i,j);                          %Zero Padding the input image.
img = double(img);                                          %Converting the image to double,
%img = imgaussfilt(img,1);                                 %Gaussian filtering the image.
imgret = img; 
[g(:,:,1) ,g(:,:,2)] = gradient(img);
[gmag, gdir] = imgradient(img);
normfact = 255*sqrt(2);
gmag = gmag/normfact;

t(:,:,1) = -g(:,:,2);
t(:,:,2) = +g(:,:,1);

t_next = zeros(sz1(1),sz1(2),2);
t(:,:,1) = t(:,:,1)/normfact;
t(:,:,2) = t(:,:,2)/normfact;
for run = 1:6
    for i= winby2+1:winby2+sz(1)
        for j = winby2+1:winby2+sz(2)
            sum_x = 0;
            sum_y = 0;
            norm = 0;
            for w_i = -winby2:winby2
                for w_j = -winby2:winby2
                    ws = 1;
                    wm = (tanh(gmag(i+w_i,j+w_j) - gmag(i,j)) + 1)/2;
                    wd = t(i,j,1)*t(i+w_i,j+w_j,1) + t(i,j,2)*t(i+w_i,j+w_j,2);
                    sum_x = sum_x + t(i+w_i,j+w_j,1)*ws*wm*wd;
                    sum_y = sum_y + t(i+w_i,j+w_j,2)*ws*wm*wd;
                    norm = norm + ws*wm*wd;
                end
            end
            t_next(i,j,1) = sum_x/normfact;
            t_next(i,j,2) = sum_y/normfact;
        end
    end
    t = t_next;
end
t_next1 = t;
[tmag1 ,tdir1] = imgradient(t_next1(:,:,1),t_next1(:,:,2));
imshow(tdir1);
