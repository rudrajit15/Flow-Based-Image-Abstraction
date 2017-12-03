function [img2, Tx, Ty] = myETF(img, mu, sigma)
img = imgaussfilt(img,sigma);
s = size(img);
s2_1 = s(1) + 2*mu;
s2_2 = s(2) + 2*mu;
s2 = [s2_1 s2_2];
img2 = zeros(s2);
i = 1:s(1);
j = 1:s(2);
img2(i+mu,j+mu) = img(i,j);
[Gx,Gy] = imgradientxy(img2);
[Gmag,Gdir] = imgradient(Gx,Gy);
max_mag = max(Gmag(:));
Tx = -Gy/max_mag;
Ty = Gx/max_mag;
[Tmag,Tdir] = imgradient(Tx,Ty);
for p = 1:1
    Tx2 = zeros(size(Tx));
    Ty2 = zeros(size(Ty));
    for i = 1+mu:s(1)
        for j = 1+mu:s(2)
            sum_x = 0;
            sum_y = 0;
            for w_i = -mu:mu
                for w_j = -mu:mu
                    %ws = 1;
                    exponent = ((w_i*w_i) + (w_j*w_j))/(2*mu*mu);
                    ws = exp(-exponent);
                    wm = (Tmag(i+w_i,j+w_j) - Tmag(i,j) + 1)/2;
                    wd = Tx(i,j)*Tx(i+w_i,j+w_j) + Ty(i,j)*Ty(i+w_i,j+w_j);
                    sum_x = sum_x + Tx(i+w_i,j+w_j)*ws*wm*wd;
                    sum_y = sum_y + Ty(i+w_i,j+w_j)*ws*wm*wd;
                end
            end
            k = (2*mu+1)*(2*mu+1);
            Tx2(i,j) = sum_x;
            Ty2(i,j) = sum_y;
        end
    end
    Tx = Tx2;
    Ty = Ty2;
    Tmag = sqrt(Tx.*Tx+Ty.*Ty);
    max_mag = max(Tmag(:));
    Tmag = Tmag/max_mag;
    Tx = Tx/max_mag;
    Ty = Ty/max_mag;
end


    
    
