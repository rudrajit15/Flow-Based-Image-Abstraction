function [img_out_color] = myBFL2_color(img,sigma_e,r_e,sigma_g,r_g)
%img = imread('boat.png');
S = 8;
T = 4;
max_ST = max(S,T);
img_out = rgb2gray(img);
sc = size(img);
sc2 = sc + [2*max_ST+2 2*max_ST+2 2*max_ST+2];
img2 = zeros(sc2);
i = 1:sc(1);
j = 1:sc(2);
img2(i+max_ST,j+max_ST,1) = img(i,j,1);
img2(i+max_ST,j+max_ST,2) = img(i,j,2);
img2(i+max_ST,j+max_ST,3) = img(i,j,3);
for p = 1:5
    s = size(img_out);
    [img_out2, Tx, Ty] = myETF(img_out,max_ST+1,0.2);
    [Tmag, Tdir] = imgradient(Tx,Ty);
    Gx = Ty;
    Gy = -Tx;
    [Gmag, Gdir] = imgradient(Gx,Gy);
    s2 = size(img_out2);
    if(mod(p,2) == 1)
        C_e2 = zeros(sc2);
        for i = 1:s2(1)-2*max_ST-1
            for j = 1:s2(2)-2*max_ST-1
                sum_u_R = 0;
                sum_u_G = 0;
                sum_u_B = 0;
                sum_n = 0;
                c_i = i;
                c_j = j;
                i2 = i;
                j2 = j;
                for k = 1:(2*S+1)
                    exponent_s = ((k-S-1)*(k-S-1))/(2*sigma_e*sigma_e);
                    g_s = exp(-exponent_s);
                    if(-22.5 <= Tdir(i2,j2) <= 22.5)
                        diff_r_R = img2(i2,j2,1) - img2(i2+1,j2,1);
                        diff_r_G = img2(i2,j2,2) - img2(i2+1,j2,2);
                        diff_r_B = img2(i2,j2,3) - img2(i2+1,j2,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_e*r_e);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i2+1,j2,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i2+1,j2,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i2+1,j2,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                        i2 = i2+1;

                    elseif(22.5 <= Tdir(i2,j2) <= 67.5)
                        diff_r_R = img2(i2,j2,1) - img2(i2+1,j2+1,1);
                        diff_r_G = img2(i2,j2,2) - img2(i2+1,j2+1,2);
                        diff_r_B = img2(i2,j2,3) - img2(i2+1,j2+1,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_e*r_e);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i2+1,j2+1,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i2+1,j2+1,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i2+1,j2+1,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                        i2 = i2+1;
                        j2 = j2+1;

                    elseif(67.5 <= Tdir(i2,j2) <= 112.5)
                        diff_r_R = img2(i2,j2,1) - img2(i2,j2+1,1);
                        diff_r_G = img2(i2,j2,2) - img2(i2,j2+1,2);
                        diff_r_B = img2(i2,j2,3) - img2(i2,j2+1,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_e*r_e);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i2,j2+1,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i2,j2+1,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i2,j2+1,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                        j2 = j2+1;

                    elseif(112.5 <= Tdir(i2,j2) <= 157.5)
                        diff_r_R = img2(i2,j2,1) - img2(i2-1,j2+1,1);
                        diff_r_G = img2(i2,j2,1) - img2(i2-1,j2+1,2);
                        diff_r_B = img2(i2,j2,1) - img2(i2-1,j2+1,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_e*r_e);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i2-1,j2+1,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i2-1,j2+1,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i2-1,j2+1,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                        i2 = i2-1;
                        j2 = j2+1;

                    elseif(Tdir(i2,j2) >= 157.5 || Tdir(i2,j2) <= -157.5)
                        diff_r_R = img2(i2,j2,1) - img2(i2-1,j2,1);
                        diff_r_G = img2(i2,j2,2) - img2(i2-1,j2,2);
                        diff_r_B = img2(i2,j2,3) - img2(i2-1,j2,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_e*r_e);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i2-1,j2,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i2-1,j2,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i2-1,j2,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                        i2 = i2-1;

                    elseif(-157.5 <= Tdir(i2,j2) <= -112.5)
                        diff_r_R = img2(i2,j2,1) - img2(i2-1,j2-1,1);
                        diff_r_G = img2(i2,j2,2) - img2(i2-1,j2-1,2);
                        diff_r_B = img2(i2,j2,3) - img2(i2-1,j2-1,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_e*r_e);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i2-1,j2-1,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i2-1,j2-1,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i2-1,j2-1,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                        i2 = i2-1;
                        j2 = j2-1;

                    elseif(-112.5 <= Tdir(i2,j2) <= -67.5)
                        diff_r_R = img2(i2,j2,1) - img2(i2,j2-1,1);
                        diff_r_G = img2(i2,j2,2) - img2(i2,j2-1,2);
                        diff_r_B = img2(i2,j2,3) - img2(i2,j2-1,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_e*r_e);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i2,j2-1,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i2,j2-1,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i2,j2-1,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                        j2 = j2-1;

                    elseif(-67.5 <= Tdir(i2,j2) <= -22.5)
                        diff_r_R = img2(i2,j2,1) - img2(i2+1,j2-1,1);
                        diff_r_G = img2(i2,j2,2) - img2(i2+1,j2-1,2);
                        diff_r_B = img2(i2,j2,3) - img2(i2+1,j2-1,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_e*r_e);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i2+1,j2-1,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i2+1,j2-1,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i2+1,j2-1,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                        i2 = i2+1;
                        j2 = j2-1;
                    end

                    if(k == S+1)
                        c_i = i2;
                        c_j = j2;
                    end
                end
                C_e2(c_i,c_j,1) = sum_u_R/sum_n;
                C_e2(c_i,c_j,2) = sum_u_G/sum_n;
                C_e2(c_i,c_j,3) = sum_u_B/sum_n;
            end
        end
        C_e = zeros(sc);
        i = 1:s(1);
        j = 1:s(2);
        C_e(i,j,1) = C_e2(i+max_ST+1,j+max_ST+1,1);
        C_e(i,j,2) = C_e2(i+max_ST+1,j+max_ST+1,2);
        C_e(i,j,3) = C_e2(i+max_ST+1,j+max_ST+1,3);
        img_out_color = C_e;
        
    else
        C_g2 = zeros(sc2);
        for i = 1+max_ST:s2(1)-max_ST-1
            for j = 1+max_ST:s2(2)-max_ST-1
                sum_u_R = 0;
                sum_u_G = 0;
                sum_u_B = 0;
                sum_n = 0;
                if(-22.5 <= Gdir(i,j) <= 22.5)
                    for q = -T:T
                        g_s = exp(-(q*q)/(2*sigma_g*sigma_g));
                        diff_r_R = img2(i+q,j,1) - img2(i,j,1);
                        diff_r_G = img2(i+q,j,2) - img2(i,j,2);
                        diff_r_B = img2(i+q,j,3) - img2(i,j,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_g*r_g);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i+q,j,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i+q,j,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i+q,j,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                    end

                elseif(22.5 <= Gdir(i,j) <= 67.5)
                    for q = -T:T
                        g_s = exp(-(q*q)/(2*sigma_g*sigma_g));
                        diff_r_R = img2(i+q,j+q,1) - img2(i,j,1);
                        diff_r_G = img2(i+q,j+q,2) - img2(i,j,2);
                        diff_r_B = img2(i+q,j+q,3) - img2(i,j,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_g*r_g);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i+q,j+q,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i+q,j+q,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i+q,j+q,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                    end

                elseif(67.5 <= Gdir(i,j) <= 112.5)
                    for q = -T:T
                        g_s = exp(-(q*q)/(2*sigma_g*sigma_g));
                        diff_r_R = img2(i,j+q,1) - img2(i,j,1);
                        diff_r_G = img2(i,j+q,2) - img2(i,j,2);
                        diff_r_B = img2(i,j+q,2) - img2(i,j,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_g*r_g);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i,j+q,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i,j+q,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i,j+q,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                    end

                elseif(112.5 <= Gdir(i,j) <= 157.5)
                    for q = -T:T
                        g_s = exp(-(q*q)/(2*sigma_g*sigma_g));
                        diff_r_R = img2(i-q,j+q,1) - img2(i,j,1);
                        diff_r_G = img2(i-q,j+q,2) - img2(i,j,2);
                        diff_r_B = img2(i-q,j+q,3) - img2(i,j,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_g*r_g);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i-q,j+q,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i-q,j+q,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i-q,j+q,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                    end

                elseif(Gdir(i,j) >= 157.5 || Gdir(i,j) <= -157.5)
                    for q = -T:T
                        g_s = exp(-(q*q)/(2*sigma_g*sigma_g));
                        diff_r_R = img2(i-q,j,1) - img2(i,j,1);
                        diff_r_G = img2(i-q,j,2) - img2(i,j,2);
                        diff_r_B = img2(i-q,j,3) - img2(i,j,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_g*r_g);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i-q,j,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i-q,j,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i-q,j,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                    end

                elseif(-157.5 <= Gdir(i,j) <= -112.5)
                    for q = -T:T
                        g_s = exp(-(q*q)/(2*sigma_g*sigma_g));
                        diff_r_R = img2(i-q,j-q,1) - img2(i,j,1);
                        diff_r_G = img2(i-q,j-q,2) - img2(i,j,2);
                        diff_r_B = img2(i-q,j-q,3) - img2(i,j,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_g*r_g);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i-q,j-q,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i-q,j-q,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i-q,j-q,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                    end

                elseif(-112.5 <= Gdir(i,j) <= -67.5)
                    for q = -T:T
                        g_s = exp(-(q*q)/(2*sigma_g*sigma_g));
                        diff_r_R = img2(i,j-q,1) - img2(i,j,1);
                        diff_r_G = img2(i,j-q,2) - img2(i,j,2);
                        diff_r_B = img2(i,j-q,3) - img2(i,j,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_g*r_g);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i,j-q,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i,j-q,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i,j-q,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                    end

                elseif(-67.5 <= Gdir(i,j) <= -22.5)
                    for q = -T:T
                        g_s = exp(-(q*q)/(2*sigma_g*sigma_g));
                        diff_r_R = img2(i+q,j-q,1) - img2(i,j,1);
                        diff_r_G = img2(i+q,j-q,2) - img2(i,j,2);
                        diff_r_B = img2(i+q,j-q,3) - img2(i,j,3);
                        diff = diff_r_R*diff_r_R + diff_r_G*diff_r_G + diff_r_B*diff_r_B;
                        exponent_r = (diff)/(2*r_g*r_g);
                        g_r = exp(-exponent_r);
                        sum_u_R = sum_u_R + img2(i+q,j-q,1)*g_s*g_r;
                        sum_u_G = sum_u_G + img2(i+q,j-q,2)*g_s*g_r;
                        sum_u_B = sum_u_B + img2(i+q,j-q,3)*g_s*g_r;
                        sum_n = sum_n + g_s*g_r;
                    end
                end
                C_g2(i,j,:) = sum_u_R/sum_n;
                C_g2(i,j,2) = sum_u_G/sum_n;
                C_g2(i,j,3) = sum_u_B/sum_n;
            end
        end
        C_g = zeros(sc);
        i = 1:s(1);
        j = 1:s(2);
        C_g(i,j,1) = C_g2(i+max_ST+1,j+max_ST+1,1);
        C_g(i,j,2) = C_g2(i+max_ST+1,j+max_ST+1,2);
        C_g(i,j,3) = C_g2(i+max_ST+1,j+max_ST+1,3);
        img_out_color = C_g;
    end
    img_out = rgb2gray(img_out_color);
end

%subplot(1,2,1),imshow(img);
%subplot(1,2,2),imshow(uint8(img_out_color));
%imwrite(uint8(img_out_color),'boat_BFL.png');













