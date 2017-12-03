function [img_out] = myBFL2(img,sigma_e,r_e,sigma_g,r_g)
%img = imread('4.2.04.tiff');
%img = rgb2gray(img);
img_out = img;
for p = 1:2
    s = size(img_out);
    %S = 3*sigma_e;
    %T = 3*sigma_g;
    S = 8;
    T = 4;
    max_ST = max(S,T);
    [img_out2, Tx, Ty] = myETF(img_out,max_ST+1,0.2);
    [Tmag, Tdir] = imgradient(Tx,Ty);
    Gx = Ty;
    Gy = -Tx;
    [Gmag, Gdir] = imgradient(Gx,Gy);
    s2 = size(img_out2);
    C_e2 = zeros(s2);
    C_g2 = zeros(s2);
    for i = 1:s2(1)-2*max_ST-1
        for j = 1:s2(2)-2*max_ST-1
     %for i = 1+max_ST:s2(1)-max_ST-1
      %  for j = 1+max_ST:s2(2)-max_ST-1
            sum_u = 0;
            sum_n = 0;
            c_i = i;
            c_j = j;
            i2 = i;
            j2 = j;
            for k = 1:(2*S+1)
                exponent_s = ((k-S-1)*(k-S-1))/(2*sigma_e*sigma_e);
                g_s = exp(-exponent_s);
                if(-22.5 <= Tdir(i2,j2) <= 22.5)
                    diff_r = img_out2(i2,j2) - img_out2(i2+1,j2);
                    exponent_r = (diff_r*diff_r)/(2*r_e*r_e);
                    g_r = exp(-exponent_r);
                    sum_u = sum_u + img_out2(i2+1,j2)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                    i2 = i2+1;

                elseif(22.5 <= Tdir(i2,j2) <= 67.5)
                    diff_r = img_out2(i2,j2) - img_out2(i2+1,j2+1);
                    exponent_r = (diff_r*diff_r)/(2*r_e*r_e);
                    g_r = exp(-exponent_r);
                    sum_u = sum_u + img_out2(i2+1,j2+1)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                    i2 = i2+1;
                    j2 = j2+1;

                elseif(67.5 <= Tdir(i2,j2) <= 112.5)
                    diff_r = img_out2(i2,j2) - img_out2(i2,j2+1);
                    exponent_r = (diff_r*diff_r)/(2*r_e*r_e);
                    g_r = exp(-exponent_r);
                    sum_u = sum_u + img_out2(i2,j2+1)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                    j2 = j2+1;

                elseif(112.5 <= Tdir(i2,j2) <= 157.5)
                    diff_r = img_out2(i2,j2) - img_out2(i2-1,j2+1);
                    exponent_r = (diff_r*diff_r)/(2*r_e*r_e);
                    g_r = exp(-exponent_r);
                    sum_u = sum_u + img_out2(i2-1,j2+1)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                    i2 = i2-1;
                    j2 = j2+1;

                elseif(Tdir(i2,j2) >= 157.5 || Tdir(i2,j2) <= -157.5)
                    diff_r = img_out2(i2,j2) - img_out2(i2-1,j2);
                    exponent_r = (diff_r*diff_r)/(2*r_e*r_e);
                    g_r = exp(-exponent_r);
                    sum_u = sum_u + img_out2(i2-1,j2)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                    i2 = i2-1;

                elseif(-157.5 <= Tdir(i2,j2) <= -112.5)
                    diff_r = img_out2(i2,j2) - img_out2(i2-1,j2-1);
                    exponent_r = (diff_r*diff_r)/(2*r_e*r_e);
                    g_r = exp(-exponent_r);
                    sum_u = sum_u + img_out2(i2-1,j2-1)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                    i2 = i2-1;
                    j2 = j2-1;

                elseif(-112.5 <= Tdir(i2,j2) <= -67.5)
                    diff_r = img_out2(i2,j2) - img_out2(i2,j2-1);
                    exponent_r = (diff_r*diff_r)/(2*r_e*r_e);
                    g_r = exp(-exponent_r);
                    sum_u = sum_u + img_out2(i2,j2-1)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                    j2 = j2-1;

                elseif(-67.5 <= Tdir(i2,j2) <= -22.5)
                    diff_r = img_out2(i2,j2) - img_out2(i2+1,j2-1);
                    exponent_r = (diff_r*diff_r)/(2*r_e*r_e);
                    g_r = exp(-exponent_r);
                    sum_u = sum_u + img_out2(i2+1,j2-1)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                    i2 = i2+1;
                    j2 = j2-1;
                end

                if(k == S+1)
                    c_i = i2;
                    c_j = j2;
                end
            end
            C_e2(c_i,c_j) = sum_u/sum_n;
        end
    end

    for i = 1+max_ST:s2(1)-max_ST-1
        for j = 1+max_ST:s2(2)-max_ST-1
            sum_u = 0;
            sum_n = 0;
            if(-22.5 <= Gdir(i,j) <= 22.5)
                for q = -T:T
                    g_s = normpdf(q,0,sigma_g);
                    diff_r = img_out2(i+q,j) - img_out2(i,j);
                    g_r = normpdf(diff_r,0,r_g);
                    sum_u = sum_u + img_out2(i+q,j)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                end

            elseif(22.5 <= Gdir(i,j) <= 67.5)
                for q = -T:T
                    g_s = normpdf(q,0,sigma_g);
                    diff_r = img_out2(i+q,j+q) - img_out2(i,j);
                    g_r = normpdf(diff_r,0,r_g);
                    sum_u = sum_u + img_out2(i+q,j+q)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                end

            elseif(67.5 <= Gdir(i,j) <= 112.5)
                for q = -T:T
                    g_s = normpdf(q,0,sigma_g);
                    diff_r = img_out2(i,j+q) - img_out2(i,j);
                    g_r = normpdf(diff_r,0,r_g);
                    sum_u = sum_u + img_out2(i,j+q)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                end

            elseif(112.5 <= Gdir(i,j) <= 157.5)
                for q = -T:T
                    g_s = normpdf(q,0,sigma_g);
                    diff_r = img_out2(i-q,j+q) - img_out2(i,j);
                    g_r = normpdf(diff_r,0,r_g);
                    sum_u = sum_u + img_out2(i-q,j+q)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                end
                
            elseif(Gdir(i,j) >= 157.5 || Gdir(i,j) <= -157.5)
                for q = -T:T
                    g_s = normpdf(q,0,sigma_g);
                    diff_r = img_out2(i-q,j) - img_out2(i,j);
                    g_r = normpdf(diff_r,0,r_g);
                    sum_u = sum_u + img_out2(i-q,j)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                end

            elseif(-157.5 <= Gdir(i,j) <= -112.5)
                for q = -T:T
                    g_s = normpdf(q,0,sigma_g);
                    diff_r = img_out2(i-q,j-q) - img_out2(i,j);
                    g_r = normpdf(diff_r,0,r_g);
                    sum_u = sum_u + img_out2(i-q,j-q)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                end

            elseif(-112.5 <= Gdir(i,j) <= -67.5)
                for q = -T:T
                    g_s = normpdf(q,0,sigma_g);
                    diff_r = img_out2(i,j-q) - img_out2(i,j);
                    g_r = normpdf(diff_r,0,r_g);
                    sum_u = sum_u + img_out2(i,j-q)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                end

            elseif(-67.5 <= Gdir(i,j) <= -22.5)
                for q = -T:T
                    g_s = normpdf(q,0,sigma_g);
                    diff_r = img_out2(i+q,j-q) - img_out2(i,j);
                    g_r = normpdf(diff_r,0,r_g);
                    sum_u = sum_u + img_out2(i+q,j-q)*g_s*g_r;
                    sum_n = sum_n + g_s*g_r;
                end
            end
            C_g2(i,j) = sum_u/sum_n;
        end
    end
    C_e = zeros(s);
    C_g = zeros(s);
    i = 1:s(1);
    j = 1:s(2);
    C_e(i,j) = C_e2(i+max_ST+1,j+max_ST+1);
    C_g(i,j) = C_g2(i+max_ST+1,j+max_ST+1);
    
    if(mod(p,2) == 1)
        C = C_e;
    else
        C = C_g;
    end
    img_out = C;
end

%imwrite(uint8(img_out),'FlowBasedBilateralFilter.jpg');














