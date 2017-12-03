function [] = combine()
img = imread('../data/house.tiff');
img_gray = rgb2gray(img);
img1 = LineExtraction4(img_gray,0.25,0.7);
img2 = myBFL2_color(img,2,10,0.3,10);
img3(:,:,1) = img1;
img3(:,:,2) = img1;
img3(:,:,3) = img1;
imgr = (img3.*img2)/255;
imshow(uint8(imgr));
imwrite(uint8(imgr),'../images/house_cartoon.jpg');
%imshow(img1);
%imshow(uint8(img2));
%subplot(1,2,1),imshow(img);
%subplot(1,2,2),imshow(uint8(imgr));