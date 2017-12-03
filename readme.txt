Brief Information about the codes.

myETF.m is code for Edge Tangent Flow used in myBFL2_color.m and myBFL2.m(Flow Based Bilateral Filter code for RGB and grayscale images respectively).
ETF2.m is code for Edge Tangent Flow used in LineExtraction4.m and LineExtraction2.m(we are using the former, but the latter can be also used).
combine.m calls myBFL2.m and LineExtraction4.m(or LineExtraction2.m) and displays the final cartoonified image.

How to run the codes.
1.Open the combine.m file.
2.In the 2nd line of this file write the name of the image which you want to convert to a cartoon.
3.In the 4th line you need to provide the 2 parameters sigma_m and sigma_c.
4.Now just run the combine.m to get the cartooned image.
5.If you just want to see the result of the Line Extraction or Flow Based Bilateral Filter then you must =>
 comment these lines :5,6,7,8,9,10,11 (for Line Extraction)
                     :4,6,7,8,9,10,11 (for Flow Based Bilateral Filter). 
 uncomment these lines : 12(for Line Extraction)
                         13(for FLow Based Bilateral Filter).
