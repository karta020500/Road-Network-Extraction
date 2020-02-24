# Road-Network-Extraction
This project develops a systematic procedure to extract roads in different scales depending on the resolution of images. The developed procedure to extract roads and reconstruct the network can be divided into two steps. The first step, obtaining texture parameters including contrast, entropy and homogeneity of each pixel using gray level co-occurrence matrix (GLCM); afterwards, support vector machine (SVM) is used as a classifier to identify road pixels. Finally, linear feature detection of the road network with radon transform is performed to trace road pixel.

# Installation
This project has the following dependencies:

-*Numpy pip install numpy*

-*OpenCV Python apt-get install python-opencv*

-*GDAL conda install --GDAL*

-*Scikit-learn pip install -U scikit-learn*

# Results

|  <img src="https://github.com/karta020500/Road-Network-Extraction/blob/master/SVM_Classification/Results/classification.png" width = "500" height = "330" /> | <img src="https://github.com/karta020500/Road-Network-Extraction/blob/master/Radon_line_tracing/Results/final.jpg" width = "500" height = "330" />  | 
|:-------:|:-----:|
|Classification|Road network completeness|

# Conclusions
1. This study presents an effective road extraction technique from high-resolution imagery in urban areas.
2. The estimated classification achieves accuracy higher than 90% and kappa coefficient over 85%.
3. Road network evaluation shows that the completeness higher than 85% and RMSD of 4.82 meters. 
4. The main omission seems to happen in areas with strong illumination changes or occlusion.


