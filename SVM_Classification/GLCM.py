# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 22:49:30 2018

@author: User
"""
from osgeo import gdal #, gdalconst
import numpy as np
import matplotlib.pyplot as plt

#from scipy import stats
from skimage import io
from skimage import data
from skimage.feature import greycomatrix, greycoprops

filename = "C:\\Users\\User\\Desktop\\Reaserch Data\\second\\SVM-for-road-network-extraction-from-HRRSI-master\\Text_2\\Ortho_19.tif"
outFileName = 'C:\\Users\\User\\Desktop\\Reaserch Data\\second\\SVM-for-road-network-extraction-from-HRRSI-master\\ASM_45_copy.tif'

img = gdal.Open(filename, gdal.GA_ReadOnly)
band = img.GetRasterBand(1)
geo_transform = img.GetGeoTransform()
proj = img.GetProjectionRef()
array = np.array(band.ReadAsArray())
[rows, cols] = array.shape  

nan = np.isnan(array)
np.where(np.isnan(array))
np.nan_to_num(array)


nodata = band.GetNoDataValue()

#array = array.astype((np.float)

#array[(array < 1) | (array > 4095)] = 0

where_are_NaNs = np.isnan(array)
array[where_are_NaNs] = nodata

if nodata is not None:
    array = np.ma.masked_equal(array, nodata)

print('min: %s, max: %s, mean: %s, std: %s' %
      (array.min(), array.max(), array.mean(), array.std()))


array_min= array.min()
array_max = array.max()
array_mean = array.mean()


img_scaled = np.interp(array, (array.min(), array.max()), (1, 64))

print('min: %s, max: %s, mean: %s, std: %s' %
      (img_scaled.min(), img_scaled.max(), img_scaled.mean(), img_scaled.std()))



img_scaled = img_scaled.astype((np.uint8))
def itemfreq(a): 
    items, inv = np.unique(a, return_inverse=True) 
    freq = np.bincount(inv) 
    return np.array([items, freq]).T 

def offset(length, angle):
    """Return the offset in pixels for a given length and angle"""
    dv = length * np.sign(-np.sin(angle)).astype(np.int32)
    dh = length * np.sign(np.cos(angle)).astype(np.int32)
    return dv, dh

def crop(img, center, win):
    """Return a square crop of img centered at center (side = 2*win + 1)"""
    row, col = center
    side = 2*win + 1
    first_row = row - win
    first_col = col - win
    last_row = first_row + side    
    last_col = first_col + side
    return img[first_row: last_row, first_col: last_col]

def cooc_maps(img, center, win, d=[1], theta=[0], levels=256):
    """
    Return a set of co-occurrence maps for different d and theta in a square 
    crop centered at center (side = 2*w + 1)
    """
    shape = (2*win + 1, 2*win + 1, len(d), len(theta))
    cooc = np.zeros(shape=shape, dtype=np.int32)
    row, col = center
    Ii = crop(img, (row, col), win)
    for d_index, length in enumerate(d):
        for a_index, angle in enumerate(theta):
            dv, dh = offset(length, angle)
            Ij = crop(img, center=(row + dv, col + dh), win=win)
            cooc[:, :, d_index, a_index] = encode_cooccurrence(Ii, Ij, levels)
    return cooc

def encode_cooccurrence(x, y, levels=256):
    """Return the code corresponding to co-occurrence of intensities x and y"""
    return x*levels + y

def decode_cooccurrence(code, levels=256):
    """Return the intensities x, y corresponding to code"""
    return code//levels, np.mod(code, levels)    

def compute_glcms(cooccurrence_maps, levels=256):
    """Compute the cooccurrence frequencies of the cooccurrence maps"""
    Nr, Na = cooccurrence_maps.shape[2:]
    glcms = np.zeros(shape=(levels, levels, Nr, Na), dtype=np.float64)
    for r in range(Nr):
        for a in range(Na):
            table = itemfreq(cooccurrence_maps[:, :, r, a])
           # table = np.unique(cooccurrence_maps[:, :, r, a],return_counts=True)
            codes = table[:, 0]
            freqs = table[:, 1]/float(table[:, 1].sum())
            i, j = decode_cooccurrence(codes, levels=levels)
            glcms[i, j, r, a] = freqs
    return glcms

def compute_props(glcms, props=('contrast',)):
    """Return a feature vector corresponding to a set of GLCM"""
    Nr, Na = glcms.shape[2:]
    features = np.zeros(shape=(Nr, Na, len(props)))
    for index, prop_name in enumerate(props):
        features[:, :, index] = greycoprops(glcms, prop_name)
    return features.ravel()

def haralick_features(img, win, d, theta, levels, props):
    """Return a map of Haralick features (one feature vector per pixel)"""
    rows, cols = img.shape
    margin = win + max(d)
    arr = np.pad(img, margin, mode='reflect')
    n_features = len(d) * len(theta) * len(props)
    feature_map = np.zeros(shape=(rows, cols, n_features), dtype=np.float64)
    for m in xrange(rows):
        for n in xrange(cols):
            coocs = cooc_maps(arr, (m + margin, n + margin), win, d, theta, levels)
            glcms = compute_glcms(coocs, levels)
            feature_map[m, n, :] = compute_props(glcms, props)
    return feature_map


d = (1, 2)

theta = (0, np.pi/4, np.pi/2, 3*np.pi/4)

#props = ('contrast', 'homogeneity', 'energy')

props = ('contrast','homogeneity', 'energy')

levels = 64

win = 2

feature_map = haralick_features(img_scaled, win, d, theta, levels, props)


outFile_dir = 'C:\\Users\\User\\Desktop\\Reaserch Data\\second\\SVM-for-road-network-extraction-from-HRRSI-master\\'
for i in range(0, 24):
    filename = "GLCM_" + str(i) + ".tif"
    driver = gdal.GetDriverByName("GTiff")
    out_data = driver.Create(outFile_dir + filename, cols, rows, 1, gdal.GDT_Float64)
    out_data.SetGeoTransform(img.GetGeoTransform())
    out_data.SetProjection(img.GetProjection())
    out_data.GetRasterBand(1).WriteArray(feature_map[:, :, i])
    out_data.GetRasterBand(1).SetNoDataValue(nodata)
    out_data.FlushCache()  # saves to disk!!
    out_data = None

