# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#import arcpy as arc
#import pandas as pd
from osgeo import gdal #, gdalconst
import os
#import cv2
import numpy as np
import matplotlib.pyplot as plt

#from scipy import stats
from skimage import io
from skimage import data
from sklearn import metrics
#from skimage.feature import greycomatrix, greycoprops
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC, LinearSVC
from sklearn import preprocessing
from skimage.morphology import erosion, dilation, opening, closing
from skimage.morphology import disk
#from sklearn.cluster import KMeans
#import matplotlib.image as mpimg

file_dir = '/Users/nick/Desktop/Coding/Machine Learning/Test_two/Input/GLCM/'
outFileName = '/Users/nick/Desktop/Coding/Machine Learning/Test_two/Output/test.tif'
array = []
count = 0
for i in range(0, 12):
    filename = "GLCM_" + str(i) + ".tif"
    img = gdal.Open(file_dir + filename, gdal.GA_ReadOnly)
    band = img.GetRasterBand(1)
    array.append(band.ReadAsArray())
    count = count+1

'''
geo_transform_GL = img.GetGeoTransform()
proj_GL = img.GetProjectionRef()
bands_data_GL = []
for b in range(1, img.RasterCount+1):
    band = img.GetRasterBand(b)
    bands_data_GL.append(band.ReadAsArray())
'''
    
'''
array = np.dstack(bands_data_GL)
array = np.delete(array, (0), axis=0)
array = np.delete(array, (0), axis=1)
'''
'''
nan = np.isnan(array)
np.where(np.isnan(array))
np.nan_to_num(array)
'''
'''
def build_band_stack(image_dataset, num_bands):
    band_list = []
    for band in range(num_bands):
        b = image_dataset.GetRasterBand(band+1)
        b_array = b.ReadAsArray(0, 0)
        band_list.append(b_array)

    array = np.dstack(band_list)

    return array

build_band_stack
'''
'''
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


#img_scaled = cv.normalize(array, dst=None, alpha=1, beta=256, norm_type=cv.NORM_MINMAX)
'''
'''
#arr_out = np.where((arr < arr_mean), 65535, arr)
driver = gdal.GetDriverByName("GTiff")
out_data = driver.Create(outFileName, rows, cols, 1, gdal.GDT_UInt16)
#out_data.SetGeoTransform(array.GetGeoTransform())
#out_data.SetProjection(array.GetProjection())  # sets same projection as input
out_data.GetRasterBand(1).WriteArray(array)
out_data.GetRasterBand(1).SetNoDataValue(65535)  # if you want these values transparent
out_data.FlushCache()  # saves to disk!!
out_data = None
band = None
ds = None
'''

#plt.imshow(img_scaled, cmap="gray", vmin=1, vmax=256)
'''
plt.imshow(array, cmap="RGB", vmin=1, vmax=4095)
plt.show()
'''

'''

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

def encode_cooccurrence(x, y, levels=64):
    """Return the code corresponding to co-occurrence of intensities x and y"""
    return x*levels + y

def decode_cooccurrence(code, levels=64):
    """Return the intensities x, y corresponding to code"""
    return code//levels, np.mod(code, levels)    

def compute_glcms(cooccurrence_maps, levels=64):
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


d = (1,2)

theta = (0, np.pi/4, np.pi/2, 3*np.pi/4)

props = ('contrast', 'homogeneity', 'energy')

levels = 64

win = 2
'''
#feature_map = haralick_features(img_scaled, win, d, theta, levels, props)

'''
distances = [1, 2, 3]
angles = [0, np.pi/4, np.pi/2, 3*np.pi/4]
properties = ['energy', 'homogeneity']

glcm = greycomatrix(img_scaled, 
                    distances=distances, 
                    angles=angles,
                    symmetric=True,
                    normed=True)
'''
#----------------------------------------------------------------------------
'''
feats = np.hstack([greycoprops(glcm, prop).ravel() for prop in properties])
#sarraster = band.ReadAsArray()
#  sarraster is satellite image, testraster will receive texture
testraster = np.copy(img_scaled)
testraster[:] = 0

for i in range(testraster.shape[0]):
    print i,
    for j in range(testraster.shape[1]):

        #  windows needs to fit completely in image
        if i < 3 or j < 3:
            continue
        if i > (testraster.shape[0] - 4) or j > (testraster.shape[0] - 4):
            continue

        #  Calculate GLCM on a 7x7 window
        glcm_window = img_scaled[i-3: i+4, j-3: j+4]
        glcm = greycomatrix(glcm_window, [1], [0],  symmetric=True, normed=True)

        #  Calculate contrast and replace center pixel
        contrast = greycoprops(glcm, 'contrast')
        testraster[i, j] = contrast
sarplot = plt.imshow(testraster, cmap='gray', vmin=0, vmax=255)
'''

#-----------------------------------------------------------------------------
'''
X = img_scaled.reshape(rows*cols, 1)
classes = {'vegetation': 0, 'building': 1, 'water': 2, 'road': 3}
n_classes = len(classes)
palette = np.uint8([[0, 255, 0], [255, 0, 0], [0, 0, 255], [0, 0, 0]])
'''
'''
supervised = n_classes*np.ones(shape=(rows, cols), dtype=np.int)

supervised[200:220, 150:170] = classes['building']
supervised[650:700, 550:600] = classes['vegetation']
supervised[750:800, 550:600] = classes['water']
supervised[750:800, 70:100]  = classes['road']

y = supervised.ravel()
train = np.flatnonzero(supervised < 4)
test = np.flatnonzero(supervised == 4)


clf = SVC(kernel='rbf')
clf.fit(X[train], y[train])
y[test] = clf.predict(X[test])
supervised = y.reshape(rows, cols)
io.imshow(palette[supervised])
'''
'''
kmeans = KMeans(n_clusters=n_classes, random_state=0).fit(X)
unsupervised = kmeans.labels_.reshape(rows, cols)
io.imshow(palette[unsupervised])
'''
#-----------------------------------------------------------------------------
COLORS = [
    "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
    "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
    "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
    "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
    "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"
]


def create_mask_from_vector(vector_data_path, cols, rows, geo_transform, projection, target_value=1,
                            output_fname='', dataset_format='MEM'):
    """
    Rasterize the given vector (wrapper for gdal.RasterizeLayer). Return a gdal.Dataset.
    :param vector_data_path: Path to a shapefile
    :param cols: Number of columns of the result
    :param rows: Number of rows of the result
    :param geo_transform: Returned value of gdal.Dataset.GetGeoTransform (coefficients for
                          transforming between pixel/line (P,L) raster space, and projection
                          coordinates (Xp,Yp) space.
    :param projection: Projection definition string (Returned by gdal.Dataset.GetProjectionRef)
    :param target_value: Pixel value for the pixels. Must be a valid gdal.GDT_UInt16 value.
    :param output_fname: If the dataset_format is GeoTIFF, this is the output file name
    :param dataset_format: The gdal.Dataset driver name. [default: MEM]
    """
    
    data_source = gdal.OpenEx(vector_data_path, gdal.OF_VECTOR)
    if data_source is None:
        print("File read failed: %s", vector_data_path)
    layer = data_source.GetLayer(0)
    driver = gdal.GetDriverByName(dataset_format)
    target_ds = driver.Create(output_fname, cols, rows, 1, gdal.GDT_UInt16)
    target_ds.SetGeoTransform(geo_transform)
    target_ds.SetProjection(projection)
    gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[target_value])
    return target_ds


def vectors_to_raster(file_paths, rows, cols, geo_transform, projection):
    """
    Rasterize, in a single image, all the vectors in the given directory.
    The data of each file will be assigned the same pixel value. This value is defined by the order
    of the file in file_paths, starting with 1: so the points/poligons/etc in the same file will be
    marked as 1, those in the second file will be 2, and so on.
    :param file_paths: Path to a directory with shapefiles
    :param rows: Number of rows of the result
    :param cols: Number of columns of the result
    :param geo_transform: Returned value of gdal.Dataset.GetGeoTransform (coefficients for
                          transforming between pixel/line (P,L) raster space, and projection
                          coordinates (Xp,Yp) space.
    :param projection: Projection definition string (Returned by gdal.Dataset.GetProjectionRef)
    """
    labeled_pixels = np.zeros((rows, cols))
    for i, path in enumerate(file_paths):
        label = i+1
        print("Processing file %s: label (pixel value) %i", path, label)
        ds = create_mask_from_vector(path, cols, rows, geo_transform, projection,
                                     target_value=label)
        band = ds.GetRasterBand(1)
        a = band.ReadAsArray()
        print("Labeled pixels: %i", len(a.nonzero()[0]))
        labeled_pixels += a
        ds = None
    return labeled_pixels


def write_geotiff(fname, data, geo_transform, projection, data_type=gdal.GDT_Byte):
    """
    Create a GeoTIFF file with the given data.
    :param fname: Path to a directory with shapefiles
    :param data: Number of rows of the result
    :param geo_transform: Returned value of gdal.Dataset.GetGeoTransform (coefficients for
                          transforming between pixel/line (P,L) raster space, and projection
                          coordinates (Xp,Yp) space.
    :param projection: Projection definition string (Returned by gdal.Dataset.GetProjectionRef)
    """
    driver = gdal.GetDriverByName('GTiff')
    rows, cols = data.shape
    dataset = driver.Create(fname, cols, rows, 1, data_type)
    dataset.SetGeoTransform(geo_transform)
    dataset.SetProjection(projection)
    band = dataset.GetRasterBand(1)
    band.WriteArray(data)

    ct = gdal.ColorTable()
    for pixel_value in range(len(classes)+1):
        color_hex = COLORS[pixel_value]
        r = int(color_hex[1:3], 16)
        g = int(color_hex[3:5], 16)
        b = int(color_hex[5:7], 16)
        ct.SetColorEntry(pixel_value, (r, g, b, 255))
    band.SetColorTable(ct)

    metadata = {
        'TIFFTAG_COPYRIGHT': 'CC BY 4.0',
        'TIFFTAG_DOCUMENTNAME': 'classification',
        'TIFFTAG_IMAGEDESCRIPTION': 'Supervised classification.',
        'TIFFTAG_MAXSAMPLEVALUE': str(len(classes)),
        'TIFFTAG_MINSAMPLEVALUE': '0',
        'TIFFTAG_SOFTWARE': 'Python, GDAL, scikit-learn'
    }
    dataset.SetMetadata(metadata)

    dataset = None  # Close the file
    return





raster_data_path =  '/Users/nick/Desktop/Coding/Machine Learning/Test_two/Input/Ortho_19.tif'
train_data_path = '/Users/nick/Desktop/Coding/Machine Learning/Test_two/Input/Traning_1'
validation_data_path = '/Users/nick/Desktop/Coding/Machine Learning/Test_two/Input/Validate'
output_fname = '/Users/nick/Desktop/Coding/Machine Learning/Test_two/Output/try.tif'
#method = "random-forest"
method = "svm"


gdal.UseExceptions()

#logger.debug("Reading the input: %s", raster_data_path)
raster_dataset = gdal.Open(raster_data_path, gdal.GA_ReadOnly)
#except RuntimeError as e:
#    report_and_exit(str(e))


geo_transform = raster_dataset.GetGeoTransform()
proj = raster_dataset.GetProjectionRef()
bands_data = []
#source image
for b in range(1, raster_dataset.RasterCount+1):
    band = raster_dataset.GetRasterBand(b)
    bands_data.append(band.ReadAsArray())

for c in range(0, count):    
    bands_data.append(array[c])


bands_data = np.dstack(bands_data)
nodata = band.GetNoDataValue()


'''
where_are_NaNs = np.isnan(bands_data)
array[where_are_NaNs] = 65535
'''



if nodata is not None:
    bands_data = np.ma.masked_equal(bands_data, nodata)
    
rows, cols, n_bands = bands_data.shape
# A sample is a vector with all the bands data. Each pixel (independent of its position) is a
# sample.
n_samples = rows*cols



#logger.debug("Process the training data")
files = [f for f in os.listdir(train_data_path) if f.endswith('.shp')]
classes = [f.split('.')[0] for f in files]
shapefiles = [os.path.join(train_data_path, f) for f in files if f.endswith('.shp')]

labeled_pixels = vectors_to_raster(shapefiles, rows, cols, geo_transform, proj)

is_train = np.nonzero(labeled_pixels)

training_labels = labeled_pixels[is_train]
training_samples = bands_data[is_train]
training_samples_scaled = preprocessing.scale(training_samples)
training_samples_normalized = preprocessing.normalize(training_samples, norm='l2')

flat_pixels = bands_data.reshape((n_samples, n_bands))

#
# Perform classification
#
CLASSIFIERS = {
    # http://scikit-learn.org/dev/modules/generated/sklearn.ensemble.RandomForestClassifier.html
    "random-forest": RandomForestClassifier(n_jobs= 2, n_estimators=500, class_weight='balanced' ),
    # http://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html
    "svm" : SVC(C=10, gamma=0.000003)
}

classifier = CLASSIFIERS[method]
#logger.debug("Train the classifier: %s", str(classifier))
classifier.fit(training_samples, training_labels).score(training_samples, training_labels)

#logger.debug("Classifing...")
result = classifier.predict(flat_pixels)
#classifier.fit(training_samples, training_labels).score(flat_pixels, result)
# Reshape the result: split the labeled pixels into rows to create an image
classification = result.reshape((rows, cols))
write_geotiff(output_fname, classification, geo_transform, proj)
io.imshow(classification)
#plt.scatter(training_samples[:,0], training_samples[:,1], c=training_labels)
#plt.show

'''
plt.scatter(training_samples_scaled[:,2], training_samples_scaled[:,3], c=training_labels)
'''
if validation_data_path:
    
    shapefiles = [os.path.join(validation_data_path, "%s.shp" % c) for c in classes]
    verification_pixels = vectors_to_raster(shapefiles, rows, cols, geo_transform, proj)
    for_verification = np.nonzero(verification_pixels)
    verification_labels = verification_pixels[for_verification]
    predicted_labels = classification[for_verification]
    metrics_1= metrics.confusion_matrix(verification_labels, predicted_labels)
    target_names = ['Class %s' % s for s in classes]
   
    metrics_2 = metrics.classification_report(verification_labels, predicted_labels,
                                              target_names=target_names)
    
    metrics_3 = metrics.accuracy_score(verification_labels, predicted_labels)
#------------------------------------------------------------------------------