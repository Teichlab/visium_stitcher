"""Stitch multiple Visium slides together"""
__version__ = "0.0.3"

from matplotlib.image import imread
from scipy.spatial import cKDTree
import xmltodict
import anndata as ad
import numpy as np
import cv2

def transform_finder(xmlpath, field="@file_path"):
    '''
    Load affine transformation matrices from a TrakEM2 XML file. Returns a dictionary 
    with the values of the specified field (to help determine which sample/image the 
    transform corresponds to) as the keys, and their corresponding transforms as the 
    values. Store the matrices in  ``.uns['transform']`` of the appropriate slide's 
    ``AnnData`` prior to calling ``visium_stitcher.stitch()``.
    
    Input
    -----
    xmlpath : ``filename``
        XML output of aligning the Visium hires images within TrakEM2.
    field : ``str``, optional (default: ``"@file_path"``)
        The XML field to be used to label the individual layers in the transform 
        dictionary. The default is ``"@file_path"``, the path to each input image provided 
        to TrakEM2. Another possibility is ``"@title"``, storing the layer name.
    '''
    #get the XML
    with open(xmlpath, "r") as fid:
        #this gobbles the whole file, which is what xmltodict wants on input
        xml = xmltodict.parse(fid.read())
    transforms = {}
    #this is where the various images reside in the XML structure
    for entry in xml['trakem2']['t2_layer_set']['t2_layer']['t2_patch']:
        #get out the transform, then turn it from matrix() to [] around the numbers
        transform_list = eval(entry["@transform"].replace("matrix(","[").replace(")","]"))
        #turn the list to a 2D numpy array, placing the elements in rows first (order F)
        transform_matrix = np.array(transform_list).reshape((2,3),order="F")
        #stash based on the desired field to name on
        transforms[entry[field]] = transform_matrix
    return transforms

def stitch(adatas, image=None, dist_fact = 1.5):
    '''
    Perform a series of affine transformations on individual Visium sample images/spot 
    coordinates, then merge them into a single multi-sample ``AnnData``. The output object 
    has the merged ``hires`` image stored in ``.uns['spatial']['joint']``, with the scale 
    factors factored into each sample's spot coordinates and set to 1 in the object. The 
    spot diameter is derived from the first object, and multiplied by the corresponding 
    scale factor.
    
    Input
    -----
    adatas : list of ``AnnData``
        Must have a 2x3 affine transformation matrix, as defined `here <https://theailearner.com/tag/cv2-warpaffine/>`_, 
        stored in ``.uns['transform']``. The raw feature count matrix should be loaded for 
        accurate overlap evaluation. The samples will be given image/spot priority based 
        on the order of the objects in the list.
    image : ``filepath``, optional (Default: ``None``)
        Optional path to previously prepared overlapped tissue image to use in the merged 
        object.
    dist_fact : ``float``, optional (Default: 1.5)
        For the first object in the list, ``med_dist`` is defined as the median of the 
        distances between each spot and its closest neighbour. For subsequent samples, 
        spots that fall within ``dist_fact*med_dist`` of a prior non-overlapping spot are 
        deemed to be overlap.
    '''
    #STEP ONE - transform the spots, determine overlaps based on sample order, infer canvas size
    adata = None
    canvas = np.array([0,0])
    main_spots = None
    spot_dist = None

    #start with a go through for coordinates/canvas size stuff
    for obj in adatas:
        #copy the object to be able to muck around on it in peace
        sample = obj.copy()
        #extract the scale factor for the hires image, and get scaled coordinates of the spots
        scalef = list(sample.uns['spatial'].values())[0]['scalefactors']["tissue_hires_scalef"]
        spatial_scaled = sample.obsm["spatial"]*scalef
        #transform spots, matching cv2.warpAffine() logic
        #https://theailearner.com/tag/cv2-warpaffine/
        #add an extra column of ones, and transpose the coordinates to match the example
        #spatial[:,0] is columns(x) and spatial[:,1] is rows(y)
        spatial_scaled_ones = np.hstack((spatial_scaled, np.ones((spatial_scaled.shape[0],1)))).T
        #transform the coordinates, then transpose them to match obsm orientation
        sample.obsm['spatial'] = np.matmul(sample.uns['transform'], spatial_scaled_ones).T
        #perform a second transformation - just the corners of the image
        #to get a feel for canvas size for image merging later
        shape = list(sample.uns['spatial'].values())[0]['images']['hires'].shape[:2]
        #the four corners are defined like so (and then transposed)
        #shape[0] is rows(y) and shape[1] is columns(x), so have to reverse them
        #to match the [x,y,1] orientation in the example
        dummy = np.array([[0,0,1],[shape[1],0,1],[0,shape[0],1],[shape[1],shape[0],1]]).T
        #transform the corners, then transpose to match obsm orientation again
        edges = np.matmul(sample.uns['transform'], dummy).T
        #get the maximum encountered coordinates, doing a ceiling of the spot locations
        canvas = np.maximum(canvas, np.max(np.ceil(edges), axis=0))
        #spot overlap, then stick the objects together
        if adata is None:
            #no prior spots. we are the starting main spots
            main_spots = sample.obsm['spatial'].copy()
            #get the distance threshold to use for subsequent overlap computations
            ckd = cKDTree(main_spots)
            ckdout = ckd.query(x=main_spots, k=2, workers=-1)
            #the distances are the first element of this tuple
            #we want the median distance to the second neighbour, as each spot is going to find itself as the closest
            med_dist = np.median(ckdout[0][:,1])
            spot_dist = dist_fact * med_dist
            #stash single sample object as starting point. we have no overlaps here
            sample.obs["overlap"] = False
            adata = sample.copy()
        else:
            #get the nearest neighbour within the existing main spots for each spot in this sample
            query_spots = sample.obsm['spatial']
            ckd = cKDTree(main_spots)
            ckdout = ckd.query(x=query_spots, k=1, workers=-1)
            #spots that are closer than the spot_dist are an overlap
            sample.obs["overlap"] = ckdout[0] < spot_dist
            #spots that are further than the spot_dist are new main spots
            main_spots = np.vstack((main_spots, query_spots[ckdout[0] >= spot_dist,:]))
            #join objects
            adata = ad.concat([adata, sample])
    #STEP TWO - transform the image and paste it together, using sample order as priority
    #well, unless the user just gives us an image on input.
    if image is not None:
        img = imread(image)
    else:
        img = None
        #the identified canvas size needs to be turned to integers
        #warpAffine expects it in cols(x),rows(y) order, which we already have
        img_size = canvas.astype(int)

        #image stuff
        for obj in adatas:
            #mask the non-spot area with zeroes. extract the image
            simg = list(obj.uns['spatial'].values())[0]['images']['hires'].copy()
            #get the size factor
            sf = list(obj.uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef']
            #get the radius of the spots - there's a diameter
            #halve it and scale it by the size factor
            sr = np.ceil(list(obj.uns['spatial'].values())[0]['scalefactors']['spot_diameter_fullres']*0.5*sf).astype(int)
            #get the box defined by the minimum/maximum coordinates of the spots
            #(don't forget to reverse as spots are x,y while the image array is y,x)
            #these are pre-transformation, multiply them by the size factor
            #and then account for the spot radius
            mins = np.floor(np.min(obj.obsm['spatial'], axis=0)[::-1] * sf).astype(int) - sr
            maxes = np.ceil(np.max(obj.obsm['spatial'], axis=0)[::-1] * sf).astype(int) + sr
            #now we can mask the image
            simg[:mins[0], :, :] = 0
            simg[:, :mins[1], :] = 0
            simg[maxes[0]:, :, :] = 0
            simg[:, maxes[1]:, :] = 0
            #transform the image
            imgt = cv2.warpAffine(simg,obj.uns['transform'], img_size)
            #fill in all-black (empty) areas on the glued together image
            if img is None:
                img = imgt.copy()
            else:
                #this masks the entries in the array based on the sum of the colour axis
                #it gets all the coordinates right, and then gets all the colour values for them
                img[np.sum(img, axis=2)==0] = imgt[np.sum(img, axis=2)==0]
    #STEP THREE - populate the resulting object's .uns['spatial'] so it can be used for things
    adata.uns['spatial'] = dict()
    library_id = "joint"
    adata.uns["spatial"][library_id] = dict()
    #just the high resolution image
    adata.uns["spatial"][library_id]['images'] = dict()
    adata.uns["spatial"][library_id]['images']["hires"] = img.copy()
    #the coordinates are already scaled to account for possible differences between the samples
    #set the scale factor to 1 as it's factored in on a per sample level
    adata.uns["spatial"][library_id]['scalefactors'] = dict()
    adata.uns["spatial"][library_id]['scalefactors']['tissue_hires_scalef'] = 1
    #we need the spot diameter for plotting purposes
    #sc.pl.spatial() multiplies it by the scale factor to get spot sizes
    #and our scale factor is now 1, so extract a spot diameter and multiply it by the corresponding scale factor
    #this way the dots will be the correct size in the resulting object
    scalef = list(adatas[0].uns['spatial'].values())[0]['scalefactors']['tissue_hires_scalef']
    diam = list(adatas[0].uns['spatial'].values())[0]['scalefactors']['spot_diameter_fullres']
    adata.uns["spatial"][library_id]['scalefactors']['spot_diameter_fullres'] = diam*scalef
    #and we're done!
    return adata
