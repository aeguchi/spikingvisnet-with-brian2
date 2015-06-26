#!/usr/bin/env python
 
import numpy as np
import cv2
from Parameters import *
import matplotlib.pyplot as plt
import math
from numpy import mean
 
def build_filters():
    filters = []
    for lamda in lamdaList:
        for theta in thetaList:
            for psi in psiList:
                # cv2.getGaborKernel(ksize, sigma, theta, lambda, gamma, psi, ktype)
                #ksize = int(lamda)+1
                sigma = float(lamda/np.pi*math.sqrt(math.log(2)/2)*(2**bw+1)/(2**bw-1))
                kern = cv2.getGaborKernel((ksize, ksize), sigma, theta, lamda, gamma, psi, ktype=cv2.CV_32F)
                #kern /= 1.5 * kern.sum()
                #kern /= 0.2 * kern.sum()
                filters.append(kern)
    return filters


 
def process(img, filters):
    filteredImages = [];
    for kern in filters:
        fimg = cv2.filter2D(img, cv2.CV_8UC3, kern)
        fimg = fimg/mean(fimg);
        filteredImages.append(fimg)
    return filteredImages



def imageLoad(img_fn):
    import sys
    
    img_fn = 'test2.png'
    
    img = cv2.imread(img_fn)
    if img is None:
        print 'Failed to load image file:', img_fn
        sys.exit(1)
    
    filters = build_filters()
    
    res = process(img, filters)
    
    
    for r in res:
        cv2.imshow('result', r)
        cv2.waitKey(0)
        cv2.destroyAllWindows()
        
    return res;
        


# if __name__ == '__main__':
#     import sys
#     
#     print __doc__
#     try:
#         img_fn = sys.argv[1]
#     except:
#         img_fn = 'test2.png'
#     
#     img = cv2.imread(img_fn)
#     if img is None:
#         print 'Failed to load image file:', img_fn
#         sys.exit(1)
#     
#     filters = build_filters()
#     
#     res = process(img, filters)
#     
#     
#     for r in res:
#         cv2.imshow('result', r)
#         cv2.waitKey(0)
#         cv2.destroyAllWindows()
#         
    
