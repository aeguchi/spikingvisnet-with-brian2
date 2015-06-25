#!/usr/bin/env python
 
import numpy as np
import cv2
 
def build_filters():
    filters = []
    ksize = 3#31 the size of the Gabor kernel. If ksize = (a, b), we then have a Gabor kernel of size a x b pixels. As with many other convolution kernels, ksize is preferably odd and the kernel is a square (just for the sake of uniformity).
#    for theta in np.arange(0, np.pi, np.pi / 16):
#        #cv2.getGaborKernel(ksize, sigma, theta, lambda, gamma, psi, ktype)
#        kern = cv2.getGaborKernel((ksize, ksize), 4.0, theta, 10.0, 0.5, 0, ktype=cv2.CV_32F)
#        kern /= 1.5 * kern.sum()
#        filters.append(kern)
#    return filters

    sigma = 1.5 #spatial bandwidth:  the standard deviation of the Gaussian function used in the Gabor filter.
    theta = np.pi/2#0 #orientation: the orientation of the normal to the parallel stripes of the Gabor function.
    lamda = 2.0 #wavelength: the wavelength of the sinusoidal factor in the above equation.
    gamma = 0.5 #aspect ratio:  the spatial aspect ratio.
    psi = 0 #phase shift:  phase offset.
    #(ksize, ksize), 4.0, theta, 10.0, 0.5, 0
    
    #cv2.getGaborKernel(ksize, sigma, theta, lambda, gamma, psi, ktype)
    kern = cv2.getGaborKernel((ksize, ksize), sigma, theta, lamda, gamma, psi, ktype=cv2.CV_32F)
    #kern /= 1.5 * kern.sum()
    kern /= 0.5 * kern.sum()
    filters.append(kern)
    return filters


 
def process(img, filters):
    accum = np.zeros_like(img)
    for kern in filters:
        fimg = cv2.filter2D(img, cv2.CV_8UC3, kern)
        np.maximum(accum, fimg, accum)
    return accum
 
if __name__ == '__main__':
    import sys
    
    print __doc__
    try:
        img_fn = sys.argv[1]
    except:
        img_fn = 'test.png'
    
    img = cv2.imread(img_fn)
    if img is None:
        print 'Failed to load image file:', img_fn
        sys.exit(1)
    
    filters = build_filters()
    
    res1 = process(img, filters)
    cv2.imshow('result', res1)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
