 
import cv2
from Parameters import *
import math
 
def build_filters():
    filters = []
    for lamda in lamdaList:
        for theta in thetaList:
            for psi in psiList:

                sigma = float(lamda/math.pi*math.sqrt(math.log(2)/2)*(2**bw+1)/(2**bw-1))
                kern = cv2.getGaborKernel((ksize, ksize), sigma, theta, lamda, gamma, psi, ktype=cv2.CV_32F)
                filters.append(kern)

    return filters


 
def process(img, filters):
    filteredImages = [];
    for kern in filters:
        fimg = cv2.filter2D(img, cv2.CV_8UC3, kern)
        fimg = fimg;#/mean(fimg);
        filteredImages.append(cv2.resize(fimg,(layerGDim,layerGDim),interpolation = cv2.INTER_AREA))

        # plt.imshow(cv2.resize(fimg,(layerGDim,layerGDim),interpolation = cv2.INTER_AREA  ))
        # plt.show()
        
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
    
