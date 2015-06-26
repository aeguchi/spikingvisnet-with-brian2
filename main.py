from Parameters import *
from imageImport import *
import sys
import pylab, glob

#construct NetworkState

#Training
trainingImages = glob.iglob("./images/training/*.png")
for img_fn in trainingImages:
    #load gabor filtered ImageSurface
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
    
    #To-DO: trainNetworkWith
    
#Testing
