from Parameters import *
from imageImport import *
import sys
import pylab, glob

#construct NetworkState

#Training
trainingImages = sorted(glob.iglob("./images/training/*.png"))

img_fns = []
for img_fn in trainingImages:
    img_fns.append(img_fn)

if len(img_fns)!=nStim*nTrans:
    print 'Error: the number of images files does not match',len(img_fns);
    sys.exit(1)

#trainingImages = glob.iglob("./images/training/*.png")
for i_stim in range(1, nStim):
        for i_trans in range(1, nTrans):
            i=(i_stim-1)*nTrans+i_trans-1;
#for img_fn in trainingImages:
            img_fn = img_fns[i]
                
            print img_fn
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
