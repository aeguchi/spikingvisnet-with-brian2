from Parameters import *
from imageImport import *
import sys

#construct NetworkState

#Training
#load gabor filtered ImageSurface
    
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
    
