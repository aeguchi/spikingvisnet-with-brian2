from Parameters import *
import scipy
import scipy.misc
import scipy.signal
import os
import sys
import numpy as np
import numpy.linalg;
import math
import pylab




class GaborFilter(object):

    """Gabor filter"""

    def __init__(self):
        self.filtering()

    def filtering(self):
        # file
        psi = [0, np.pi, -np.pi / 2, np.pi / 2];
        scale = lamdaList;
        orient = thetaList;
        # gamma
        # set
        paddingColor = 128;
        
        fileList_train = np.genfromtxt(os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/fileList_train.txt", dtype='str');
        numObj_train = int(fileList_train[0]);
        numTrans_train = int(fileList_train[1]);

        for index_obj in range(numObj_train):
            for index_trans in range(numTrans_train):
                print "obj: " + str(index_obj) + ", trans: " + str(index_trans);
                index_img = index_obj * numTrans_train + index_trans;
                file = os.path.split(os.path.realpath(__file__))[0] + "/images/" + imageFolder + "/train/" + fileList_train[index_obj * numTrans_train + index_trans + 2];
                # img = cv2.imread(img_fn)
                I = scipy.misc.imread(file, flatten=True)
                fsize = len(I) * len(I);
                originalImageD = np.array(([len(I), len(I)]));
                paddedImageD = 2 * originalImageD;
                imageSize = len(I);
                ss = len(scale);
                ors = len(orient);
                ps = len(psi);
                # FI = np.array((ss, ors, ps));
                # GF = np.array((ss, ors, ps));
                # minF = np.zeros((ss, ors, ps));
                # maxF = np.zeros((ss, ors, ps));

        
                # print I;
        
                print 'Filtering image: ' + file;
                print 'Orientations: ' + str(np.degrees(orient)) + ' (degrees)';
                print 'Wavelengths: ' + str(scale) + ' (pixels)';
                print 'Phases: ' + str(np.degrees(psi)) + ' (degrees)';
        
                # pylab.imshow(I);
                # pylab.show();
                
        
             
                # Create a new image of twice the dimensions
                # and copy the original image into the center.
                # This resulting image is then convolved, and
                # we only keep convolution result values for the
                # part of the image where the original image 
                # is located.
         
                # Make padded image, set to provided background color
                paddedImage = paddingColor * np.ones((paddedImageD[0], paddedImageD[1]));
         
                # Copy original image into center of padded image
                top = paddedImageD[0] / 2 - imageSize / 2 - 1;
                left = paddedImageD[1] / 2 - imageSize / 2 - 1;
                paddedImage[left:left + imageSize, top:top + imageSize] = I;
        
                # pylab.imshow(paddedImage,interpolation='none');
                # pylab.show();
        
        
                for p in range(ps):
                    for s in range(ss):
                        for o in range(ors):
                            # GF[s, o, p] = self.gabor_fn(psi[p], scale[s], orient[o]);
                            gf = self.gabor_fn(psi[p], scale[s], orient[o]);
                            gf = gf - numpy.mean(gf);
                            gf = gf / numpy.linalg.norm(gf[:]);
                
                            # Convolve padded image with Gabor filter
                            tmp = scipy.signal.convolve2d(paddedImage, gf, 'same');
                         # Copy out part of padded image convolution that corresponds to
                         # the original image
                            fi = tmp[(paddedImageD[0] - originalImageD[0]) / 2:(paddedImageD[0] - originalImageD[0]) / 2 + originalImageD[0], (paddedImageD[1] - originalImageD[1]) / 2:(paddedImageD[1] - originalImageD[1]) / 2 + originalImageD[1]]

                         # Normalize with selfconvolution of gabor filter to control for
                         # varying size of gabor filer for diffrent parameter
                         # combinations
                            fi = fi / np.linalg.norm(fi[:]);
                         
                           # Cut away negatives
                            fi[fi < 0] = 0;
                            
                            # Display minimum and maximum convolved image values
                            minF = np.min(np.min(fi));
                            maxF = np.max(np.max(fi));
                            # pylab.figure();
                            pylab.subplot(ss, ors, ((s) * ors) + o + 1);
                            pylab.imshow(fi, interpolation='none');
                pylab.show()
    

    def gabor_fn(self, psi, lmd, theta):
        # bw    = bandwidth, (1)
        # gamma = aspect ratio, (0.5)
        # psi   = phase shift, (0)
        # lambda= wave length, (>=2) (pixels)
        # theta = angle in rad, [0 pi)
        # The output is always real. For a complex Gabor filter, you may use (0 phase) + j(pi/2 phase).
        
        # if nargin<6; plot=false; end
        
        sigma = lmd / math.pi * np.sqrt(np.log(2) / 2) * (2 ** bw + 1) / (2 ** bw - 1);
        sigma_x = sigma;
        sigma_y = sigma / gamma;
        
        sz = math.floor(8 * np.max([sigma_y, sigma_x]));
        if np.mod(sz, 2) == 0:
            sz = sz + 1
         
        x, y = np.meshgrid(range(-int(math.floor(sz / 2)), int(math.floor(sz / 2) + 1)), range(int(math.floor(sz / 2)), int(-math.floor(sz / 2) - 1), -1));
        # [x y] = meshgrid(-fix(sz / 2):fix(sz / 2), fix(sz / 2):-1:fix(-sz / 2));
        # x (right +)
        # y (up +)
        
        # Rotation 
        x_theta = x * np.cos(theta) + y * np.sin(theta);
        y_theta = -x * np.sin(theta) + y * np.cos(theta);
         
        gb = np.multiply(np.exp(-.5 * (np.power(x_theta, 2) / sigma_x ** 2 + np.power(y_theta, 2) / sigma_y ** 2)), np.cos(2 * np.pi / lmd * x_theta + psi));
        
        # pylab.imshow(gb,interpolation='none');
        # pylab.show();
        
        return gb
        # if (plot)
        #    figure(); imshow(gb / 2 + 0.5);  # Rescale from [-1,1] to [0,1]
        #    figure(); meshc(x, y, gb);
        # end
        
        # http://www.mathworks.com/matlabcentral/fileexchange/23253
