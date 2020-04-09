import numpy as np
from PIL import Image
import time
import cv2 as cv


prev_frame =  ''
next_fram = ''

im1 = np.array(Image.open(prev_frame))
im1 = im1.astype(float) / 256.0


im2 = np.array(Image.open(curr_frame))
im2 = im2.astype(float) / 256.0


np.savetext('', u, delimiter='\n')
np.savetext('', v, delimiter= '\n')

flow = np.concatenate((u[..., None], v[..., None]), axix = 2)



#plot using OpenCV
hsv = np.zeros(im1.shape, dtype = np.uint8)
hsv[:,:,0] = 255
hsv[:,:,1] = 255
mag, ang = cv.cartToPolar(flow[..., 0], flow[..., 1])
hsv[..., 1] = ang * 180/ np.pi /2
hsv[..., 2] = cv.normalize(mag, None, 0, 255, cv.NORM_MINMAX)

rgb = cv.cvtColor(hsv, cv.COLOR_HSV2BGR)