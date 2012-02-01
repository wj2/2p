"""
Image registration
"""

from numpy.fft import rfft2,irfft2
import numpy as np
from matplotlib import pyplot as plt
from scipy.ndimage.interpolation import shift
import scipy.ndimage
#import time

def balanced_mod(x,d):
    x = np.asarray(x)
    d = np.asarray(d)
    shift = (d-1)//2
    return (x+shift)%d - shift

def get_avg_img(path):
    """load data and calculate the mean image"""
    try:
        X = np.load(path)
    except:
        X = path
    M = np.array(np.mean(X,axis=0))
    return M

def test_balanced_mod():
    assert balanced_mod(10,5) == 0
    assert balanced_mod(3,6) == 3
    assert balanced_mod(4,6) == -2
"""    
def register_img(img,template):
    return the amount that image is shifted relative to template.
    you can correct the shift with     
    new_img = np.roll(np.roll(img,-shift[0],0),-shift[1],1)
    corr = irfft2(rfft2(img)*rfft2(template).conj())
    return balanced_mod(np.unravel_index(corr.argmax(),corr.shape),corr.shape)
"""
def register_imgs(imgs,template):
    "save some time by only taking fft of template once"
    rfft2_template_conj = rfft2(template).conj()
    shifts = []
    for img in imgs:
        corr = irfft2(rfft2(img)*rfft2_template_conj)
        shifts.append(balanced_mod(np.unravel_index(corr.argmax(),corr.shape),corr.shape))
    return shifts

def correct_jitter(imgs):
    shifts = register_imgs(imgs,get_avg_img(imgs))
    for i in range(len(imgs)):
        imgs[i] = shift(imgs[i],-shifts[i])
    return imgs,get_avg_img(imgs)
"""
def show_stack(stack):
    for items in stack:
        plt.imshow(items)
        time.sleep(.1)
    return 0
"""
"""
to shift the images by their shift coords
from scipy import ndimage
scipy.ndimage.interpolation.shift(img,-shift)
"""


"""
trans_avg_img = rfft2(create_average_image(imgs)).conj())
for img in imgs:
for i in range(4):

    corr = []
    corr.append(irfft2(rfft2(img)*trans_avg_img))
    np.unravel_index(corr.argmax(),corr.shape)

        
"""
