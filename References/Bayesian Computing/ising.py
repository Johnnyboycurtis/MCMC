# just ignore this - it's for using Python 3 division
from __future__ import division
# module used for timing
import time

# package used for plotting
import matplotlib as mpl
# subpackage used for a MATLAB-like interface to plotting
import matplotlib.pyplot as plt
# package used for optimizing the execution of numerical expressions
import numexpr as ne
# package for numerical Python - this is arguably the single most important 3rd party Python package
import numpy as np
# subpackage of scipy which is used for scientific Python - this is another important and huge package
import scipy.stats
# just ignore this - we're using it for easy image loading
import skimage

# some setup of matplotlib
mpl.rcParams['image.interpolation'] = 'none'
mpl.rcParams['figure.figsize'] = (16.0, 8.0)


# function for initializing the image by using nearest neighbour on the pixels to [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
def x_init(img):
    # round values to [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    img = np.round(img * 10) / 10
    # replace 0.0 values by 0.1 values
    img[img < 0.1] = 0.1
    # replace 1.0 values by 0.9 values
    img[img > 0.9] = 0.9
    return img


# naive loopy version of the function
# NOTE: this only performs 1 iteration
def loop_ising(img, beta=100):
    # create a nearest neighbour version of the image padded with 1 pixel with value 0.0
    x = np.zeros(np.array(img.shape) + 2)
    x[1:-1, 1:-1] = x_init(img)

    # this is used in the g function to obtain the density of the given Gaussian distribution
    z_dist = scipy.stats.norm(loc=0 , scale=img.std())

    def g(x_s, x_nes, z_s):
        # the following should be straight-forward with x_nes being [x_ne_north, x_ne_west, x_ne_east, x_ne_west]
        return z_dist.pdf(z_s - x_s) * np.exp(beta * np.sum(x_s == x_nes))

    # loop over the rows of the image (not including the padding)
    for i in range(1, img.shape[0] + 1):
        # loop over the pixels of the row (not including the padding)
        for j in range(1, img.shape[1] + 1):
            # draw a proposal value from our proposal distribution
            x_sp = np.random.choice([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
            # x_nes = [x_ne_north, x_ne_west, x_ne_east, x_ne_west]
            x_nes = np.array([x[i, j-1], x[i+1, j], x[i, j+1], x[i-1, j]])

            alpha = np.min([1, g(x_sp, x_nes, img[i-1, j-1]) / g(x[i, j], x_nes, img[i-1, j-1])])
            
            if np.random.rand() < alpha:
                x[i, j] = x_sp

    return x[1:-1, 1:-1]


def vector_ising(img, beta=100, iterations=100):
    # create a nearest neighbour version of the image padded with 1 pixel with value 0.0
    A = np.zeros(np.array(img.shape) + 2)
    A[1:-1, 1:-1] = x_init(img)

    # x is the image (centre part of A - NOT a copy)
    x = A[1:-1, 1:-1]
    # x_ne are the neighbour images, i.e., the image shifted 1 up, 1 left, 1 right, and 1 down, respectively (again NOT copies)
    x_ne = (A[:-2, 1:-1], A[1:-1, :-2], A[1:-1, 2:], A[2:, 1:-1])
    # z is the original image
    z = img

    # the density function of a normal random variable is
    #     1 / (sqrt(2 * pi) * standard_deviation) * exp(-1 / (2 * variance) * (x - mean)^2)
    # for which we define the constants pdf_c1 and pdf_c2 so that the expression becomes
    #     pdf_c1 * exp(pdf_c2 * (x - mean)^2)
    pdf_c1 = 1 / (np.sqrt(2 * np.pi) * img.std())
    pdf_c2 = -1 / (2 * img.var())

    # ignore the "pdf_c1=pdf_c1, pdf_c2=pdf_c2, beta=beta" part - due to Python and numexpr mechanics it is required here
    def g(x, x_ne, z, pdf_c1=pdf_c1, pdf_c2=pdf_c2, beta=beta):
        # "x == x_ne[i]" creates a logical array for which "+" would mean logical "OR" - by using "1 * " it is converted to a numerical array
        beta_sum = 1 * (x == x_ne[0]) + (x == x_ne[1]) + (x == x_ne[2]) + (x == x_ne[3])
        # ne.evaluate simply optimizes the expression which has been rewritten to the given form from
        #     pdf_c1 * exp(pdf_c2 * (x - z)^2) * exp(beta * sum_of_x_equals_neighbours)
        return ne.evaluate('pdf_c1 * exp(pdf_c2 * (x - z)**2 + beta * beta_sum)')

    # this might seem a bit magical - it simply creates a way of masking the "chess-board structure" with mask1 and the inverse with mask2
    mask = np.zeros((img.shape[0], img.shape[1] + 1), dtype=np.bool)
    mask[::2] = 1
    mask1 = mask[:, :-1]
    mask2 = np.invert(mask1)

    for i in range(iterations):
        # first update one colour of the chess-board then the other
        for mask in (mask1, mask2):
            # draw proposal values from our proposal distribution
            x_p = np.random.choice([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], size=x.shape)
            # get rid of alpha and use H instead as this does not change the functionality
            h = g(x_p, x_ne, z) / g(x, x_ne, z)
            # draw random values for comparison with alpha to decide to accept or reject proposal values
            p = np.random.random(size=h.shape)

            # the values whish should be updated are those of the chess-board which we are working on (mask) and which we would accept (p < h)
            mask = (p < h) * mask
            x[mask] = x_p[mask]

    return x


# load image with values from 0 to 255
img = plt.imread('thorningS.ppm')
# rescale image to have values from 0.0 to 1.0
img = skimage.img_as_float(img)
plt.figure()
# plot the original image
plt.imshow(img, cmap='gray')
plt.figure()
# plot the initialized image (nearest neighbour applied)
plt.imshow(x_init(img), cmap='gray')
plt.show()


start_time = time.time()
img_processed = vector_ising(img, beta=4, iterations=100)
end_time = time.time()

print('Processing took {} seconds'.format(end_time - start_time))
plt.imshow(img_processed, cmap='gray')
plt.show()
