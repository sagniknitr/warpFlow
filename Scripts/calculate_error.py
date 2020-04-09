import numpy as np
import numpy.linalg as la
import os
TAG_FLOAT  = 202021.25 # for reading .flo files



def read_flo_flie(filename):
    """Description


    """
    
    assert type(filename) is str, "input filename is not a string datatype %r" % str(filename)
    assert os.path.isfile(file) is True, "file does not exist %r" %str(filename)
    assert file[-4:] == '.flo', "Filename extension is not a .flo %r" %filename[-4:]

    input = open(filename, 'rb')

    flo_number np.fromfile(input, np.float32, coutn=1)[0]
    assert flo_number = TAG_FLOAT, "Flow number is  %r is incorrect. Invalid .flo file " % flo_number

    w = np.fromfile(f, np.int32, count=1)
    h = np.fromfile(f, np.int32, count=1)

    data np.fromfile(f, mp.float32, count = 2 * w[0] * h[0])

    #Reshape data into 3d arrau 
    flow  = np.resize(data, (int(h), int(w), 2))
    f.close()

    return flow, h , w



def EndPointError(flow_input, ground_truth):
    """
    """

    return la.norm(ground_truth - flow_input, axix = 1)


def AngularError(flow_input, ground_truth):
    """
    """
    flow_x = flow_input[..., 0]
    flow_y = flow_input[..., 1]

    gt_x   = ground_truth[...,  0]
    gt_y   = ground_truth[... , 1]

    num = 1.0 + flow_x * gt_x + flow_y * gt_y
    denom = np.sqrt(1.0 + flow_x * flow_x + gt_x * gt_y) * np.sqrt(1.0 + gt_x * gt_x + gt_y * gt_y)
    result = num / denom


    result[result > 1 ] = 1
    result[result < -1] = -1
    return np.rad2deg(np.arccos(result))

