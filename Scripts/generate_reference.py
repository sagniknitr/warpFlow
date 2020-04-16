import file_parser as fp
import cv2
import sys
import os
import util as ut
import numpy as np




def evaluate_parameter(config_item):
    prev_img        = cv2.imread(config_item["files"]["prevImg"])
    curr_img        = cv2.imread(config_item["files"]["currImg"])
    flow_method     = config_item["parameter"]["flow_method"]
    estimate_base   = config_item["files"]["estimatepath"]  + "/"
    
    if os.path.exists(estimate_base) == False:
       os.makedirs(estimate_base)

    if os.path.exists(config_item["files"]["estflow"]):
        return
    #  compute optical flow
    if  flow_method.find("dual") >= 0:
        dual_proc = cv2.DualTVL1OpticalFlow_create(config_item["parameter"]["tau"],
                                                   config_item["parameter"]["lambda"],
                                                   config_item["parameter"]["theta"],
                                                   config_item["parameter"]["nscales"],
                                                   config_item["parameter"]["warps"])
        est_flow = np.zeros(shape=(prev_img.shape[0], prev_img.shape[1],2), dtype=np.float32)
        dual_proc.calc(cv2.cvtColor(prev_img, cv2.COLOR_BGR2GRAY), cv2.cvtColor(curr_img, cv2.COLOR_BGR2GRAY), est_flow)
    #
    elif flow_method.find("farneback") >= 0:
        est_flow = cv2.calcOpticalFlowFarneback(cv2.cvtColor(prev_img, cv2.COLOR_BGR2GRAY),
                                                cv2.cvtColor(curr_img, cv2.COLOR_BGR2GRAY),
                                                None, 0.5, 3, 15, 3, 5, 1.2, 0)
    elif flow_method.find("PLK") >= 0:
        prev_pts = list()
        for r in range(prev_img.shape[0]):
            for c in range(prev_img.shape[1]):
                prev_pts.append((c,r))
        prev_pts = np.array(prev_pts, dtype=np.float32)
        parameters = get_plk_params()
        curr_pts, st, err = cv2.calcOpticalFlowPyrLK(cv2.cvtColor(prev_img, cv2.COLOR_BGR2GRAY),
                                                cv2.cvtColor(curr_img, cv2.COLOR_BGR2GRAY),
                                               prev_pts, None,
                                               parameters.winSize, parameters.maxLevel=3, parameters.criteria, parameters.20,
                                               parameters.Epsilon)
        est_flow = np.zeros(shape=(prev_img.shape[0], prev_img.shape[1],2), dtype=np.float32)
        n = 0
        flow_pts = curr_pts - prev_pts
        for r in range(prev_img.shape[0]):
            for c in range(prev_img.shape[1]):
                est_flow[r, c, :] = flow_pts[n,:]
                n = n + 1
    #here alternative optical flow methods can be applied
    #
    else:
        raise ValueError("flow method has not been implemented")

    ut.writeFlowFile(config_item["files"]["estflow"], est_flow)
    ut.drawFlowField(config_item["files"]["estflow"][:-3] + "png", est_flow)
    print("Done -> ", config_item["files"]["estflow"])

def create_dual_parameter(parameter_list):
    parameter_list.append({"flow_method": "dual",
                            "tau" : 0.25,
                            "lambda" : 0.08,
                            "theta" : 0.37,
                            "nscales" : 6,
                            "warps" : 5})
    return parameter_list

def main():
    if len(sys.argv) < 2:
        print("Please provide the first argument that is the root path of the dataset")
        return
    base = sys.argv[1]

    flow_method_list = list()
    for n in range(2, len(sys.argv)):
        flow_method_list.append(sys.argv[n])


    parameter_list = list()
    for flow_method in flow_method_list:
        if flow_method == "dual":
            parameter_list = create_dual_parameter(parameter_list)
        else:
            parameter_list.append({"flow_method": flow_method})

    for parameter in parameter_list:
        basepath_dict = {"basepath": base,
                         "images" : base + "/images/",
                         "groundtruth": basepath + "/gt_flow/",
                         "estimate": basepath + "/estimate/" + parameter["flow_method"],
                         "masks": basepath + "/masks/",
                         }

        filenames = fp.create_filename_list(basepath_dict)
        config_list = ut.create_config(parameter, filenames)

        for config_item in config_list:
            run_parameter(config_item)

