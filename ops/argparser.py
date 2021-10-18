

import parser
import argparse

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-F',type=str, help='Map file path')#File path for decoy dir
    parser.add_argument('-M', type=str,default="best_model/qa_model/Multimodel.pth",  help='Phase1 model path')
    parser.add_argument('--mode',type=int,required=True,help='Running Mode')
    parser.add_argument("--contour",type=float,default=0,help="Map contour level")
    parser.add_argument("--stride",type=int,default=1,help="Stride size for scanning maps (default:1)")
    parser.add_argument("--voxel_size",type=int,default=11,help="Phase1 input voxel size")
    parser.add_argument("--gpu",type=str,default=None,help="specify the gpu we will use")
    parser.add_argument('--batch_size', type=int, default=256, help='batch size for inference')
    parser.add_argument('--cardinality', default=32, type=int, help='ResNeXt cardinality')
    args = parser.parse_args()
    params = vars(args)
    return params