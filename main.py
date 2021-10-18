import os
from ops.argparser import argparser
from ops.os_operation import mkdir
import time

if __name__ == "__main__":
    params = argparser()
    if params['mode']==0:
        choose = params['gpu']
        if choose is not None:
            os.environ["CUDA_VISIBLE_DEVICES"] = choose
        cur_map_path = os.path.abspath(params['F'])
        model_path = os.path.abspath(params['M'])
        save_path = os.path.join(os.getcwd(), 'Predict_Result')
        mkdir(save_path)
        from data_processing.generate_trimmap import generate_trimmap
        save_path,trimmap_path = generate_trimmap(save_path,cur_map_path,params)
        print("Finished processing model input!")
        from predict.predict_trimmap import predict_trimmap
        output_path = predict_trimmap(trimmap_path, save_path, model_path, params)
        print("Our predictions are saved in %s, please have a check!"%output_path)

   