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
        pdb_path = os.path.abspath(params['P'])
        save_path = os.path.join(os.getcwd(), 'Predict_Result')
        mkdir(save_path)
        #give pdb path to save computing voxels
        from data_processing.generate_trimmap import generate_trimmap
        save_path,trimmap_path = generate_trimmap(save_path,cur_map_path,pdb_path,params)
        print("Finished processing model input!")
        from predict.predict_trimmap import predict_trimmap
        output_path = os.path.join(save_path, "prediction.txt")
        #if not os.path.exists(output_path):
        output_path = predict_trimmap(trimmap_path, save_path, model_path, params)
        print("Our predictions are saved in %s, please have a check!"%output_path)
        #further call daq score to output the final score
        
        daq_code_path = os.path.join(os.getcwd(),"predict")
        daq_code_path = os.path.join(daq_code_path,"DAQscore_colab")
        raw_score_save_path = os.path.join(save_path,"dqa_raw_score.pdb")
        os.system("chmod 777 "+daq_code_path)
        map_name = os.path.split(cur_map_path)[1].replace(".mrc", "")
        new_map_path = os.path.join(save_path,map_name+"_new.mrc")
        os.system(daq_code_path+" -i "+new_map_path+" -p "+output_path+" -Q "+str(pdb_path)+" >"+raw_score_save_path)
        from ops.process_raw_score import read_pdb_info,get_resscore,save_pdb_with_score
        window_size = params['window']
        score_save_path = os.path.join(save_path,"dqa_score_w"+str(window_size)+".pdb")
        
        score_dict = get_resscore(raw_score_save_path,window_size)
        residue_dict = read_pdb_info(pdb_path)
        save_pdb_with_score(score_dict, residue_dict,score_save_path)
   