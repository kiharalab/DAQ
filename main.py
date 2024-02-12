import os
from ops.argparser import argparser
from ops.os_operation import mkdir
import time
def compile_online(code_path):
    exe_path = os.path.join(code_path,"DAQscore_colab")
    if params.get('server') != 1:
        if os.path.exists(exe_path):
            os.remove(exe_path)
        root_path = os.getcwd()
        os.chdir(code_path)
        os.system("make")
        os.chdir(root_path)
        if not os.path.exists(exe_path):
            print("Assign score compilation failed! Please make contact with dkihara@purdue.edu!")
        assert os.path.exists(exe_path)
    return exe_path

def init_save_path(origin_map_path):
    save_path = os.path.join(os.getcwd(), 'Predict_Result')
    mkdir(save_path)
    map_name = os.path.split(origin_map_path)[1].replace(".mrc", "")
    map_name = map_name.replace("(","").replace(")","")
    save_path = os.path.join(save_path, map_name)
    mkdir(save_path)
    return save_path
import subprocess
import time
def get_gpu_memory_usage():
    result = subprocess.run(['nvidia-smi', '--query-gpu=memory.used', '--format=csv,nounits,noheader'], stdout=subprocess.PIPE)
    output = result.stdout.decode('utf-8').strip()
    memory_usage = [int(x) for x in output.split('\n')]
    return memory_usage
if __name__ == "__main__":
    params = argparser()
    if params['mode']==0:
        choose = params['gpu']
        if choose is not None:
            os.environ["CUDA_VISIBLE_DEVICES"] = choose
        cur_map_path = os.path.abspath(params['F'])
        #process the map path if it's ending with .gz
        if ".gz"==cur_map_path[-3:]:
            from ops.os_operation import unzip_gz
            cur_map_path = unzip_gz(cur_map_path)
        model_path = os.path.abspath(params['M'])
        pdb_path = os.path.abspath(params['P'])
        if params['output'] is None:
            save_path = init_save_path(cur_map_path)
        else:
            save_path=params['output']
            mkdir(save_path)
        #give pdb path to save computing voxels
        from data_processing.generate_trimmap import generate_trimmap
        save_path,trimmap_path = generate_trimmap(save_path,cur_map_path,pdb_path,params)
        print("Finished processing model input!")
        from predict.predict_trimmap import predict_trimmap
        output_path = os.path.join(save_path, "prediction.txt")
        #if not os.path.exists(output_path):
       # gpu = get_gpu_memory_usage()
       # free_gpu =12288 - max(gpu)
        #while free_gpu<=6000:
        #    print("waiting for available gpu")
        #    time.sleep(60)
        #    gpu = get_gpu_memory_usage()
         #   free_gpu =12288 - max(gpu)
        output_path = predict_trimmap(trimmap_path, save_path, model_path, params)
        print("Our predictions are saved in %s, please have a check!"%output_path)
        #further call daq score to output the final score
        daq_code_path = os.path.join(os.getcwd(),"assign_score")
        exe_path = compile_online(daq_code_path)
        
        raw_score_save_path = os.path.join(save_path,"daq_raw_score.pdb")
        os.system("chmod 777 "+exe_path)
        map_name = os.path.split(cur_map_path)[1].replace(".mrc", "")
        map_name = map_name.replace("(","").replace(")","")
        new_map_path = os.path.join(save_path,map_name+"_new.mrc")
        print("exe_path")
        print(exe_path+" -i "+new_map_path+" -p "+output_path+" -Q "+str(pdb_path)+" >"+raw_score_save_path)
        os.system(exe_path+" -i "+new_map_path+" -p "+output_path+" -Q "+str(pdb_path)+" >"+raw_score_save_path)

        #smooth the score to give the final output
        from ops.process_raw_score import read_pdb_info,get_resscore,save_pdb_with_score,read_chain_set,get_resscore_atom
        window_size = params['window']
        window_size_atom = 9
        score_save_path = os.path.join(save_path,"daq_score_w"+str(window_size)+".pdb")
        score_save_path_atom = os.path.join(save_path,"daq_score_atom_w"+str(window_size_atom)+".pdb")
        chain_list = read_chain_set(pdb_path)
        print("total different chains:",chain_list)
        for chain_name in chain_list:
            score_chain_save_path = os.path.join(save_path,"daq_score_w"+str(window_size)+"_"+str(chain_name)+".pdb")
            score_dict = get_resscore(raw_score_save_path,window_size,chain_name)
            residue_dict = read_pdb_info(pdb_path,chain_name)
            save_pdb_with_score(score_dict, residue_dict,score_chain_save_path)

            score_chain_save_path_atom = os.path.join(save_path,"daq_score_atom_w"+str(window_size_atom)+"_"+str(chain_name)+".pdb")
            score_dict_atom = get_resscore_atom(raw_score_save_path,window_size_atom,chain_name)
            residue_dict_atom = read_pdb_info(pdb_path,chain_name)
            save_pdb_with_score(score_dict_atom, residue_dict_atom, score_chain_save_path_atom)
        #concat all chain visualization together
        with open(score_save_path,'w') as wfile:
            for chain_name in chain_list:
                score_chain_save_path = os.path.join(save_path,"daq_score_w"+str(window_size)+"_"+str(chain_name)+".pdb")
                with open(score_chain_save_path,'r') as rfile:
                    line = rfile.readline()
                    while line:
                        wfile.write(line)
                        line = rfile.readline()

        with open(score_save_path_atom,'w') as wfile:
            for chain_name in chain_list:
                score_chain_save_path_atom = os.path.join(save_path,"daq_score_atom_w"+str(window_size_atom)+"_"+str(chain_name)+".pdb")
                with open(score_chain_save_path_atom,'r') as rfile:
                        line = rfile.readline()
                        while line:
                            wfile.write(line)
                            line = rfile.readline()
        print("Please check result here: %s"%score_save_path)
        print("Please check result here: %s"%score_save_path_atom)

        print("Please open it in pymol and visualize it by putting the following command to Pymol:")
        print("-"*100)
        print("spectrum b, red_white_blue,  all, -1,1")
        print("-"*100)
    elif params['mode']==1:
        #input dir and then process
        choose = params['gpu']
        if choose is not None:
            os.environ["CUDA_VISIBLE_DEVICES"] = choose
        cur_map_path = os.path.abspath(params['F'])
        #process the map path if it's ending with .gz
        if ".gz"==cur_map_path[-3:]:
            from ops.os_operation import unzip_gz
            cur_map_path = unzip_gz(cur_map_path)
        model_path = os.path.abspath(params['M'])

        pdb_dir = os.path.abspath(params['P'])
        if params['output'] is None:
            save_path = init_save_path(cur_map_path)
        else:
            save_path=params['output']
            mkdir(save_path)
        from ops.pdb_utils import  concat_pdbs
        pdb_path = os.path.join(save_path,"concat.pdb")
        concat_pdbs(pdb_dir,pdb_path)

        from data_processing.generate_trimmap import generate_trimmap
        save_path,trimmap_path = generate_trimmap(save_path,cur_map_path,pdb_path,params)
        print("Finished processing model input!")
        from predict.predict_trimmap import predict_trimmap
        output_path = os.path.join(save_path, "prediction.txt")
        #if not os.path.exists(output_path):
        output_path = predict_trimmap(trimmap_path, save_path, model_path, params)
        print("Our predictions are saved in %s, please have a check!"%output_path)
        #further call daq score to output the final score
        daq_code_path = os.path.join(os.getcwd(),"assign_score")
        exe_path = compile_online(daq_code_path)
        map_name = os.path.split(cur_map_path)[1].replace(".mrc", "")
        map_name = map_name.replace("(","").replace(")","")
        new_map_path = os.path.join(save_path,map_name+"_new.mrc")

        #smooth the score to give the final output
        from ops.process_raw_score import read_pdb_info,get_resscore,save_pdb_with_score,read_chain_set
        window_size = params['window']
        listfiles = [x for x in os.listdir(pdb_dir) if ".pdb" in x]
        root_save_path = save_path
        for pdb_item in listfiles:
            pdb_path = os.path.join(pdb_dir,pdb_item)
            save_path = os.path.join(root_save_path,pdb_item)
            mkdir(save_path)
            raw_score_save_path = os.path.join(save_path,"daq_raw_score.pdb")
            os.system(exe_path+" -i "+new_map_path+" -p "+output_path+" -Q "+str(pdb_path)+" >"+raw_score_save_path)
            score_save_path = os.path.join(save_path,"daq_score_w"+str(window_size)+".pdb")
            chain_list = read_chain_set(pdb_path)
            print("%s total different chains :"%pdb_item,chain_list)
            for chain_name in chain_list:
                score_chain_save_path = os.path.join(save_path,"daq_score_w"+str(window_size)+"_"+str(chain_name)+".pdb")
                score_dict = get_resscore(raw_score_save_path,window_size,chain_name)
                residue_dict = read_pdb_info(pdb_path,chain_name)
                save_pdb_with_score(score_dict, residue_dict,score_chain_save_path)
            #concat all chain visualization together
            with open(score_save_path,'w') as wfile:
                for chain_name in chain_list:
                    score_chain_save_path = os.path.join(save_path,"daq_score_w"+str(window_size)+"_"+str(chain_name)+".pdb")
                    with open(score_chain_save_path,'r') as rfile:
                        line = rfile.readline()
                        while line:
                            wfile.write(line)
                            line = rfile.readline()
            print("Please check result here: %s"%score_save_path)
        print("Please open it in pymol and visualize it by putting the following command to Pymol:")
        print("-"*100)
        print("spectrum b, red_white_blue,  all, -1,1")
        print("-"*100)




