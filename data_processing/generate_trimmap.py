import os
from ops.os_operation import mkdir
import shutil


def generate_trimmap(save_path,origin_map_path,input_pdb_path,params):

    map_name = os.path.split(origin_map_path)[1].replace(".mrc", "")
    save_path = os.path.join(save_path, map_name)
    mkdir(save_path)
    cur_map_path = os.path.join(save_path,map_name+".mrc")
    if not os.path.exists(cur_map_path):
        shutil.copy(origin_map_path,cur_map_path)

    new_map_path = cur_map_path[:-4]+"_new.mrc"
    from process_map.Reform_Map_Voxel import Reform_Map_Voxel, Reform_Map_Voxel_Final
    raise_exception_flag = False
    try:
        Reform_Map_Voxel(cur_map_path, new_map_path)
    except:
        raise_exception_flag = True
        Reform_Map_Voxel_Final(cur_map_path, new_map_path)
    finally:
        print("-" * 100)
        print("resizing finished!")
        if raise_exception_flag:
            print("if it's a map that has different sizes in different axis, we suggest to use eman2 e2proc3d.py to change it's to symmetrical")
            print("Example command:")
            print("e2proc3d.py input.mrc output.mrc --clip max_dim_size")
        print("-" * 100)
    trimmap_path = os.path.join(save_path, map_name + ".trimmap")
    factor = params['stride']
    voxel_size = params['voxel_size']
    #contour_level = params['contour']

    half_voxel_size = int((voxel_size - 1) / 2)
    code_path = os.path.join(os.getcwd(), 'process_map')
    code_path = os.path.join(code_path, 'gen_trimmap')
    run_code_path = os.path.join(code_path, 'TrimMapAtom')
    root_path = os.getcwd()
    os.chdir(code_path)
    if os.path.exists(run_code_path):
        os.remove(run_code_path)
    os.system("make clean")
    os.system("rm *.o")
    os.system("make")  # compile the code
    os.chdir(root_path)
    commandline = run_code_path + ' -i ' + new_map_path + ' -p '+str(input_pdb_path)+\
                      ' -v ' + str(half_voxel_size) + \
                      ' -s ' + str(factor) + ' -L 0.005 >' + trimmap_path
    #if not os.path.exists(trimmap_path) or os.path.getsize(trimmap_path)<=10000:
    os.system(commandline)
    return save_path,trimmap_path
