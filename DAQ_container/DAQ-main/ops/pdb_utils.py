import os

def concat_pdbs(input_dir,output_file):
    listfiles = [x for x in os.listdir(input_dir) if ".pdb" in x]
    with open(output_file,'w') as wfile:
        for item in listfiles:
            cur_file = os.path.join(input_dir,item)
            with open(cur_file,'r') as rfile:
                for line in rfile:
                    wfile.write(line)
    return output_file

