
import torch
from torch import nn
import torch.nn.functional as F
import os
import numpy as np
from ops.os_operation import  mkdir
from models.resnet import resnet18 as resnet18_multi
import time

def predict_trimmap(trimmap_path, save_path,model_path, params):
    model = resnet18_multi(sample_size=params['voxel_size1'])
    model = model.cuda()
    # model = nn.DataParallel(model, device_ids=None)

    state_dict = torch.load(model_path)
    model = nn.DataParallel(model, device_ids=None)
    model.load_state_dict(state_dict['state_dict'])
    model.eval()

    voxel_size1 = params['voxel_size1']
    batch_size = params['batch_size']
    output_line = ""
    end_time = time.time()
    with open(trimmap_path, 'r') as file:
        line = file.readline()
        origin = np.zeros(3)
        width = np.zeros(3)
        input_use = np.zeros([batch_size, 1, voxel_size1, voxel_size1, voxel_size1])

        Position_List = []
        count_example = 0
        while line:
            #
            line = line.strip()
            if line[:7] == "#orgXYZ":
                split_result = line.split()
                for k in range(3):
                    origin[k] = float(split_result[k + 1])
            if line[:11] == "#New orgXYZ":
                split_result = line.split()
                for k in range(3):
                    origin[k] = float(split_result[k + 2])
            if line[:10] == "#Grid size":
                split_result = line.split()
                for k in range(3):
                    width[k] = float(split_result[k * 2 + 3])
            if line[:6] == "AtomID":
                current_position = np.zeros(3)
                relative_position = np.zeros(3)
                split_result0 = line.split(":")
                voxel = split_result0[-1]
                density_list = voxel.split(",")
                assert len(density_list) == voxel_size1 ** 3 + 1
                try:
                    for index in range(len(density_list) - 1):
                        x = index % (voxel_size1)
                        y = int(((index - x) % (voxel_size1 ** 2)) / voxel_size1)
                        z = int((index - x - y * voxel_size1) / (voxel_size1 ** 2))
                        input_use[count_example, 0, x, y, z] = float(density_list[index])
                except:
                    print("!!!residue information wrong:")
                    print(density_list)
                    # continue
                    continue
                split_result1 = line.split(",")
                split_result2 = int(split_result1[2].split(":")[1])

                split_result1 = split_result1[:-len(density_list)]
                current_position[0] = int(split_result1[3].split(":")[1])
                current_position[1] = int(split_result1[4])
                current_position[2] = int(split_result1[5])
                tmp_coord_key = ""
                for kk in range(3):
                    tmp_coord_key += str(int(current_position[kk])) + ","

                Position_List.append(current_position)
                count_example += 1
            if count_example == batch_size:
                input_data = torch.from_numpy(input_use).float()
                # print(input_data.size())
                input_data = input_data.cuda()
                with torch.no_grad():
                    pred1, pred2, pred3 = model(input_data)
                    pred1 = F.softmax(pred1, dim=1)
                    pred2 = F.softmax(pred2, dim=1)
                    pred3 = F.softmax(pred3, dim=1)

                pred1 = pred1.cpu().detach().numpy()
                pred2 = pred2.cpu().detach().numpy()
                pred3 = pred3.cpu().detach().numpy()
                for i in range(len(Position_List)):
                    tmp_key = Position_List[i]
                    wline = ""
                    for kk in range(3):
                        wline += str(int(tmp_key[kk])) + ","
                    tmp_coord_key = wline

                    pred_prob1 = pred1[i]
                    pred_prob2 = pred2[i]
                    pred_prob3 = pred3[i]
                    for kk in range(len(pred_prob1)):
                        wline += str(pred_prob1[kk]) + ","
                    for kk in range(len(pred_prob2)):
                        wline += str(pred_prob2[kk]) + ","
                    for kk in range(len(pred_prob3)):
                        wline += str(pred_prob3[kk]) + ","
                    output_line +=wline+"\n"

                # print(Prediction_Dict)
                # exit()
                print("batch processing time %f"%(time.time()-end_time))
                end_time = time.time()
                count_example = 0
                Position_List = []
                input_use = np.zeros([batch_size, 1, voxel_size1, voxel_size1, voxel_size1])

            line = file.readline()
        if count_example != 0:
            input_use = input_use[:count_example]
            input_data = torch.from_numpy(input_use).float()
            # print(input_data.size())
            input_data = input_data.cuda()
            with torch.no_grad():

                pred1, pred2, pred3 = model(input_data)
                pred1 = F.softmax(pred1, dim=1)
                pred2 = F.softmax(pred2, dim=1)
                pred3 = F.softmax(pred3, dim=1)

            pred1 = pred1.cpu().detach().numpy()
            pred2 = pred2.cpu().detach().numpy()
            pred3 = pred3.cpu().detach().numpy()
            for i in range(len(Position_List)):
                tmp_key = Position_List[i]
                wline = ""
                for kk in range(3):
                    wline += str(int(tmp_key[kk])) + ","
                tmp_coord_key = wline

                pred_prob1 = pred1[i]
                pred_prob2 = pred2[i]
                pred_prob3 = pred3[i]


                for kk in range(len(pred_prob1)):
                    wline += str(pred_prob1[kk]) + ","
                for kk in range(len(pred_prob2)):
                    wline += str(pred_prob2[kk]) + ","
                for kk in range(len(pred_prob3)):
                    wline += str(pred_prob3[kk]) + ","
                output_line+=wline+"\n"
    prediction_path = os.path.join(save_path, "prediction.txt")
    with open(prediction_path, 'w') as file:
        file.write(output_line)