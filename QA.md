1 TrimMapAtom file missing
```
mrc.c: In function ‘readmrc’:
mrc.c:193:11: error: ‘for’ loop initial declarations are only allowed in C99 mode
for(int i = 0; iNumVoxels; ++i)
…
sh: /scratch/g/gregp/daq/DAQ/process_map/gen_trimmap/TrimMapAtom: No such file or directory
```
A: That error typically raised when you run main.py outside of DAQ directory. The change of root directory will disable the compliation of our volume reading script.

2 CUDA error
```
UserWarning: 
    There is an imbalance between your GPUs. You may want to exclude GPU 1 which
    has less than 75% of the memory or cores of GPU 0. You can do so by setting
    the device_ids argument to DataParallel, or by setting the CUDA_VISIBLE_DEVICES
    environment variable.
  warnings.warn(imbalance_warn.format(device_ids[min_pos], device_ids[max_pos]))
terminate called after throwing an instance of 'std::runtime_error'
  what():  NCCL Error 1: unhandled cuda error
Aborted (core dumped)
```
A: That is typically because of available GPU id is not correct. Please use export CUDA_VISIBLE_DEVICES=xx to set the GPU with cuda support.

3 What if If the entire model is fitted poorly? 
```
In that case, simply all the residues in the model will have a negative score, indicating that each residue in the model may be incorrectly assigned.
```
