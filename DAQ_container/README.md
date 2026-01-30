# DAQ Docker

This directory contains the Docker configuration for running DAQ-score with GPU support.

## Prerequisites

- [Docker](https://docs.docker.com/get-docker/)
- [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) (for GPU acceleration)

## Usage

### 1. Build the image

From the project root:
```bash
docker build -t daq-score DAQ_container/
```

### 2. Run the container

Mount your data to the `/DAQ-main/inputs` directory inside the container.

```bash
docker run --rm --gpus all \
  -v $(pwd):/DAQ-main/inputs \
  daq-score \
  -F inputs/example/2566_3J6B_9.mrc \
  -P inputs/example/3J6B_9.pdb \
  --mode 0 \
  --output inputs/test_output
```

### Arguments

- `-F`: Path to the Map file (.mrc).
- `-P`: Path to the PDB file or directory.
- `--mode`: `0` for single PDB, `1` for directory of PDBs.
- `--output`: Directory to save results.

