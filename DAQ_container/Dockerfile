FROM nvidia/cuda:11.7.1-cudnn8-runtime-ubuntu22.04

COPY DAQ-main DAQ-main

COPY requirements.txt .

RUN apt update && \
    apt install --no-install-recommends -y build-essential python3 python3-pip && \
    apt clean && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir --upgrade -r requirements.txt


# WORKDIR /DAQ-main/assign_score
# RUN make

# WORKDIR /DAQ-main/process_map/gen_trimmap
# RUN make

WORKDIR /DAQ-main
RUN mkdir -p inputs

ENTRYPOINT ["python3", "main.py"]