# Intelli-NGS
***Intelli-NGS:Intelligent NGS, a deep neural network-based artificial intelligence to delineate good and bad variant calls from IonTorrent sequencer data***
Preprint: DOI https://doi.org/10.1101/2019.12.17.879403

If you've a TensorFlow compatible NVIDIA chip configured with CUDA and CUDNN libraries then kindly use:

pip3 install requirements_gpu.txt --user

Otherwise, to install with CPU support only, use:

pip3 install requirements.txt --user

Once installed, kindly use:

python3 Intelli-NGS.py <hg38/hg19> <Model.h5> <Vars_1.vcf> <Vars_2.vcf> ...<Vars_n.vcf>

Thank you for using the software, kindly cite it!

For any suggestions/ queries, feel free to contact me at:

# Aditya Singh

**aditya.onco@gmail.com**
