#!/usr/bin/env python3
# coding: utf-8

# Human error check
import sys
if len(sys.argv) <= 3:
    print(f"Kindly check the passed arguments.\nOnly {len(sys.argv)-1} arguments given, expected at least 3\nUSAGE:\nIoNN <hg38/hg19> <model.h5> <Vars_1.vcf> <Vars_2.vcf> ...<Vars_n.vcf>")
    exit()

# Dependencies
import tensorflow as tf
import pandas as pd
import os
import allel
import myvariant
mv = myvariant.MyVariantInfo()
hg = "hg38"

# My name
software = "IntelliNGS v0.03a by Aditya Singh"

# Picking params
hg = sys.argv[1]
if "hg38" in str(hg).lower():
    assembly = "hg38"
elif "hg37" in str(hg).lower():
    hg = "hg37"
else:
    print(f"You forgot to set or wrongly set the human genome assembly version\n'{sys.argv[1]}'\nChoose either hg38 or hg37\nTry again!")
    exit()
#________________________________________________

if os.path.isfile(sys.argv[2]) and str(sys.argv[2]).endswith(".h5"):
    model = tf.keras.models.load_model(sys.argv[2])
    print("Model loaded successfully!")
    print(model.summary())
else:
    print(f"Model file incorrect!\n'{sys.argv[2]}' not found/ incorrect!")
    exit()
# Module
def IntelliNGS(file, genome_ver):
    os.system(f'bcftools norm -m -both {file} > {file[:-4]}_SPLIT.vcf')
    file = f'{file[:-4]}_SPLIT.vcf'
    params = ['QUAL',"AF","AO","DP","FAO","FDP","FDVR","FRO","FSAF","FSAR","FSRF","FSRR","FWDB","FXX","HRUN","LEN","MLLD","PB","PBP","QD","RBI","REFB","REVB","RO","SAF","SAR","SRF","SRR","SSEN","SSEP","SSSB","STB","STBP","TYPE","VARB"]
    variations = allel.read_vcf(file)
    if variations is None:
        print(f'The file {file} has no variants.\nSkipping...')
        return
    name = ''
    name = str(variations['samples'][0])
    vcf = allel.read_vcf(file, fields=["variants/QUAL","variants/AF","variants/AO","variants/DP","variants/FAO","variants/FDP","variants/FDVR","variants/FRO","variants/FSAF","variants/FSAR","variants/FSRF","variants/FSRR","variants/FWDB","variants/FXX","variants/HRUN","variants/LEN","variants/MLLD","variants/PB","variants/PBP","variants/QD","variants/RBI","variants/REFB","variants/REVB","variants/RO","variants/SAF","variants/SAR","variants/SRF","variants/SRR","variants/SSEN","variants/SSEP","variants/SSSB","variants/STB","variants/STBP","variants/TYPE","variants/VARB"])
    df = pd.DataFrame(columns=params)
    for i in params:
        df[i] = vcf[f'variants/{i}']
    zygo = []
    for z in variations['calldata/GT']:
        if 0 in z:
            zygo.append('Heterozygous')
        else:
            zygo.append('Homozygous')
    df["TYPE"].replace(to_replace="snp", value=0, inplace=True)
    df["TYPE"].replace(to_replace="mnp", value=1, inplace=True)
    df["TYPE"].replace(to_replace="ins", value=2, inplace=True)
    df["TYPE"].replace(to_replace="del", value=3, inplace=True)
    df["TYPE"].replace(to_replace="complex", value=4, inplace=True)
    results = model.predict_classes(df.values)
    results_prob = model.predict_proba(df.values)
    df['IntelliNGS'] = results
    output = pd.DataFrame()
    output['Chromosome'] = variations['variants/CHROM']
    output['Position'] = variations['variants/POS']
    output['Reference'] = variations['variants/REF']
    output['Alternate'] = variations['variants/ALT'][:,0]
    output['ID'] = variations['variants/ID']
    output['Zygosity'] = zygo
    output['HGVS'] = output['Chromosome'].apply(str) + ":g." + output['Position'].apply(str) + output['Reference'].apply(str)+'>'+output['Alternate'].apply(str)
    output['VarSome'] = output['Chromosome'].apply(str) + "-" + output['Position'].apply(str) + '-' + output['Reference'].apply(str)+'-'+output['Alternate'].apply(str)
    output[["QUAL","DP"]] = df[["QUAL","DP"]]
    output['IntelliNGS'] = df['IntelliNGS']
    output['Fail_Prob'] = results_prob[:,0]
    output['Pass_Prob'] = results_prob[:,1]
    ann = mv.getvariants(output["HGVS"], assembly=genome_ver, as_dataframe=True)
    output['IntelliNGS'] = output['IntelliNGS'].replace(0, 'FAIL')
    output['IntelliNGS'] = output['IntelliNGS'].replace(1, 'PASS')
    output.reset_index(drop=True, inplace=True)
    ann.reset_index(drop=True, inplace=True)
    output = pd.concat([output, ann], axis=1)
    pass_out = output[output['IntelliNGS'].str.contains("PASS")]
    pass_out.to_excel(f'{file[:-4]}_{name}_IntelliNGS_ONLY_PASS.xlsx', index=False, sheet_name='IntelliNGS_Predictions_All_Pass')
    output.to_excel(f'{file[:-4]}_{name}_IntelliNGS.xlsx', index=False, sheet_name='IntelliNGS_Predictions')
    return

# Parser
for i in sys.argv[3:]:
    if os.path.isfile(i):
        IntelliNGS(i, hg)
    else:
        print(f"{i} not found on system.\nCheck the address and try again\nSkipping..")
print(f"All processes finished.\nThank you for using {software}\nFeel free to contact the developer:\nAditya Singh\naditya.onco@gmail.com")