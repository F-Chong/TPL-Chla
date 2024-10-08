# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 15:13:38 2023

@author: Administrator
"""
import numpy as np
import pandas as pd
from osgeo import gdal
import os
import sys
import time

FilePath = r"J:\979个湖泊处理\Landsat8" 
file_names = os.listdir(FilePath)
SavePath = r"J:\979个湖泊处理\Landsat8-Chla反演"  

def image_open(img):
    data = gdal.Open(img)
    if data is None:
        print("图像无法读取")
    return data

start_time = time.time() 
k = 0
n = 0
for m in file_names:
    k = k + 1
    father_file_path = FilePath + "/" + m

    new_file_path = SavePath + "/" + m + "_Chla_result"
    if os.path.exists(new_file_path): 
        continue
    else:
        os.mkdir(new_file_path) 

    try:
        son_file_names = os.listdir(father_file_path)
    except NotADirectoryError:
        print("该目录下没有文件夹")
        break
    else:
        if son_file_names == []:
            print("文件夹 " + m + " 为空")
            print("--------------------------------------------------------------------------------------")
            continue
        else:
            print("现在载入文件夹 " + m)
            j = 0
            for i in son_file_names: 
                QZ = os.path.splitext(i)[0]  
                HZ = os.path.splitext(i)[1] 
                if (HZ == ".tif"):
                    j = j + 1
                    image = father_file_path + "/" + i

                    data = image_open(image)
                    Coastal = data.GetRasterBand(1).ReadAsArray().astype(np.float32)
                    Blue = data.GetRasterBand(2).ReadAsArray().astype(np.float32)
                    Green = data.GetRasterBand(3).ReadAsArray().astype(np.float32)
                    Red = data.GetRasterBand(4).ReadAsArray().astype(np.float32)
                    Nir = data.GetRasterBand(5).ReadAsArray().astype(np.float32)
                    Swir = data.GetRasterBand(6).ReadAsArray().astype(np.float32)
                    Swir1 = data.GetRasterBand(7).ReadAsArray().astype(np.float32)

                    Rrs_data = pd.DataFrame({'Coastal':Coastal.reshape(-1),
                                             'Blue':Blue.reshape(-1),
                                             'Green':Green.reshape(-1),
                                             'Red':Red.reshape(-1),
                                             'Nir':Nir.reshape(-1),
                                             'Swir':Swir.reshape(-1),
                                             'Swir1':Swir1.reshape(-1)})
                    Rrs_data_cpu = Rrs_data.multiply(0.0000275).add(-0.2)

                    B1 = Rrs_data_cpu.loc[:, 'Coastal'].values.reshape(-1)
                    B2 = Rrs_data_cpu.loc[:, 'Blue'].values.reshape(-1)
                    B3 = Rrs_data_cpu.loc[:, 'Green'].values.reshape(-1)
                    B4 = Rrs_data_cpu.loc[:, 'Red'].values.reshape(-1)
                    B5 = Rrs_data_cpu.loc[:, 'Nir'].values.reshape(-1)
                    B6 = Rrs_data_cpu.loc[:, 'Swir'].values.reshape(-1)
                    B7 = Rrs_data_cpu.loc[:, 'Swir1'].values.reshape(-1)

                    cloud_mask = B3 > 0.08
                    B5[cloud_mask] = 1  

                    MNDWI = (B3 - B6) / (B3 + B6)  
                    water = MNDWI > 0
                    water_mask = ~water

                    Mask = water_mask + cloud_mask
                    B1[Mask] = 1
                    B2[Mask] = 2
                    B3[Mask] = 1
                    B4[Mask] = 0
                    B5[Mask] = 1
                    B6[Mask] = 1
                    B7[Mask] = 1

                    FAI = pd.DataFrame(
                        {'Rrc(NIR)-RrcL(NIR)': B5 - (B4 + (B6 - B4) * (859 - 645) / (1240 - 645))})  
                    FAI = FAI.values.reshape(-1)

                    FAI_mask = FAI > -0.003

                    GI = pd.DataFrame(
                        {'GI': 0.299*B4+0.587*B3+0.114*B2})  
                    GI = GI.values.reshape(-1)
                    GI_mask = GI < 0.01

                    MASK = Mask + FAI_mask + GI_mask

                    Chla = pd.DataFrame({"Ln(Chla)": (B3 - B2) / (B3 + B2)})

                    Chla_result = Chla.multiply(3.6166).add(-0.3282)

                    Chla_Result_last = np.exp(Chla_result)
                    
                    Chla_jiangwei = Chla_Result_last.values.reshape(-1)
                    maskycz =  Chla_jiangwei > 32.0
                    MASK = MASK + maskycz

                    Chla_Result_last[MASK] = 0

                    Chla_Result_last = Chla_Result_last.values

                    band1 = Chla_Result_last.reshape(Blue.shape)

                    A_output = gdal.GetDriverByName("GTiff")
                    output_0 = A_output.Create(new_file_path + "/" + QZ + "_Chla_result.tif", band1.shape[1],
                                               band1.shape[0], bands=1, eType=gdal.GDT_Float32)
                    output_0.SetProjection(data.GetProjection())
                    output_0.SetGeoTransform(data.GetGeoTransform())

                    band1 = output_0.GetRasterBand(1).WriteArray(band1)
                    output_0 = None
                    print(image + " 已完成，这是 " + m + " 文件夹中的第 " + str(j) + " 个文件")

            n = n + j

    print("文件夹 " + m + " 目录下所有结果已生成,这是第 " + str(k) + " 个文件夹")

print("--------------------------------------------------------------------------------------")
print("所有文件夹已处理完毕！")
end_time = time.time() 
print("文件夹总数为：" + str(k))
print("处理文件个数为：" + str(n))
print("处理时间:%d" % (end_time - start_time))  
sys.exit()







