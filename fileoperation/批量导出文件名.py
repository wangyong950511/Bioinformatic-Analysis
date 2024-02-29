import os
from openpyxl import Workbook

folder_path = '/Users/wangyong/Documents/CASH/CASHrawdata/Cell/YHRS-2024.2.29.nosync/RG+FFAs'  # 替换成你想读取的文件夹路径

# 检查路径是否存在并且是一个目录
if os.path.exists(folder_path) and os.path.isdir(folder_path):
    # 获取文件夹中的所有文件名
    file_names = os.listdir(folder_path)

    # 创建一个新的Excel工作簿
    workbook = Workbook()
    sheet = workbook.active

    # 将文件名写入Excel表格的第一列
    for index, file_name in enumerate(file_names, start=1):
        sheet.cell(row=index, column=1, value=file_name)

    # 保存Excel文件
    workbook.save('file_names.xlsx')
else:
    print("指定的路径不存在或不是一个有效的文件夹路径。")
