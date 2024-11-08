import os
import openpyxl
import re

# 自然排序的辅助函数
def natural_sort_key(s):
    # 匹配字符串中的数字，用于自然排序
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]

def export_filenames_to_excel(folder_path, output_excel):
    # 获取文件夹中所有文件的文件名，并按自然顺序排序
    filenames = os.listdir(folder_path)
    filenames.sort(key=natural_sort_key)  # 使用自然排序

    # 创建新的 Excel 工作簿和工作表
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = "Filenames"

    # 写入表头
    sheet["A1"] = "Filename"

    # 写入每个文件的文件名，去除后缀
    for index, filename in enumerate(filenames, start=2):
        # 去掉文件的后缀名
        base_filename = os.path.splitext(filename)[0]
        sheet[f"A{index}"] = base_filename

    # 保存 Excel 文件
    workbook.save(output_excel)
    print(f"文件名（无后缀，按自然顺序）已导出到 {output_excel}")

# 示例：指定文件夹路径和导出文件路径
folder_path = '/Users/wangyong/Desktop/大样本/Tiff_Raw'  # 替换为你的文件夹路径
output_excel = '/Users/wangyong/Desktop/文件名列表.xlsx'  # 替换为你希望导出的文件路径

export_filenames_to_excel(folder_path, output_excel)
