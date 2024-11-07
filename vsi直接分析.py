import os
import bioformats
import javabridge
import tifffile
import sys  # 导入 sys 模块以使用 sys.exit()


# 启动 Java 虚拟机
javabridge.start_vm(class_path=bioformats.JARS)

# 定义读取 VSI 文件并导出第一个通道的函数
def extract_channel_1(vsi_file, output_folder):
    with bioformats.ImageReader(vsi_file) as reader:
        # 获取图像的通道数
        image_count = reader.rdr.getSizeC()  # 获取通道数量

        # 打印通道数量
        print(f"Number of channels in {vsi_file}: {image_count}")

        # 只读取通道1（索引为0）
        if image_count > 0:
            image = reader.read(c=0, rescale=False)  # 读取第一个通道

            # 检查读取的图像数据是否正常
            if image is None:
                print(f"Warning: Channel 1 could not be read for {vsi_file}.")
                return

            # 获取文件名并去掉扩展名
            base_filename = os.path.splitext(os.path.basename(vsi_file))[0]

            # 保存为 TIFF 文件，命名为原文件名去掉 .vsi 后缀
            channel_filename = os.path.join(output_folder, f"{base_filename}.tiff")
            tifffile.imwrite(channel_filename, image)
            print(f"Saved channel 1 to {channel_filename}")

# 批量处理文件夹中的所有 VSI 文件，导出第一个通道
def extract_channel_1_from_folder(vsi_folder, output_folder):
    # 创建输出文件夹（如果不存在）
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 遍历文件夹中的所有 VSI 文件
    for file_name in os.listdir(vsi_folder):
        if file_name.lower().endswith('.vsi'):
            vsi_file = os.path.join(vsi_folder, file_name)
            extract_channel_1(vsi_file, output_folder)



# 示例：从指定文件夹导出第一个通道的图像
vsi_folder = '/Users/wangyong/Desktop/大样本/CASH'  # 替换为你的 VSI 文件夹路径
output_folder = '/Users/wangyong/Desktop/大样本/CASH_tiff'  # 替换为你的目标输出文件夹路径

extract_channel_1_from_folder(vsi_folder, output_folder)



