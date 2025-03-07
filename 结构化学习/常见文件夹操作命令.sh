

# 解压tar
mkdir -p Rawdata
tar -xvf *.tar -C Rawdata/



# 解压gz
cd Rawdata
gunzip *.gz
