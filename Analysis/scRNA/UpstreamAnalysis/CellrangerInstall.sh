#!/bin/bash

# Download cellranger form https://www.10xgenomics.com/support/software/cell-ranger/downloads
#chose your ture url
# ex
curl -o cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1701529203&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=YBofZ-etPeRA3wVfDtpgadpx-MrqIQ38gZ8nbSLr1Dm5XIEuNnNpoVwuk4V9NPDI8vP6UtNJ9MdnwMe8pdJrTsiUofAqPK0dFzHH~Tyxt7foa9UeUOud0luLBehqQCVkrydA3uH6XcPhI-Dlxt97q3gYXgbJ5Rc4vWv4-5kVbfieVo3HFqD7TyPY5CRlGNonA93LTWcq3moQkCF-vtidPqAbLfCdS706pT7Tk9kU2RXvoN3vFk0q5xHpzqXk8Qv0JFZYWcy6Ni8pf86ZkXdqV3t7v5vSQlQv4WhtOnbeK6MVv~pbySPPD~W0Jet9mkaksSjs9Vc3a3Z4hM8NQh4iQA__"

# 比对下载包是否完整
md5sum cellranger-7.1.0.tar.gz

# 下载参考值
# 鼠参考
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
# 人参考
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz

# 切换到下载包所在文件夹后解压
tar -zxvf cellranger-7.1.0.tar.gz
