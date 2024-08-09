import pandas as pd

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

data = pd.read_csv('../data/mass-spec.csv')

# 启用 pandas 和 R 之间的数据帧转换
pandas2ri.activate()
# 加载R的 `load` 函数
robjects.r['load']('../data/matrisome.list.rda')
# 读取R对象，假设matrisome_list 是R中的对象名
matrisome_list = robjects.globalenv['matrisome.list']

human_matrisome = pandas2ri.rpy2py(matrisome_list.rx2('human'))
mouse_matrisome = pandas2ri.rpy2py(matrisome_list.rx2('mouse'))
c_elegans_matrisome = pandas2ri.rpy2py(matrisome_list.rx2('c.elegans'))
zebrafish_matrisome = pandas2ri.rpy2py(matrisome_list.rx2('zebrafish'))
drosophila_matrisome = pandas2ri.rpy2py(matrisome_list.rx2('drosophila'))