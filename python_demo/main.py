import pandas as pd

from python_demo.gui import matri_bar, matri_flow, matri_ring
from python_demo.util import matriannotate, matrianalyze


def run():
    data = pd.read_csv('../data/mass-spec.csv')

    # 'human'    'mouse'    'c.elegans'    'zebrafish'    'drosophila'
    ann = matriannotate(data=data, gene_column="Gene Symbol", species="human")
    tbl = matrianalyze(ann)
    # matri_bar(ann)
    # matri_flow(ann)
    matri_ring(ann)

if __name__ == '__main__':
    run()
