import pandas as pd

from python_demo.gui import matri_bar
from python_demo.util import matriannotate, matrianalyze


def run():
    data = pd.read_csv('../data/mass-spec.csv')

    ann = matriannotate(data=data, gene_column="Gene Symbol", species="human")
    tbl = matrianalyze(ann)

    matri_bar(ann)


if __name__ == '__main__':
    run()
