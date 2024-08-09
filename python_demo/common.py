import pandas as pd


def data_check1(data):
    if data is None:
        print("no data provided, execution stops")
        return

    if not isinstance(data, pd.DataFrame):
        print("data should be in data.frame format, execution stops")
        return


def data_check2(data):
    data_check1(data)

    if 'workflow' not in data.attrs or data.attrs['workflow'] != "matrisomeannotatoR":
        print("data should be annotated first, execution stops")
        return
