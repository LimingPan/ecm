import pandas as pd


def matrianalyze(data=None):
    if data is None:
        print("no data provided, execution stops")
        return

    if not isinstance(data, pd.DataFrame):
        print("data should be in data.frame format, execution stops")
        return

    if len(data.columns) < 1:
        print("data should be annotated first, execution stops")
        return

    if 'workflow' not in data.attrs or data.attrs['workflow'] != "matrisomeannotatoR":
        print("data should be annotated first, execution stops")
        return

    n = "Annotated Gene"
    df = data.copy()

    # Remove non-numeric characters and convert to numeric
    tr = df.apply(lambda x: pd.to_numeric(x.astype(str).str.replace(r'[^0-9.-]', '', regex=True), errors='coerce'))

    # Remove columns with all NaN values
    tr = tr.dropna(axis=1, how='all')

    nmtr = tr.columns.tolist()

    print(f'nmtr: {nmtr}')

    # Merge with the original data
    tr = pd.concat([df[[n]], tr], axis=1)
    tr[n] = pd.to_numeric(tr[n], errors='coerce')

    # Aggregate data
    df2_2 = df.iloc[:, :3]
    bf = pd.merge(df2_2, tr, left_on=n, right_on=n).drop_duplicates()

    a = bf.groupby('Annotated Matrisome Division').sum()
    b = bf.groupby('Annotated Matrisome Category').sum()
    z = pd.concat([a, b])

    z.index.name = "Matrisome Annotation"
    z.attrs['workflow'] = "matrisomeanalyzeR"

    return z


def matriannotate(data=None, gene_column=None, species=None):
    if data is None:
        print("no data provided, execution stops")
        return

    if not isinstance(data, pd.DataFrame):
        print("data should be in data.frame format, execution stops")
        return

    if gene_column is None:
        print("a column indicating gene IDs must be provided, execution stops")
        return

    if species is None:
        print("no species provided, execution stops")
        return

    df = data.copy()
    n = gene_column

    # Define the matrisome lists for different species
    if species == "human":
        k = matrisome_list["human"]
    elif species == "mouse":
        k = matrisome_list["mouse"]
    elif species == "c.elegans":
        k = matrisome_list["c.elegans"]
        k['family'] = k['family'].replace("ECM-affiliated", "ECM-affiliated Proteins")
    elif species == "zebrafish":
        k = matrisome_list["zebrafish"]
        k['category'] = k['category'].replace("not.available", "Non-matrisome")
        k['family'] = k['family'].replace({"": "Non-matrisome", "not.available": "Non-matrisome"})
    elif species == "drosophila":
        k = matrisome_list["drosophila"]
        k['category'] = k['category'].replace({
            "Homologs/Orthologs to Mammalian Matrisome-Associated Genes": "Matrisome-associated",
            "Homologs/Orthologs to Mammalian Core Matrisome Genes": "Core matrisome"
        }).fillna("Drosophila matrisome")
        k['family'] = k['family'].replace("ECM-affiliated", "ECM-affiliated Proteins")

    # Merge dataframes on gene column
    df2 = pd.merge(k, df, left_on="gene", right_on=n, how="right").drop_duplicates()

    # Rename columns
    df2.rename(columns={
        'gene': 'Annotated Gene',
        'division': 'Annotated Matrisome Division',
        'category': 'Annotated Matrisome Category'
    }, inplace=True)

    # Fill missing values
    df2['Annotated Gene'] = df2['Annotated Gene'].replace("", "gene name missing in original data")
    df2['Annotated Matrisome Division'].fillna("Non-matrisome", inplace=True)
    df2['Annotated Matrisome Category'].fillna("Non-matrisome", inplace=True)
    df2.fillna("", inplace=True)

    # Set workflow attribute
    df2.attrs['workflow'] = "matrisomeannotatoR"

    return df2