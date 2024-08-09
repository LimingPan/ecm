import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri


def get_matrisome_data():
    # Enable conversion between pandas and R data frames
    pandas2ri.activate()
    # Load the R `load` function
    robjects.r['load']('../data/matrisome.list.rda')
    # Read the R object, assuming matrisome_list is the object name in R
    matrisome_list = robjects.globalenv['matrisome.list']

    human_matrisome = pandas2ri.rpy2py(matrisome_list.rx2('human'))
    mouse_matrisome = pandas2ri.rpy2py(matrisome_list.rx2('mouse'))
    c_elegans_matrisome = pandas2ri.rpy2py(matrisome_list.rx2('c.elegans'))
    zebrafish_matrisome = pandas2ri.rpy2py(matrisome_list.rx2('zebrafish'))
    drosophila_matrisome = pandas2ri.rpy2py(matrisome_list.rx2('drosophila'))

    return {'human': human_matrisome,
            'mouse': mouse_matrisome,
            'c.elegans': c_elegans_matrisome,
            'zebrafish': zebrafish_matrisome,
            'drosophila': drosophila_matrisome}


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
    matrisome_list = get_matrisome_data()

    # Select the corresponding matrisome data based on species
    if species == "human":
        k = matrisome_list['human']
    elif species == "mouse":
        k = matrisome_list['mouse']
    elif species == "c.elegans":
        k = matrisome_list['c.elegans']
        k['family'] = k['family'].apply(lambda x: "ECM-affiliated Proteins" if x == "ECM-affiliated" else x)
    elif species == "zebrafish":
        k = matrisome_list['zebrafish']
        k['category'] = k['category'].apply(lambda x: "Non-matrisome" if x == "not.available" else x)
        k['family'] = k['family'].apply(lambda x: "Non-matrisome" if x in ["", "not.available"] else x)
    elif species == "drosophila":
        k = matrisome_list['drosophila']
        k['category'] = k['category'].apply(
            lambda x: "Matrisome-associated" if x == "Homologs/Orthologs to Mammalian Matrisome-Associated Genes" else (
                "Core matrisome" if x == "Homologs/Orthologs to Mammalian Core Matrisome Genes" else "Drosophila matrisome"))
        k['family'] = k['family'].apply(lambda x: "ECM-affiliated Proteins" if x == "ECM-affiliated" else x)
    else:
        print(f"Species {species} is not recognized, execution stops")
        return

    # Merge data frames
    left_key = 'gene'
    df2 = pd.merge(k, df, left_on=left_key, right_on=n, how='right').drop_duplicates()
    df2 = df2.drop(columns=left_key)  # Drop the key column from the left data frame
    # Reorder the key column 'gene' to be the first column
    columns = [n] + [col for col in df2.columns if col != n]
    df2 = df2[columns]

    # Handle column names
    # Rename only the first 3 columns
    new_column_names = ['Annotated Gene', 'Annotated Matrisome Division', 'Annotated Matrisome Category']
    if len(df2.columns) > 3:
        new_column_names.extend(df2.columns[3:])

    # Update column names
    df2.columns = new_column_names

    # Fill missing values
    df2['Annotated Gene'] = df2['Annotated Gene'].replace("", "gene name missing in original data")
    df2.update(df2[['Annotated Matrisome Division']].fillna("Non-matrisome"))
    df2.update(df2[['Annotated Matrisome Category']].fillna("Non-matrisome"))
    df2.fillna("", inplace=True)

    # Set attributes
    df2.attrs['workflow'] = "matrisomeannotatoR"

    return df2


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
    tr = tr.loc[:, ~tr.isna().any()]

    nmtr = tr.columns.tolist()

    # Merge with the original data
    tr = pd.concat([df[[n]], tr], axis=1)

    # Ensure columns are numeric, avoiding issues with duplicate column names
    tr.iloc[:, 1:] = tr.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')

    # Aggregate data
    df2_2 = df.iloc[:, :3]
    bf = pd.merge(df2_2, tr, left_on=n, right_on=n).drop_duplicates()

    # Prepare columns for aggregation: "Annotated Matrisome Division" and columns from the 4th column onward
    division_cols = ['Annotated Matrisome Division'] + bf.columns[3:].tolist()
    category_cols = ['Annotated Matrisome Category'] + bf.columns[3:].tolist()

    # Aggregate based on "Annotated Matrisome Division"
    a = bf[division_cols].groupby('Annotated Matrisome Division').sum()

    # Aggregate based on "Annotated Matrisome Category"
    b = bf[category_cols].groupby('Annotated Matrisome Category').sum()

    # Concatenate the results
    z = pd.concat([a, b])

    z.index.name = "Matrisome Annotation"
    z.attrs['workflow'] = "matrisomeanalyzeR"

    return z
