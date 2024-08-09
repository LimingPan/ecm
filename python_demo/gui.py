import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def matri_bar(data=None, print_plot=True):
    if data is None:
        print("no data provided, execution stops")
        return

    if not isinstance(data, pd.DataFrame):
        print("data should be in data.frame format, execution stops")
        return

    if 'workflow' not in data.attrs or data.attrs['workflow'] != "matrisomeannotatoR":
        print("graphs can only be drawn for annotated files, execution stops")
        return

    gdf = data.copy()

    # Create data frames for counts
    v1 = gdf['Annotated Matrisome Division'].value_counts().reset_index()
    v1.columns = ['Var1', 'Freq']
    v1['source'] = "Annotated Matrisome Division"

    v2 = gdf['Annotated Matrisome Category'].value_counts().reset_index()
    v2.columns = ['Var1', 'Freq']
    v2['source'] = "Annotated Matrisome Category"

    # Combine data
    d1 = pd.concat([v1, v2])

    # Define color mapping for 'Var1'
    color_map = {
        "Drosophila matrisome": "#B118DB",
        "Nematode-specific core matrisome": "#B118DB",
        "Nematode-specific matrisome-associated": "#741B47",
        "Putative Matrisome": "#B118DB",
        "Core matrisome": "#002253",
        "Matrisome-associated": "#DB3E18",
        "Apical Matrix": "#B3C5FC",
        "Cuticular Collagens": "#B3C5FC",
        "Cuticlins": "#BF9000",
        "ECM Glycoproteins": "#13349D",
        "Collagens": "#0584B7",
        "Proteoglycans": "#59D8E6",
        "ECM-affiliated Proteins": "#F4651E",
        "ECM Regulators": "#F9A287",
        "Secreted Factors": "#FFE188",
        "Non-matrisome": "#D9D9D9"
    }

    # Map colors to Var1 values
    d1['color'] = d1['Var1'].map(color_map).fillna('#D9D9D9')  # Default color if not in color_map

    # Plot for each source on a separate figure
    sources = d1['source'].unique()

    for source in sources:
        plt.figure(figsize=(12, 8))
        subset = d1[d1['source'] == source]
        p1 = sns.barplot(data=subset, x='Var1', y='Freq', hue='Var1', legend=False, palette=color_map)

        # Fix xticks and labels
        p1.set_xticks(range(len(subset['Var1'])))
        p1.set_xticklabels(subset['Var1'], rotation=90, ha='right')

        p1.set_xlabel("")
        p1.set_ylabel("Counts")
        p1.set_title(source)
        plt.tight_layout()  # Adjust layout to fit labels and titles

        if print_plot:
            plt.show()
        else:
            plt.savefig(f'{source}.png')  # Save figure if not printing

    # If print_plot is False, the function will save plots as files in the current directory.
