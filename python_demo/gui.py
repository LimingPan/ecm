import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import plotly.graph_objects as go

from python_demo.common import data_check2


def create_color_map(unique_values, default_color="#D9D9D9"):
    """
    Create a color map from unique values based on a base color map.

    Parameters:
    - unique_values: List of unique values in the data
    - base_colors: Dictionary with predefined color mappings
    - default_color: Default color if a value is not in base_colors

    Returns:
    - color_map: Dictionary mapping unique values to colors
    """
    color_map = {}

    base_colors = {
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

    for i, value in enumerate(unique_values):
        # Generate a color if not in base_colors
        if value in base_colors:
            color_map[value] = base_colors[value]
        else:
            color_map[value] = plt.cm.tab20.colors[i % len(plt.cm.tab20.colors)]
    return color_map


def matri_bar(data=None, print_plot=True):
    data_check2(data)

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

    # Generate color map based on data
    unique_values = d1['Var1'].unique()
    color_map = create_color_map(unique_values)

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


def matri_flow(data=None, print_plot=True):
    if data is None:
        print("no data provided, execution stops")
        return

    if not isinstance(data, pd.DataFrame):
        print("data should be in data.frame format, execution stops")
        return

    if 'workflow' not in data.attrs or data.attrs['workflow'] != "matrisomeannotatoR":
        print("graphs can only be drawn for annotated files, execution stops")
        return

    # Create the alluvial data
    d2 = pd.crosstab(data['Annotated Matrisome Division'], data['Annotated Matrisome Category'])
    d2 = d2.reset_index()
    d2 = d2.melt(id_vars='Annotated Matrisome Division', var_name='Annotated Matrisome Category', value_name='Freq')
    d2 = d2[d2['Freq'] > 0]
    d2 = d2.fillna("Non-matrisome")

    # Define color mapping
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

    # Mapping the colors to categories
    color_map_rev = {v: k for k, v in color_map.items()}
    d2['color'] = d2['Annotated Matrisome Category'].map(color_map).fillna('#D9D9D9')

    # Prepare the data for the Sankey diagram
    labels = list(pd.concat([d2['Annotated Matrisome Division'], d2['Annotated Matrisome Category']]).unique())
    label_map = {label: i for i, label in enumerate(labels)}

    source_indices = d2['Annotated Matrisome Division'].map(label_map)
    target_indices = d2['Annotated Matrisome Category'].map(label_map)

    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels
        ),
        link=dict(
            source=source_indices,
            target=target_indices,
            value=d2['Freq'],
            color=d2['color']
        )
    )])

    fig.update_layout(
        title_text="Matrisome Flow",
        font_size=10
    )

    if print_plot:
        fig.show()
    else:
        fig.write_image("matri_flow.png")  # Save the figure if not printing


def matri_ring(data=None, print_plot=True):
    data_check2(data)

    # Create data frame for counts
    d2 = pd.crosstab(data['Annotated Matrisome Division'], data['Annotated Matrisome Category'])
    d2 = d2.reset_index()
    d2 = d2.melt(id_vars='Annotated Matrisome Division', var_name='Annotated Matrisome Category', value_name='Freq')
    d2 = d2[d2['Freq'] > 0]

    # Define color mapping
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

    d2['color'] = d2['Annotated Matrisome Category'].map(color_map).fillna('#D9D9D9')

    # Data preparation for plotting
    fin = d2.groupby('Annotated Matrisome Category').agg({'Freq': 'sum', 'color': 'first'}).reset_index()
    fin['Percent'] = (fin['Freq'] / fin['Freq'].sum() * 100).round(1)
    fin['lab'] = fin['Annotated Matrisome Category'] + ' (' + fin['Freq'].astype(str) + ', ' + fin['Percent'].astype(
        str) + '%)'
    fin = fin.sort_values(by='Freq')
    fin['ymax'] = fin['Freq'].cumsum()
    fin['ymin'] = fin['ymax'] - fin['Freq']

    # Plotting
    fig, ax = plt.subplots(figsize=(12, 8), subplot_kw=dict(aspect="equal"))

    wedges, _ = ax.pie(
        fin['Freq'],
        colors=fin['color'],
        wedgeprops=dict(width=0.5, edgecolor='w'),
        startangle=90,
        autopct=None
    )

    # Add legend with percentages
    ax.legend(
        wedges,
        fin['lab'],
        title="Annotated Matrisome Categories",
        loc="center left",
        bbox_to_anchor=(1, 0, 0.5, 1),
        fontsize=10
    )

    plt.title("Annotated Matrisome Categories")
    plt.tight_layout()  # Adjust layout to fit labels and legends

    # Save the figure if not printing
    if print_plot:
        plt.show()
    else:
        plt.savefig('matri_ring_with_legend.png')


def matri_star(data=None, print_plot=True):
    if data is None:
        print("no data provided, execution stops")
        return

    if not isinstance(data, pd.DataFrame):
        print("data should be in data.frame format, execution stops")
        return

    if 'workflow' not in data.attrs or data.attrs['workflow'] != "matrisomeannotatoR":
        print("graphs can only be drawn for annotated files, execution stops")
        return

    # Create data frame for counts
    d2 = pd.crosstab(data['Annotated Matrisome Division'], data['Annotated Matrisome Category'])
    d2 = d2.reset_index()
    d2 = d2.melt(id_vars='Annotated Matrisome Division', var_name='Annotated Matrisome Category', value_name='Freq')
    d2 = d2[d2['Freq'] > 0]

    # Define color mapping
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

    d2['color'] = d2['Annotated Matrisome Category'].map(color_map).fillna('#D9D9D9')

    # Data preparation for plotting
    fin = d2.groupby('Annotated Matrisome Category').agg({'Freq': 'sum', 'color': 'first'}).reset_index()
    fin['Percent'] = (fin['Freq'] / fin['Freq'].sum() * 100).round(1)
    fin['lab'] = fin['Annotated Matrisome Category'] + ' (' + fin['Freq'].astype(str) + ')'
    fin = fin.sort_values(by='Freq')

    # Check if Non-matrisome needs scaling
    max_freq = fin['Freq'].max()
    non_matrisome_freq = fin.loc[fin['Annotated Matrisome Category'] == 'Non-matrisome', 'Freq'].values
    if len(non_matrisome_freq) > 0 and non_matrisome_freq[0] > 4 * max_freq:
        print("warning: the Non-matrisome bar was scaled down to aid visualization")
        fin.loc[fin['Annotated Matrisome Category'] == 'Non-matrisome', 'Freq'] = 2 * max_freq

    # Prepare for star plot
    num_vars = len(fin)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    angles += angles[:1]  # Complete the loop
    frequencies = fin['Freq'].tolist()
    frequencies += frequencies[:1]  # Complete the loop
    colors = fin['color'].tolist()
    colors += colors[:1]  # Complete the loop
    labels = fin['lab'].tolist()
    labels += labels[:1]  # Complete the loop

    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw=dict(polar=True))

    # Draw the bars
    bars = ax.bar(angles, frequencies, color=colors, width=0.4, edgecolor='w', linewidth=1, zorder=2)

    # Draw the white circle in the middle to create the hollow effect
    centre_circle = plt.Circle((0, 0), 0.5, color='white', fc='white', linewidth=0)
    fig.gca().add_artist(centre_circle)

    # Add labels
    for i in range(num_vars):
        angle = angles[i]
        label = labels[i]
        rotation = np.degrees(angle) - 90
        ax.text(angle, frequencies[i] + 5, label, ha='center', va='center', rotation=rotation, rotation_mode='anchor', fontsize=10)

    # Draw the lines
    for i in range(num_vars):
        ax.plot([angles[i], angles[i]], [0, frequencies[i]], color='black', linewidth=0.5, linestyle='--', zorder=1)

    # Remove the frame and set limits
    ax.spines['polar'].set_visible(False)
    ax.set_yticklabels([])
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(fin['Annotated Matrisome Category'].tolist(), fontsize=10)

    # Add legend
    legend_elements = [Patch(facecolor=color_map[label], label=label) for label in fin['Annotated Matrisome Category']]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.1, 1), loc='upper left', fontsize=10)

    plt.title("Annotated Matrisome Categories")
    plt.tight_layout()

    if print_plot:
        plt.show()
    else:
        plt.savefig('matri_star_plot.png')