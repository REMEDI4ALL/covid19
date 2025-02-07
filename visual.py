import math
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats


def plot_compounds_interactive(compounds, dmso=None, n_components=None):
    drug_names = {
        id_: name_ for (id_, name_) in zip(compounds['batch_id'], compounds['name'])
    }

    drugs = compounds[~compounds['batch_id'].isin(['non-inf', 'DMSO'])]
    fig = px.line(
        drugs.sort_values(by="conc"), x="conc", y="distance",
        color='name',
        markers=True,
        width=1200, height=400,
        error_y='distance_q_75', error_y_minus='distance_q_75',
        category_orders={"name": drug_names.values()},
        hover_name="name", 
        hover_data={"name": False, "count_nuclei": ':.1f', "conc": True, "distance": ':.3f'},
    )

    # Non-infected cells: confidence intervals based on the Chi-squared with df=`n_components`
    if n_components is not None:
        for i, conf in enumerate([0.5, 0.75, 0.95]):
            label = dict(
                text="Non-infected", font=dict(size=15, color="black"), textposition="middle left"
            ) if i == 0 else None

            left, right = stats.chi2.interval(confidence=conf, df=n_components)
            fig.add_hrect(
                type="rect",
                y0=math.sqrt(left), y1=math.sqrt(right),
                fillcolor="royalblue", opacity=0.8*(1-conf)**0.5,
                layer="below", line_width=0, label=label
            )

    # DMSO
    if dmso is not None:
        conc = range(0, 31, 1)
        median = dmso['distance'].iloc[0]
        perc_25 = dmso['distance_q_25'].iloc[0]
        perc_75 = dmso['distance_q_75'].iloc[0]
        fig.add_hrect(
            type="rect",
            y0=median-perc_25, y1=median+perc_75,
            fillcolor="grey", opacity=0.25,
            layer="below", line_width=0,
            label=dict(
                text="DMSO", font=dict(size=15, color="black"), textposition="middle left"
            )
        )

    conc = compounds['conc'].unique()
    fig.update_xaxes(title_text="concentration", gridcolor='lightgrey', type="log", tickvals=conc)
    fig.update_yaxes(title_text="distance", gridcolor='lightgrey', type="log", tickvals=[1, 5, 10, 100])
    fig.update_layout(
        title_text="Distance to non-infected cell distribution",
        plot_bgcolor='white'
    )

    fig.show()
