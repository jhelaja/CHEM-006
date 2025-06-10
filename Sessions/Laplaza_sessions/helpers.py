# Helper functions! For the purpose of this tutorial, ignore these!

# Import necessary libraries
import base64  # For encoding images for Plotly hover
from io import BytesIO  # For handling image data in memory

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
import umap.umap_ as umap
from dash import Input, Output, dcc, html, no_update
from PIL import Image, ImageDraw
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Draw import MolToImage, rdMolDraw2D  # For generating images
from sklearn.cluster import DBSCAN, AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler


def smiles_to_base64_image(
    smiles_string: str, width: int = 300, height: int = 300
) -> str:
    """
    Converts a SMILES string to a base64 encoded PNG image data URI.

    Args:
        smiles_string: The SMILES string of the molecule.
        width: The width of the output image in pixels.
        height: The height of the output image in pixels.

    Returns:
        A string representing the base64 encoded PNG image,
        formatted as a data URI (e.g., "data:image/png;base64,...").
        Returns a placeholder image URI if the SMILES string is invalid
        or image generation fails.

    Usage with Plotly:
        The string returned by this function should be used as the `src`
        attribute of an HTML `<img>` tag within Plotly's `hovertemplate`.
        For example:
        hovertemplate="<img src='%{customdata[0]}' alt='Molecule'>"
        where customdata[0] would be the output of this function.
    """
    if not smiles_string:
        # print("Warning: Empty SMILES string provided.")
        # Return a placeholder for empty SMILES
        img = Image.new("RGB", (width, height), color=(230, 230, 230))
        draw = ImageDraw.Draw(img)
        draw.text((10, 10), "No SMILES", fill=(100, 100, 100))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{img_str}"

    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        # print(f"Warning: Could not parse SMILES: {smiles_string}")
        img = Image.new("RGB", (width, height), color=(220, 220, 220))  # Light grey
        draw = ImageDraw.Draw(img)
        draw.text((10, 10), "Invalid SMILES", fill=(0, 0, 0))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{img_str}"

    try:
        # Generate 2D coordinates if they don't exist or to ensure a standard depiction
        if not mol.GetNumConformers() or Chem.rdDepictor.Compute2DCoords(mol) == -1:
            # If Compute2DCoords fails or no conformers, RDKit might still draw it,
            # or it might raise an error later. This attempts to ensure coords exist.
            # For very problematic SMILES, this might not be enough, and DrawMolecule could fail.
            Chem.rdDepictor.Compute2DCoords(
                mol
            )  # Attempt again or rely on DrawMolecule's internal handling
    except Exception as e_coords:
        # print(f"Warning: Could not generate 2D coordinates for {smiles_string}: {e_coords}")
        # Fallback to a placeholder if coordinate generation fails catastrophically
        img = Image.new("RGB", (width, height), color=(220, 220, 220))
        draw = ImageDraw.Draw(img)
        draw.text((10, 10), "Coord Error", fill=(0, 0, 0))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{img_str}"

    # Use MolDraw2DCairo for direct PNG output
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    # Optional: Configure drawing options
    # drawer.drawOptions().addAtomIndices = True
    # drawer.drawOptions().kekulize = True

    try:
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        png_data = drawer.GetDrawingText()  # This returns PNG bytes for MolDraw2DCairo
    except Exception as e_draw:
        # print(f"Error drawing molecule {smiles_string}: {e_draw}")
        img = Image.new("RGB", (width, height), color=(210, 210, 210))
        draw = ImageDraw.Draw(img)
        draw.text((10, 10), "Draw Error", fill=(0, 0, 0))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{img_str}"

    base64_encoded_image = base64.b64encode(png_data).decode("utf-8")
    return f"data:image/png;base64,{base64_encoded_image}"


def plot_scatter_plotly_2d(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    color_col: str = None,
    symbol_col: str = None,
    size_col: str = None,  # For data-driven marker size
    text_col_data: pd.Series = None,  # Accepts df.index or a df column for text/hover (renamed from text_col in notebook for clarity)
    title: str = "2D Scatter Plot",
    # Parameters from notebook usage
    labels_dict: dict = None,
    marker_size: int = None,  # Fixed marker size
    text_position: str = "top center",
    smiles_for_hover_image: bool = False,  # If True, text_col_data is assumed to be SMILES (renamed from smiles_col_for_image)
    show_text_on_plot: bool = False,  # Default to NOT showing text if text_col_data is provided
    add_zerolines: bool = True,
):
    """
    Generates a 2D scatter plot.
    Adapted to work with parameters from dim_red.ipynb.
    """
    px_args = {
        "data_frame": df,
        "x": x_col,
        "y": y_col,
        "color": color_col,
        "symbol": symbol_col,
        "title": title,
    }

    text_for_display_on_plot = None
    if text_col_data is not None and show_text_on_plot:
        text_for_display_on_plot = text_col_data
    px_args["text"] = text_for_display_on_plot

    if size_col:
        px_args["size"] = size_col

    # If marker_size is given, it overrides size_col for px.scatter marker size control
    # We will apply fixed marker_size via update_traces later if size_col is not used by px.scatter
    # or if marker_size is meant to be a fixed size regardless of data.
    # If both size_col and marker_size are present, Plotly Express 'size' takes precedence for data-driven sizing.
    # Fixed marker_size will be applied if 'size' is not in px_args.

    hovertemplate_str = None
    custom_data_for_hover = None

    if smiles_for_hover_image and text_col_data is not None:
        smiles_series = text_col_data
        # Ensure smiles_series is a list or array for iteration
        smiles_list = (
            smiles_series.tolist()
            if isinstance(smiles_series, (pd.Series, pd.Index))
            else list(smiles_series)
        )
        image_data_uris = [
            smiles_to_base64_image(s, width=150, height=150) for s in smiles_list
        ]

        # Add temporary columns for hover data to a copy of the DataFrame
        df_for_plot = df.copy()
        # Ensure text_col_data (SMILES) is correctly aligned if it's an index or series
        df_for_plot["_smiles_text_"] = smiles_list
        df_for_plot["_image_uris_"] = image_data_uris

        px_args["data_frame"] = df_for_plot
        px_args[
            "hover_name"
        ] = "_smiles_text_"  # This will be the primary text on hover
        custom_data_for_hover = ["_image_uris_"]
        px_args["custom_data"] = custom_data_for_hover
        hovertemplate_str = "<b>%{hovertext}</b><br><img src='%{customdata[0]}' alt='Molecule'><extra></extra>"
    elif text_col_data is not None:
        px_args[
            "hover_name"
        ] = text_col_data  # Use the provided text_col_data for hover_name
        hovertemplate_str = "%{hovertext}<extra></extra>"  # Show only hover_name
    else:
        px_args["hover_name"] = None  # No specific hover_name, default Plotly hover

    fig = px.scatter(**px_args)

    if marker_size is not None:
        fig.update_traces(marker=dict(size=marker_size))

    if text_for_display_on_plot is not None:
        fig.update_traces(mode="markers+text", textposition=text_position)
    else:
        fig.update_traces(mode="markers")  # Ensure no text on plot by default

    if hovertemplate_str:
        fig.update_traces(hovertemplate=hovertemplate_str)

    if labels_dict:
        layout_updates = {}
        if x_col in labels_dict:
            layout_updates["xaxis_title_text"] = labels_dict[x_col]
        if y_col in labels_dict:
            layout_updates["yaxis_title_text"] = labels_dict[y_col]
        if color_col and color_col in labels_dict:
            layout_updates["legend_title_text"] = labels_dict[color_col]
        elif color_col and color_col not in labels_dict:
            layout_updates["legend_title_text"] = color_col

        if layout_updates:
            fig.update_layout(**layout_updates)

    if not add_zerolines:
        fig.update_layout(xaxis_zeroline=False, yaxis_zeroline=False)

    return fig


def plot_scatter_plotly_3d(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    z_col: str,
    color_col: str = None,
    symbol_col: str = None,
    size_col: str = None,  # Data-driven marker size
    text_col_data: pd.Series = None,  # Accepts df.index or a df column for text/hover
    title: str = "3D Scatter Plot",
    # Parameters from notebook usage
    scene_labels: dict = None,  # For 3D axis and legend titles
    marker_size: int = None,  # Fixed marker size
    text_position: str = "top center",  # For consistency, though less effective in 3D
    smiles_for_hover_image: bool = False,  # If True, text_col_data is assumed to be SMILES
    show_text_on_plot: bool = False,  # Default to NOT showing text if text_col_data is provided
):
    """
    Generates a 3D scatter plot.
    Adapted to work with parameters from dim_red.ipynb.
    """
    px_args = {
        "data_frame": df,
        "x": x_col,
        "y": y_col,
        "z": z_col,
        "color": color_col,
        "symbol": symbol_col,
        "title": title,
    }

    text_for_display_on_plot = None
    if text_col_data is not None and show_text_on_plot:
        text_for_display_on_plot = text_col_data
    px_args["text"] = text_for_display_on_plot

    if size_col:
        px_args["size"] = size_col

    hovertemplate_str = None
    custom_data_for_hover = None

    if smiles_for_hover_image and text_col_data is not None:
        smiles_series = text_col_data
        smiles_list = (
            smiles_series.tolist()
            if isinstance(smiles_series, (pd.Series, pd.Index))
            else list(smiles_series)
        )
        image_data_uris = [
            smiles_to_base64_image(s, width=150, height=150) for s in smiles_list
        ]

        df_for_plot = df.copy()
        df_for_plot["_smiles_text_"] = smiles_list
        df_for_plot["_image_uris_"] = image_data_uris

        px_args["data_frame"] = df_for_plot
        px_args["hover_name"] = "_smiles_text_"
        custom_data_for_hover = ["_image_uris_"]
        px_args["custom_data"] = custom_data_for_hover
        hovertemplate_str = "<b>%{hovertext}</b><br><img src='%{customdata[0]}' alt='Molecule'><extra></extra>"
    elif text_col_data is not None:
        px_args["hover_name"] = text_col_data
        hovertemplate_str = "%{hovertext}<extra></extra>"
    else:
        px_args["hover_name"] = None

    fig = px.scatter_3d(**px_args)

    if marker_size is not None:
        fig.update_traces(marker=dict(size=marker_size))

    if text_for_display_on_plot is not None:
        fig.update_traces(mode="markers+text")
    else:
        fig.update_traces(mode="markers")  # Ensure no text on plot by default

    if hovertemplate_str:
        fig.update_traces(hovertemplate=hovertemplate_str)

    if scene_labels:
        scene_updates = {}
        if x_col in scene_labels:
            scene_updates["xaxis_title_text"] = scene_labels[x_col]
        if y_col in scene_labels:
            scene_updates["yaxis_title_text"] = scene_labels[y_col]
        if z_col in scene_labels:
            scene_updates["zaxis_title_text"] = scene_labels[z_col]

        layout_updates = {}
        if color_col and color_col in scene_labels:
            layout_updates["legend_title_text"] = scene_labels[color_col]
        elif color_col and color_col not in scene_labels:
            layout_updates["legend_title_text"] = color_col

        if scene_updates:
            fig.update_layout(scene=scene_updates)
        if layout_updates:
            fig.update_layout(**layout_updates)

    return fig
