import pandas as pd


def filter_damaged_taxa(df, filter_conditions):
    """Function to filter damaged taxa
    Args:
        df (panda.DataFrame): A dataframe containing metaDMG results
        filter_conditions (dict): A dictionary with filter conditions
        taxonomic_rank (str): Select the taxonomic rank to filter
    Returns:
        pandas.DataFrame: A filtered dataframe containing metaDMG results
    """

    mdmg_results = df.loc[
        (df[list(filter_conditions)] >= pd.Series(filter_conditions)).all(axis=1)
    ].copy()

    return mdmg_results


def load_mdmg_results(file_path):
    """Function to read a mdmg results file to a pandas dataframe
    Args:
        file_path (str): A file path pointing to a mdmg results file
    Returns:
        pandas.DataFrame: A pandas dataframe containing the mdmg results
    """
    mdmg_results = pd.read_csv(file_path, sep=",", index_col=None)
    mdmg_results.rename(columns={"tax_id": "reference"}, inplace=True)
    return mdmg_results


def load_fb_results(file_path):
    """Function to read a filterBAM results file to a pandas dataframe
    Args:
        file_path (str): A file path pointing to a mdmg results file
    Returns:
        pandas.DataFrame: A pandas dataframe containing the mdmg results
    """
    fb_results = pd.read_csv(file_path, sep="\t", index_col=None)
    return fb_results


def filter_fb_references(df, filter_conditions):
    fb_results = df.loc[
        (df[list(filter_conditions)] >= pd.Series(filter_conditions)).all(axis=1)
    ].copy()

    return fb_results
