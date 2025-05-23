from typing import List, Dict, Any
import pandas as pd

def get_metrics(
    df: pd.DataFrame,
    metrics: List[str] = None,
    group_by: str = None,
    group_by_values: List[str] = None,
    time_column: str = "timestamp",
    value_column: str = "value",
) -> Dict[str, Any]:
    """
    Calculate various metrics for a time series dataset.

    Args:
        df (pd.DataFrame): Input DataFrame containing time series data
        metrics (List[str], optional): List of metrics to calculate. If None, calculates all metrics.
            Available metrics: ['mean', 'std', 'min', 'max', 'range', 'cv', 'skewness', 'kurtosis']
        group_by (str, optional): Column name to group by
        group_by_values (List[str], optional): List of values to filter in the group_by column
        time_column (str): Name of the time column
        value_column (str): Name of the value column

    Returns:
        Dict[str, Any]: Dictionary containing the calculated metrics
    """
    if metrics is None:
        metrics = ['mean', 'std', 'min', 'max', 'range', 'cv', 'skewness', 'kurtosis']

    # Convert time column to datetime if it's not already
    if not pd.api.types.is_datetime64_any_dtype(df[time_column]):
        df[time_column] = pd.to_datetime(df[time_column])

    # Filter by group if specified
    if group_by and group_by_values:
        df = df[df[group_by].isin(group_by_values)]

    # Calculate time-based metrics
    time_metrics = {}
    if 'mean' in metrics:
        time_metrics['mean'] = df[value_column].mean()
    if 'std' in metrics:
        time_metrics['std'] = df[value_column].std()
    if 'min' in metrics:
        time_metrics['min'] = df[value_column].min()
    if 'max' in metrics:
        time_metrics['max'] = df[value_column].max()
    if 'range' in metrics:
        time_metrics['range'] = time_metrics.get('max', df[value_column].max()) - time_metrics.get('min', df[value_column].min())
    if 'cv' in metrics:
        time_metrics['cv'] = time_metrics.get('std', df[value_column].std()) / time_metrics.get('mean', df[value_column].mean())

    # Calculate statistical metrics
    if 'skewness' in metrics:
        time_metrics['skewness'] = df[value_column].skew()
    if 'kurtosis' in metrics:
        time_metrics['kurtosis'] = df[value_column].kurtosis()

    # Calculate time-based metrics
    time_based_metrics = {}
    if group_by:
        for group_value in group_by_values:
            group_df = df[df[group_by] == group_value]
            time_based_metrics[group_value] = {
                'mean': group_df[value_column].mean(),
                'std': group_df[value_column].std(),
                'min': group_df[value_column].min(),
                'max': group_df[value_column].max(),
                'range': group_df[value_column].max() - group_df[value_column].min(),
                'cv': group_df[value_column].std() / group_df[value_column].mean(),
                'skewness': group_df[value_column].skew(),
                'kurtosis': group_df[value_column].kurtosis()
            }

    return {
        'overall_metrics': time_metrics,
        'group_metrics': time_based_metrics if group_by else None
    } 