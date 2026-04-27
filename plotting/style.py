"""
Shared plotting style constants and helper functions.
"""
# Accessible color palette (first 4 preserved from the original project palette).
DEFAULT_COLORS = [
    '#0072B2',  # blue
    '#009E73',  # green
    '#D55E00',  # vermilion
    '#E69F00',  # amber
    '#CC79A7',  # pink
    '#56B4E9',  # sky blue
    '#F0E442',  # yellow
    '#000000',  # black
]

DEFAULT_DPI = 400
DEFAULT_FIGSIZE = (12, 6)
DEFAULT_FACECOLOR = '#f7f7fa'


def get_color_cycle(n=None):
    """
    Returns a list of colors. If n is specified and exceeds the palette size,
    cycles through the palette.
    """
    if n is None:
        return DEFAULT_COLORS[:]
    return [DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i in range(n)]


def apply_style(ax, remove_spines=None):
    """
    Apply shared style settings to a matplotlib Axes object.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to style.
    remove_spines : list of str, optional
        Which spines to remove. Default: ['top', 'right'].
    """
    if remove_spines is None:
        remove_spines = ['top', 'right']

    ax.grid(True, linestyle='--', alpha=0.35)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_facecolor(DEFAULT_FACECOLOR)

    for spine in remove_spines:
        ax.spines[spine].set_visible(False)


def prettify_label(text):
    """
    Sanitise a string for display by replacing underscores with spaces.

    Parameters
    ----------
    text : str
        Raw label text (may contain underscores from filenames or config keys).

    Returns
    -------
    str
        Human-friendly label with spaces instead of underscores.
    """
    return text.replace('_', ' ')


def format_label_with_stats(label, values):
    """
    Creates a legend label with mean and std statistics embedded.

    Parameters
    ----------
    label : str
        Base label text.
    values : array-like
        Data values to compute statistics from.

    Returns
    -------
    str
        Formatted label like "label (avg=X.XX +/- Y.YY)"
    """
    import numpy as np
    label = prettify_label(label)
    mean_val = np.mean(values)
    std_val = np.std(values)
    return f"{label} (avg={mean_val:.2f} \u00b1 {std_val:.2f})"


def format_selection_subtitle(selection_str, max_length=60):
    """
    Format a selection string for display as a plot subtitle.

    Parameters
    ----------
    selection_str : str
        MDAnalysis selection string.
    max_length : int, optional
        Maximum character length before truncating.

    Returns
    -------
    str
        Formatted selection string, truncated with '…' if needed.
    """
    if len(selection_str) <= max_length:
        return selection_str
    return selection_str[:max_length - 1] + '…'


def format_selection_context(target_selection=None, ref_selection=None,
                             group_selection=None, max_length=120):
    """Build a compact context line with explicit selection strings."""
    parts = []
    if target_selection:
        parts.append(f"target={target_selection}")
    if ref_selection:
        parts.append(f"ref={ref_selection}")
    if group_selection:
        parts.append(f"group={group_selection}")

    if not parts:
        return ''

    context = ' | '.join(parts)
    return format_selection_subtitle(context, max_length=max_length)
