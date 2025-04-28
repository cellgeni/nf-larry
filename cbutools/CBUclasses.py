import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from .hamming import hamming_filter
from math import ceil
from functools import wraps


def load_barcodes(file):
    """Load barcodes CSV file

    Loads CBC,Barcodes,UMI or CBC, Barcodes CSV files into a CBUSeries or
    CBSeries respectively

    Parameters
    ----------
    file : PATH or str
        Path to the file

    Raises
    ------
    Exception
        Raised if files does not look like CBUSeries or CBseries

    Returns
    -------
    CBUSeries or CBseries

    """

    data = pd.read_csv(file)
    cols = data.columns.tolist()
    if cols == ["CBC", "Barcode", "UMI", "0"]:
        ser = pd.Series(data["0"])
        ind = pd.MultiIndex.from_frame(data.loc[:, ["CBC", "Barcode", "UMI"]])
        ser.index = ind
        return CBUSeries(ser)

    elif cols == ["CBC", "Barcode", "0"]:
        ser = pd.Series(data["0"])
        ind = pd.MultiIndex.from_frame(data.loc[:, ["CBC", "Barcode"]])
        ser.index = ind
        return CBSeries(ser)

    else:
        raise Exception("This does not look like CBUSeries or CBSeries saved")


def filter_series(series, groupby=None, labels=["CBC", "Barcode", "UMI"], min_counts=0):
    """Filter a Series object by counts

    Filters pd.Series, CBUSeries or CBSeries by counts grouped by categories (in the index): CBC, Barcode,
    UMI or combinaion of these categories

    Parameters
    ----------
    series : CBUSeries or CBSeries
        Series to be filtered
    groupby : list or str
        list with combination of categories or a string with a single category
    labels : list
        list of category names (labels) to consider for grouping (correspond to the index)
    min_counts : int
        only entries with >= min_counts are kept

    Raises
    ------
    Exception
        Raised if groupby has invalid arguments

    Returns
    -------
    pd.Series or CBUSeries or CBSeries
        Filtered Series object

    """

    # If none or all three levels are supplied - just filter the values and return
    if (groupby == None) or sorted(groupby) == sorted(labels):
        return series[series >= min_counts]

    # Converting a string into a list
    if type(groupby) is str:
        groupby = [groupby]

    # Catching invalid groupby
    if (type(groupby) is list) and (0 < len(groupby) <= len(labels)):
        # Checking if all groupby elements are valid
        check = [i in labels for i in groupby]
        if not all(check):
            raise Exception(
                "Invalid groupby argument (should be None or combination of f{labels}"
            )
    else:
        raise Exception(
            "Invalid groupby argument (should be None or combination of f{labels}"
        )

    groupby = [x for x in labels if x in groupby]
    filtered = series.groupby(groupby).sum()
    filtered = filtered[filtered >= min_counts]

    # Getting the matching index in the original data (dropping levels that were grouped)
    todrop = [i for i in labels if i not in groupby]

    indexSUB = series.index.droplevel(level=todrop)
    return series[indexSUB.isin(filtered.index)]


def plot_groupby_hist(series, groupby, bins=50, vmax=None, title="", *args, **kwargs):
    """Plot histogram of grouped counts

    Plots a histogram of counts grouped by indicated categories or their
    combinations

    Parameters
    ----------
    series : pd.Series or CBUSeries or CBSeries
        A series object with counts
    groupby : list or str
        list with combination of categories or a string with a single category (should be in the Series index))
    bins : int
        number of bins in the histogram
    vmax : int or float
        maximum number on the X axis
    title : str
        plot title

    """

    grouped = series.groupby(groupby).sum()
    n, bins, patches = plt.hist(np.log10(grouped), bins=bins, *args, **kwargs)

    if vmax is None:
        vmax = grouped.max()
    ticks = ceil(np.log10(vmax)) + 1
    vmax = 10 ** (ticks - 1)
    logpos = np.logspace(0, np.log10(vmax), ticks)
    plt.xticks(range(ticks), logpos)
    # plt.xscale('log')
    plt.yscale("log")
    plt.title(title)
    plt.show()


def plot_barcodes_no(x):
    """Plot barcode number per cell

    Parameters
    ----------
    x : pd.Series or CBUSeries or CBSeries
        Series with number of barcodes per cell

    """

    x = x.sort_values(ascending=False)
    plt.plot(range(len(x.index)), x)
    ax = plt.gca()
    # ax.yaxis.get_major_locator().set_params(integer=True)
    plt.show()


def check_CB_index(f):
    """Decorator - checks if CBSeries index is valid

    Parameters
    ----------
    f : function
        function running on CBSeries that requires valid index

    Raises
    ------
    Exception
        CBSeries index is invalid

    """

    @wraps(f)
    def wrapper(*args, **kwargs):
        if args[0].index.names != ["CBC", "Barcode"]:
            raise Exception('The multiindex is not "CBC", "Barcode"')
        return f(*args, **kwargs)

    return wrapper


class CBSeries(pd.Series):
    """CBC-Barcode Series with counts

    Class inheriting from pd.Series, indexed by CBC and Barcode, typically
    containing UMI counts.

    Parameters
    ----------
    pd.Series : pd.Series
        pandas Series object

    Raises
    ------
    Exception
        Index invalid

    """

    def __init__(self, *args, **kwargs):
        super(CBSeries, self).__init__(*args, **kwargs)
        if not self.index.is_unique:
            raise Exception("Index is not unique. Convert to pd.Series.")

    @property
    def _constructor(self):
        if not self.index.is_unique:
            raise Exception("Index is not unique. Confert to pd.Series.")
        return CBSeries

    @check_CB_index
    def filter_by_UMI(self, groupby="Barcode", min_counts=0):
        """Filter CBSeries by UMI

        Parameters
        ----------
        groupby : list or str
            list with combination of categories or a string with a single
            category (correspond to index)
        min_counts : int
            only entires with >= min_counts are kept

        Returns
        -------
        CBSeries
            Filtered object
        """
        labels = ["CBC", "Barcode"]
        return filter_series(
            self, groupby=groupby, labels=labels, min_counts=min_counts
        )

    def plot_hist(
        self,
        groupby="Barcode",
        title="Number of UMIs per {groupby} combination",
        *args,
        **kwargs,
    ):
        """Plot histogram of counts

        Parameters
        ----------
        groupby : list or str
            List with combination of categories or a string with a sngle category to group counts (correspond to the index)
        title : str
            Plot title
        *args : tuple
            Arguments passed to matplotlib hist function
        **kwargs : dict
            Arguments passed to matplotlib hist function

        """
        plot_groupby_hist(self, groupby=groupby, *args, **kwargs)

    @check_CB_index
    def assign_barcodes(self, dispr_filter=None):
        """Assign barcodes to cell

        Loops through each cell and assigns the barcode(s) to each one

        Parameters
        ----------
        dispr_filter : float
            If None all barcodes are assigned to the cell. Otherwise all
            barcodes < dispr_filter * max_count are rejected. max_count
            correspond to the counts observed for the most abundant barcode

        Returns
        -------
        pd.DataFrame
            DataFrame with barcodes assigned per cell

        """
        df = pd.DataFrame()

        for i in self.index.get_level_values("CBC").unique():
            counts = self[i]
            if dispr_filter is not None:
                max_counts = max(counts)
                counts = counts[counts >= (dispr_filter * max_counts)]
            barcodes = tuple(counts.index.tolist())
            row = pd.DataFrame(
                {"Barcode_tuple": [barcodes], "Barcode_n": len(barcodes)}, index=[i]
            )
            df = pd.concat([df, (row)])

        def catl(x):
            """Function for concatenating strings"""
            return "-".join(x)

        df["Barcode"] = df["Barcode_tuple"].apply(func=catl)
        return df

    @check_CB_index
    def plot_barcode_no(self):
        """Plot number of barcodes per cell"""
        x = self.groupby("CBC").size()
        plot_barcodes_no(x)

    @check_CB_index
    def summary(self):
        """Print summary of observed barcodes"""
        total_umi = self.sum()
        total_CBC = len(self.index.get_level_values("CBC").unique())
        total_barcode = len(self.index.get_level_values("Barcode").unique())

        x = self.groupby("CBC").size()
        av_barcode_per_cell = x.mean()

        print(
            f"Total UMIs: {total_umi}\n"
            f"Total CBCs: {total_CBC}\n"
            f"Total Barcodes: {total_barcode}\n"
            f"Mean barcode per cell: {av_barcode_per_cell}"
        )

    def save_barcodes(self, *args, **kwargs):
        """Save CBSeries to a CSV file

        Saves to a CSV files which can be read with the load_barcodes function
        """
        self.to_csv(*args, **kwargs)


def check_CBU_index(f):
    """Decorator - checks if CBUSeries index is valid

    Parameters
    ----------
    f : function
        function running on CBUSeries that requires valid index

    Raises
    ------
    Exception
        CBUSeries index is invalid

    """

    @wraps(f)
    def wrapper(*args, **kwargs):
        if args[0].index.names != ["CBC", "Barcode", "UMI"]:
            raise Exception('The multiindex is not "CBC", "Barcode", "UMI"')
        return f(*args, **kwargs)

    return wrapper


class CBUSeries(pd.Series):
    """CBC-Barcode-UMI Series with counts

    Class inheriting from pd.Series, indexed by CBC, Barcode and UMI, typically
    containing read counts.

    Parameters
    ----------
    pd.Series : pd.Series
        pandas Series object

    Raises
    ------
    Exception
        Index invalid

    """

    def __init__(self, *args, **kwargs):
        super(CBUSeries, self).__init__(*args, **kwargs)
        self.summary_data = {}
        if not self.index.is_unique:
            raise Exception("Index is not unique (convert to pd.Series)")

    @property
    def _constructor(self):
        if not self.index.is_unique:
            raise Exception("Index is not unique (convert to pd.Series)")
        return CBUSeries

    def plot_hist(
        self,
        groupby="Barcode",
        title="Number of reads per {groupby} combination",
        *args,
        **kwargs,
    ):
        """Plot histogram of counts

        Parameters
        ----------
        groupby : list or str
            List with combination of categories or a string with a sngle category to group counts (correspond to the index)
        title : str
            Plot title
        *args : tuple
            Arguments passed to matplotlib hist function
        **kwargs : dict
            Arguments passed to matplotlib hist function

        """
        plot_groupby_hist(self, groupby=groupby, *args, **kwargs)

    @check_CBU_index
    def plot_barcode_no(self):
        """Plot number of barcodes per cell"""
        x = self.count_UMI()
        x = x.groupby("CBC").size()
        plot_barcodes_no(x)

    @check_CBU_index
    def filter_by_reads(self, groupby="Barcode", min_counts=0):
        """Filter CBUSeries by reads

        Parameters
        ----------
        groupby : list or str
            list with combination of categories or a string with a single
            category (correspond to index)
        min_counts : int
            only entires with >= min_counts are kept

        Returns
        -------
        CBUSeries
            Filtered object
        """
        labels = ["CBC", "Barcode", "UMI"]

        return filter_series(
            self, groupby=groupby, labels=labels, min_counts=min_counts
        )

    @check_CBU_index
    def filter_by_hamming(self, which="Barcode", min_distance=2):
        """Filter based on hamming distance

        Filters CBUSeries based on hamming distance for sequenced in the
        specified index level (which argument). For pairs of sequences with
        hamming distance below the threshold always the sequence with higher
        read count is kept, in case of ties both sequences are retained.

        Parameters
        ----------
        which : str
            Level of the CBU index to apply the filter on
        min_distance : int
            sequences with distance smaller than min_distance are rejected

        Returns
        --------
        CBUseries
            filtered based on hamming distance

        """

        counts = self.groupby([which]).sum()
        tokeep, toreject, ties = hamming_filter(
            counts=counts, min_distance=min_distance
        )
        if len(ties) > 0:
            print(f"Ties detected between: {ties}")

        if which == "CBC":
            return self.loc[tokeep, :, :]

        if which == "Barcode":
            return self.loc[:, tokeep, :]

        if which == "UMI":
            return self.loc[:, :, tokeep]

    @check_CBU_index
    def count_UMI(self):
        """Count UMIs

        Counts UMIs and returns a CBSeries
        """
        return CBSeries(self.groupby(["CBC", "Barcode"]).size())

    @check_CBU_index
    def summary(self):
        """Print summary of observed barcodes"""
        total_reads = self.sum()
        total_umi = len(self.index)
        total_CBC = len(self.index.get_level_values("CBC").unique())
        total_barcode = len(self.index.get_level_values("Barcode").unique())

        x = self.count_UMI()
        x = x.groupby("CBC").size()
        av_barcode_per_cell = x.mean()

        print(
            f"Total reads: {total_reads}\n"
            f"Total UMIs: {total_umi}\n"
            f"Total CBCs: {total_CBC}\n"
            f"Total Barcodes: {total_barcode}\n"
            f"Mean barcode per cell: {av_barcode_per_cell}"
        )

    def save_barcodes(self, *args, **kwargs):
        """Save CBUSeries to a CSV file

        Saves to a CSV files which can be read with the load_barcodes function
        """
        self.to_csv(*args, **kwargs)
