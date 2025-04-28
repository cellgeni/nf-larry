import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numba as nb
from scipy.spatial.distance import squareform

def sequences_to_numeric(seqs, mapping=None):
    """
    Convert an array of equal-length sequences to a 2D NumPy array of integers.
    
    By default, A, C, G, T are mapped to 0, 1, 2, 3, and ambiguous base 'N' is mapped to -1.
    Any base not found in the mapping is also treated as ambiguous (-1).
    
    Parameters
    ----------
    seqs : array-like of str
        List or array of sequences.
    mapping : dict, optional
        Mapping from nucleotide characters to integers.
    
    Returns
    -------
    np.ndarray
        2D array of shape (n_sequences, sequence_length) with dtype=np.int8.
    """
    if mapping is None:
        mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': -1}
    numeric_seqs = []
    for seq in seqs:
        # Map each base; if base not found, treat it as ambiguous (-1)
        numeric_seq = [mapping.get(base, -1) for base in seq]
        numeric_seqs.append(numeric_seq)
    return np.array(numeric_seqs, dtype=np.int8)

@nb.njit(parallel=True)
def compute_condensed_hamming_ambiguous(arr):
    """
    Compute the condensed pairwise Hamming distance between sequences, ignoring positions 
    where either sequence has an ambiguous base (encoded as -1).
    
    Parameters
    ----------
    arr : np.ndarray
        2D NumPy array of shape (n_sequences, sequence_length).
    
    Returns
    -------
    np.ndarray
        Condensed distance matrix (1D array) of raw mismatch counts.
    """
    n, L = arr.shape
    m = n * (n - 1) // 2
    dists = np.empty(m, dtype=np.int32)
    idx = 0
    for i in nb.prange(n - 1):
        for j in range(i + 1, n):
            mismatches = 0
            # Compare positions while skipping ambiguous (-1) positions.
            for k in range(L):
                a = arr[i, k]
                b = arr[j, k]
                if a == -1 or b == -1:
                    continue
                if a != b:
                    mismatches += 1
            dists[idx] = mismatches
            idx += 1
    return dists

def compute_hamming_matrix(seqs):
    """
    Compute the full pairwise Hamming distance matrix for an array of equal-length sequences.
    
    Sequences are first converted into a numeric array that handles ambiguous bases.
    The condensed distance matrix is computed using a Numba-accelerated function and then
    converted to a full square matrix.
    
    Parameters
    ----------
    seqs : array-like of str
        List or array of sequences.
    
    Returns
    -------
    np.ndarray
        Full square (n x n) Hamming distance matrix (dtype=np.int32).
    """
    # Convert sequences to numeric form (with -1 for ambiguous bases)
    arr = sequences_to_numeric(seqs)
    # Compute the condensed distance matrix with ambiguous positions skipped.
    condensed = compute_condensed_hamming_ambiguous(arr)
    # Convert condensed array to a full square matrix.
    return squareform(condensed)

def hamming_filter(counts, min_distance=3):
    """
    Filters out barcodes based on Hamming distance.
    
    For each pair of barcodes with a Hamming distance below the threshold, the barcode 
    with the lower count is rejected. In case of ties, both are recorded.
    
    Parameters
    ----------
    counts : pd.Series
        Pandas Series with sequences as the index and counts as values.
    min_distance : int, default 3
        Minimum Hamming distance threshold; pairs with distance below this value 
        are considered for filtering.
    
    Returns
    -------
    tokeep : pd.Index
        Index of sequences that passed the filter.
    toreject : list
        List of sequences to be rejected.
    ties : dict
        Dictionary mapping a sequence to a list of sequences with which it ties.
    """
    print(f"Filtering using Hamming distance threshold of {min_distance}...")
    print("using the optimised version...")
    seqs = counts.index.values.astype(str)
    # Compute pairwise Hamming distances using our custom, accelerated function.
    hdist = compute_hamming_matrix(seqs)
    
    # Plot histogram for the upper triangle (excluding diagonal) of the distance matrix.
    triu_indices = np.triu_indices_from(hdist, k=1)
    plt.hist(hdist[triu_indices], bins=50)
    plt.title("Pairwise Hamming Distance")
    plt.xlabel("Hamming Distance")
    plt.ylabel("Frequency")
    plt.show()
    
    n = len(seqs)
    toreject = []
    ties = {}
    
    for i in range(n):
        seq = seqs[i]
        seq_count = counts[seq]
        # Get indices where the Hamming distance is below the threshold (excluding self)
        below_indices = np.nonzero((hdist[i, :] < min_distance) & (np.arange(n) != i))[0]
        if below_indices.size > 0:
            below_counts = counts.iloc[below_indices]
            max_below = below_counts.max()
            if seq_count < max_below:
                toreject.append(seq)
            elif seq_count == max_below:
                ties[seq] = list(below_counts.index[below_counts == max_below])
    
    tokeep = counts.index[~counts.index.isin(toreject)]
    print(f"Passed: {len(tokeep)}\nRejected: {len(toreject)}\nTies: {len(ties)}")
    return tokeep, toreject, ties

class CBUSeries(pd.Series):
    """
    CBC-Barcode-UMI Series with counts.
    
    A Pandas Series subclass indexed by CBC, Barcode, and UMI, typically containing
    read counts.
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
    
    def filter_by_hamming(self, which="Barcode", min_distance=2):
        """
        Filter based on Hamming distance.
        
        Groups the CBUSeries by the specified index level, computes the pairwise
        Hamming distances (ignoring ambiguous positions), and filters out sequences 
        with distances below the threshold by retaining the sequence with the higher 
        read count. In the case of ties, both sequences are retained.
        
        Parameters
        ----------
        which : str
            The index level to apply the filter on ('CBC', 'Barcode', or 'UMI').
        min_distance : int
            Sequences with Hamming distance smaller than min_distance are rejected.
        
        Returns
        -------
        CBUSeries
            Filtered CBUSeries.
        """
        counts = self.groupby([which]).sum()
        tokeep, toreject, ties = hamming_filter(counts=counts, min_distance=min_distance)
        if ties:
            print(f"Ties detected between: {ties}")
        
        if which == "CBC":
            return self.loc[tokeep, :, :]
        elif which == "Barcode":
            return self.loc[:, tokeep, :]
        elif which == "UMI":
            return self.loc[:, :, tokeep]
        else:
            raise ValueError("Invalid index level specified.")
