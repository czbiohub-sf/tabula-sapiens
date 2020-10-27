
import itertools
from functools import partial
import multiprocessing
import sys
import tempfile
import time

import numpy as np
import pandas as pd

_quiet = False
_debug = False


def set_quiet(val, print_debug=False):
    global _quiet, _debug
    _quiet = bool(val)
    _debug = bool(print_debug)


# Cribbed from https://github.com/dib-lab/sourmash/blob/c7ed8bef7de9b1581b9b067517f64ac64f31f9d0/sourmash/logging.py
def notify(s, *args, **kwargs):
    "A simple logging function => stderr."
    if _quiet:
        return
    print(u'\r\033[K', end=u'', file=sys.stderr)
    print(s.format(*args, **kwargs), file=sys.stderr,
          end=kwargs.get('end', u'\n'))
    if kwargs.get('flush'):
        sys.stderr.flush()


def process_splicing(filename, gene='geneR1A_uniq', cell='cell', tissue=None, method='10x'):
    splicing_df = pd.read_parquet(filename)

    # Add cell ids that match the h5ad object
    if method == '10x':
        # regex builder: https://regex101.com/r/rE1xfN/1
        matches = splicing_df[cell].str.extractall('(?P<channel_cleaned>[\w-]+)_S\d+_L\d+_(?P<cell_barcode>[ACGT]+)')
        matches = matches.droplevel(-1)
        splicing_df = pd.concat([splicing_df, matches], axis=1)
        splicing_df['cell_id'] = splicing_df['cell_barcode'].astype(str) + '_' + splicing_df['channel_cleaned'].astype(str)
    elif method == 'ss2_tsp1':
        splicing_df['cell_id'] = splicing_df['cell'].astype(str) + '.homo.gencode.v30.ERCC.chrM'
    elif method == "ss2_tsp2":
        splicing_df['cell_id'] = splicing_df['cell']

    # Drop duplicate cell ids and gene names
    splicing_df_no_dups = splicing_df.drop_duplicates(['cell_id', gene])

    if tissue is not None:
        splicing_df_no_dups = splicing_df_no_dups.query('tissue == @tissue')

    # Don't use rows with empty gene names -- these are unannotated genes
    splicing_df_no_dups = splicing_df_no_dups.query(f'{gene} != ""')
    print(splicing_df_no_dups.shape)
    splicing_df_no_dups.head()
    # splicing2d = splicing_df_no_dups.pivot(index='cell_id', columns=gene, values='z')
    return splicing_df_no_dups


def my_nan_euclidean_metric(row_i, row_j):
    #     assert row_i.shape == row_j.shape

    i_missing = np.isnan(row_i)
    j_missing = np.isnan(row_j)

    shared = (~i_missing) & (~j_missing)
    n_shared = shared.sum()
    if n_shared == 0:
        return 0
    weight = row_i.shape[0] / n_shared

    i_shared = row_i[shared]
    j_shared = row_j[shared]

    sum_of_squares = np.sum(np.square(i_shared - j_shared))
    distance = np.sqrt(weight * sum_of_squares)
    return distance


def my_nan_manhattan_metric(row_i, row_j):
    #     assert row_i.shape == row_j.shape

    i_missing = np.isnan(row_i)
    j_missing = np.isnan(row_j)

    shared = (~i_missing) & (~j_missing)
    n_shared = shared.sum()
    if n_shared == 0:
        return 0

    i_shared = row_i[shared]
    j_shared = row_j[shared]
    distance = np.sum(np.absolute(i_shared - j_shared))
    return distance


def to_memmap(array):
    """Write a memory mapped array
    Create a memory-map to an array stored in a binary file on disk.
    Memory-mapped files are used for accessing small segments of
    large files on disk, without reading the entire file into memory.
    :param np.array array to memory map
    :return: np.array large_memmap memory mapped array
    :return: str filename name of the file that memory mapped array is written to
    """
    import numpy as np

    filename = tempfile.NamedTemporaryFile(
        prefix="array", suffix=".mmap", delete=False
    ).name
    shape = array.shape
    f = np.memmap(filename, mode="w+", shape=shape, dtype=array.dtype)
    f[:] = array[:]
    del f
    large_memmap = np.memmap(filename, dtype=array.dtype, shape=shape)
    return large_memmap, filename


def distance_args_unpack(args, metric):
    """Helper function to unpack the arguments. Written to use in pool.imap
    as it can only be given one argument."""
    row_i, row_j = args
    if metric == "euclidean":
        return my_nan_euclidean_metric(row_i, row_j)
    elif metric == "manhattan":
        return my_nan_manhattan_metric(row_i, row_j)


def get_distances_at_index(index, matrix, metric):
    """Returns similarities of all the combinations of signature at index in
    the siglist with the rest of the indices starting at index + 1. Doesn't
    redundantly calculate signatures with all the other indices prior to
    index - 1

    :param int index: generate masks from this image
    :param metric: distance metric l1 (manhattan) or l2 (euclidean)
    :param boolean ignore_abundance
        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate the angular
        similarity.
    :param boolean downsample by max_hash if True
    :param siglist list of signatures
    :return: list of similarities for the combinations of signature at index
        with rest of the signatures from index+1
    """
    startt = time.time()
    sig_iterator = itertools.product([matrix[index, :]], matrix[(index + 1) :, :])
    func = partial(distance_args_unpack, metric=metric)
    similarity_list = list(map(func, sig_iterator))
    notify(
        "comparison for index {} done in {:.5f} seconds",
        index,
        time.time() - startt,
        end="\r",
    )
    return similarity_list


def distances_parallel(matrix, n_jobs, metric="euclidean"):
    """Compare all combinations of signatures and return a matrix
    of similarities. Processes combinations parallely on number of processes
    given by n_jobs

    :param list siglist: list of signatures to compare
    :param metric: distance metric l1 (manhattan) or l2 (euclidean)
    :param boolean ignore_abundance
        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate the angular
        similarity.
    :param boolean downsample by max_hash if True
    :param int n_jobs number of processes to run the similarity calculations on
    :return: np.array similarity matrix
    """

    # Starting time - calculate time to keep track in case of lengthy siglist
    start_initial = time.time()

    # Create a memory map of the siglist using numpy to avoid memory burden
    # while accessing small parts in it
    matrix, _ = to_memmap(np.array(matrix))
    notify("Created memmapped input matrix")

    # Check that length of combinations can result in a square similarity matrix
    length_matrix = len(matrix)
    shape = (length_matrix, length_matrix)

    # Initialize with ones in the diagonal as the similarity of a signature with
    # itself is one
    distances = np.zeros(shape, dtype=np.float64)
    memmap_distances, filename = to_memmap(distances)
    notify("Initialized memmapped similarities matrix")

    # Initialize the function using func.partial with the common arguments like
    # siglist, ignore_abundance, downsample, for computing all the signatures
    # The only changing parameter that will be mapped from the pool is the index
    func = partial(get_distances_at_index, matrix=matrix, metric=metric)
    notify("Created similarity func")

    # Initialize multiprocess.pool
    pool = multiprocessing.Pool(processes=n_jobs)

    # Calculate chunk size, by default pool.imap chunk size is 1
    chunksize, extra = divmod(length_matrix, n_jobs)
    if extra:
        chunksize += 1
    notify("Calculated chunk size for multiprocessing")

    # This will not generate the results yet, since pool.imap returns a generator
    result = pool.imap(func, range(length_matrix), chunksize=chunksize)
    notify("Initialized multiprocessing pool.imap")

    # Enumerate and calculate similarities at each of the indices
    # and set the results at the appropriate combination coordinate
    # locations inside the similarity matrix
    for index, l in enumerate(result):
        startt = time.time()
        col_idx = index + 1
        for idx_condensed, item in enumerate(l):
            i = index
            j = col_idx + idx_condensed
            memmap_distances[i, j] = memmap_distances[j, i] = item
        notify(
            "Setting similarities matrix for index {} done in {:.5f} seconds",
            index,
            time.time() - startt,
            end="\r",
        )
    notify("Setting similarities completed")

    pool.close()
    pool.join()

    notify(
        "Time taken to compare all pairs parallely is {:.5f} seconds ",
        time.time() - start_initial,
    )
    return np.memmap(filename, dtype=np.float64, shape=shape)


def make2d(splicing_tidy, index='cell_id', values='scaled_z', columns='geneR1A_uniq'):
    data2d = splicing_tidy.pivot(index=index, values=values, columns=columns)
    return data2d


def compute_distances_df(splicing2d, n_jobs, metric="euclidean"):
    dists = distances_parallel(splicing2d, n_jobs=n_jobs, metric=metric)
    dists_df = pd.DataFrame(dists, index=splicing2d.index, columns=splicing2d.index)
    return dists_df
