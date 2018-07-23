import numpy as np
from itertools import combinations
from skbio.stats.ordination import pcoa
from skbio import DistanceMatrix
import matplotlib.pyplot as plt
import pandas as pd
import igraph


def single_distance(index, spn, sample_names, features, css):
    # Get sample names to be compared
    name_i = spn[index, 0]
    name_j = spn[index, 1]

    # Get their indices
    i = int(np.where(sample_names == name_i)[0])
    j = int(np.where(sample_names == name_j)[0])

    feature_union = np.where(features[:, [i, j]].sum(axis=1) > 0)[0]

    css_tmp = css[feature_union[:, None], feature_union]

    a = features[feature_union, i] / np.sum(features[feature_union, i])
    b = features[feature_union, j] / np.sum(features[feature_union, j])

    abt = np.matrix(a).T * np.matrix(b)
    aat = np.matrix(a).T * np.matrix(a)
    bbt = np.matrix(b).T * np.matrix(b)

    cssAB = np.multiply(css_tmp, abt)

    d = np.sum(cssAB) / np.max([np.sum(np.multiply(css_tmp, aat)),
                                np.sum(np.multiply(css_tmp, bbt))])

    return d


def calc_distances(features, css, sample_names):
    assert css.shape[0] == features.shape[0]
    np.fill_diagonal(css.get_values(), 1)
    css = np.matrix(css.get_values())
    spn = np.array(list(combinations(sample_names, 2)))

    distlist = []
    for I in range(spn.shape[0]):
        distlist.append(
            single_distance(index=I, spn=spn, sample_names=sample_names,
                            features=features.get_values()[:, 1:], css=css))

    return distlist


def make_distance_matrix(features_orig, css, weighted=True):
    features = features_orig.copy()
    sample_names = features.columns[1:].values

    if not weighted:
        # Convert to absence/presence
        pd.options.mode.chained_assignment = None
        subset = features.loc[:, features.columns[1:]]
        subset[subset.loc[:, :] > 0] = 1
        features.loc[:, features.columns[1:]] = subset

    dist_list = calc_distances(features, css, sample_names)

    matrix = np.zeros([len(sample_names), len(sample_names)])
    value_index = -1
    for j in range(len(sample_names)):
        for i in range(len(sample_names)):
            if i > j:
                value_index += 1
                matrix[i, j] = dist_list[value_index]
    return pd.DataFrame(data=matrix, index=sample_names, columns=sample_names)


def make_dm_symmetric(matrix_orig):
    matrix = matrix_orig.copy().get_values()
    sym_matrix = matrix_orig.copy()
    sample_names = sym_matrix.columns.tolist()
    sym_matrix = sym_matrix.get_values()
    n = sym_matrix.shape[0]
    for j in range(n):
        for i in range(n):
            if i > j:
                sym_matrix[j, i] = matrix[i, j]

    assert (sym_matrix == sym_matrix.T).all()

    return pd.DataFrame(sym_matrix, sample_names, sample_names)


def perform_pcoa(sym_matrix, title, show=True, dotsize=5, subspecies=None):
    metadata = {}
    sample_names = sym_matrix.columns.tolist()

    if subspecies is None:
        subspecies = {}
        plotby = 'concentration'
        for i in range(sym_matrix.shape[0]):
            subspecies[sample_names[i]] = None
    else:
        plotby = 'subspecies'

    for i in range(1, sym_matrix.shape[0] + 1):
        metadata[i] = {'name': sample_names[i - 1],
                       'concentration': sample_names[i - 1][0:1],
                       'subspecies': subspecies[sample_names[i - 1]]}

    df = pd.DataFrame.from_dict(metadata, orient='index')

    dm = DistanceMatrix(sym_matrix,
                        list(range(1, sym_matrix.shape[0] + 1)))

    pcoa_results = pcoa(dm)

    plt.clf()
    fig = pcoa_results.plot(df, plotby, title=title, s=dotsize)

    if show:
        fig.show()
    else:
        fig.savefig(title, dpi=1000, bbox_inches='tight')


def extract_css_matrix(edges_df):
    # Extract relevant columns
    edge_list = edges_df.get_values()[:, [0, 1, 4]]

    # Make graph object
    g = igraph.Graph()

    # Get list of feature names
    full_names_list = edge_list[:, [0]].astype(int).astype(
        str).flatten().tolist() + \
                      edge_list[:, [1]].astype(int).astype(
                          str).flatten().tolist()

    # Remove redundancies
    node_names = []
    seen = set()
    for name in full_names_list:
        if name not in seen:
            node_names.append(name)
            seen.add(name)

    # Sort feature names
    node_names = sorted(node_names)

    # Add nodes to graph
    g.add_vertices(node_names)

    # Make tuples with edges
    edge_tuples = []
    for pair in edge_list[:, [0, 1]].tolist():
        edge_tuples.append((str(int(pair[0])), str(int(pair[1]))))

    # Add edges to graph
    g.add_edges(edge_tuples)

    # Add weights to graph
    g.es['weight'] = 1
    g.es['weight'] = edge_list[:, 2].tolist()

    # Extract adjacency matrix
    matrix = g.get_adjacency(attribute='weight')
    matrix = np.array(matrix.data)

    # Remove values below 0.6
    matrix[matrix < 0.6] = 0

    matrix = pd.DataFrame(data=matrix, index=node_names, columns=node_names)

    return matrix


def similarity2dissimilarity(input_df):
    ones = np.ones(input_df.shape)
    new_values = np.array(ones - input_df)
    np.fill_diagonal(new_values, 0)
    out_df = pd.DataFrame(new_values, input_df.columns, input_df.index)

    return out_df
