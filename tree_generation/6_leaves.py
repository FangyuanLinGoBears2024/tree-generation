import networkx as nx
import matplotlib.pyplot as plt
from itertools import permutations, combinations
from collections import Counter

def generate_non_isomorphic_trees(n):
    # Generate all non-isomorphic trees of order n
    return list(nx.nonisomorphic_trees(n))

def count_leaves(tree):
    # Count the number of leaves (nodes with degree 1) in a tree
    return sum(1 for node in tree.nodes if tree.degree[node] == 1)

def has_no_degree_two(tree):
    # Check if no vertex in the tree has a degree of 2
    return all(tree.degree[node] != 2 for node in tree.nodes)

def filter_trees_with_leaves_and_degree_constraint(trees, leaf_count):
    # Filter trees with a specific number of leaves and no vertex of degree 2
    return [tree for tree in trees if count_leaves(tree) == leaf_count and has_no_degree_two(tree)]

def find_leaves(tree):
    # Find all leaves in a tree (nodes with degree 1)
    return [node for node in tree.nodes if tree.degree[node] == 1]

def compute_subtree_length(tree, selected_leaves):
    # Find the induced subgraph that includes all nodes along the paths connecting the selected leaves
    induced_nodes = set(selected_leaves)
    for u, v in combinations(selected_leaves, 2):
        path = nx.shortest_path(tree, source=u, target=v)  # Get the shortest path between any two leaves
        induced_nodes.update(path)  # Add all nodes on this path to the induced set
    
    induced_subgraph = tree.subgraph(induced_nodes)
    mst = nx.minimum_spanning_tree(induced_subgraph)
    return mst.size(weight=None)  # Return the number of edges in the MST

def compute_distribution_over_permutations(tree):
    # Find all leaves in the tree
    leaves = find_leaves(tree)

    if len(leaves) < 6:
        return None  # Not enough leaves to compute the distribution

    # We are interested in all permutations of the first 6 leaves
    selected_leaves = leaves[:6]
    all_permutations = list(permutations(selected_leaves))

    # To store the counts of all length sequences
    length_sequences = Counter()

    # Compute the lengths for each permutation
    for perm in all_permutations:
        # Compute the lengths of the subtrees for k = 2 to 6
        subtree_lengths = []
        for k in range(2, 7):
            length = compute_subtree_length(tree, perm[:k])
            subtree_lengths.append(length)
        
        # Convert to tuple to store as a hashable key
        length_tuple = tuple(subtree_lengths)
        length_sequences[length_tuple] += 1

    # Compute the count of each length tuple (distribution)
    return length_sequences

def compare_all_tree_pairs(trees):
    # Precompute distributions for all trees
    distributions = []
    for i, tree in enumerate(trees):
        print(f"Computing distribution for Tree {i+1}...")
        dist = compute_distribution_over_permutations(tree)
        distributions.append(dist)
    
    # Compare every pair of trees to check if they have identical distributions
    num_trees = len(trees)
    duplicate_pairs = []

    for i in range(num_trees):
        for j in range(i + 1, num_trees):
            print(f"Comparing Tree {i+1} and Tree {j+1}...")
            dist_i = distributions[i]
            dist_j = distributions[j]
            print(f"FIRST distribution: {dist_i}\n")
            print(f"SECOND distribution: {dist_j}\n")

            # Compare distributions
            if dist_i == dist_j:
                print(f"Trees {i+1} and {j+1} have the same distribution.")
                duplicate_pairs.append((i, j, dist_i))
    
    return duplicate_pairs

def visualize_trees(trees, tree_indices):
    # Visualize the specified trees
    for index in tree_indices:
        plt.figure(figsize=(4, 4))
        nx.draw(trees[index], with_labels=True, node_color='lightblue', font_weight='bold')
        plt.title(f"Tree {index+1}")
        plt.show()

# Generate all non-isomorphic trees for orders 7 to 11 (excluding order 1)
trees_order_7 = generate_non_isomorphic_trees(7)
trees_order_8 = generate_non_isomorphic_trees(8)
trees_order_9 = generate_non_isomorphic_trees(9)
trees_order_10 = generate_non_isomorphic_trees(10)
trees_order_11 = generate_non_isomorphic_trees(11)
trees_order_12 = generate_non_isomorphic_trees(12)

# Filter trees with exactly 6 leaves and no vertex of degree 2
trees_7_with_6_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_7, 6)
trees_8_with_6_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_8, 6)
trees_9_with_6_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_9, 6)
trees_10_with_6_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_10, 6)
trees_11_with_6_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_11, 6)
trees_12_with_6_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_12, 6)
filtered_trees_count = {
    7: len(trees_7_with_6_leaves),
    8: len(trees_8_with_6_leaves),
    9: len(trees_9_with_6_leaves),
    10: len(trees_10_with_6_leaves),
    11: len(trees_11_with_6_leaves),
    12: len(trees_12_with_6_leaves)
}
print("Number of non-isomorphic trees with 7 leaves for each order:")
for order, count in filtered_trees_count.items():
    print(f"Order {order}: {count} trees")

# Combine all trees to be checked
all_trees_with_6_leaves = (
    trees_7_with_6_leaves +
    trees_8_with_6_leaves +
    trees_9_with_6_leaves +
    trees_10_with_6_leaves +
    trees_11_with_6_leaves +
    trees_11_with_6_leaves
)

# Compare all pairs of trees
duplicates = compare_all_tree_pairs(all_trees_with_6_leaves)

if duplicates:
    print(f"\nFound {len(duplicates)} pairs of trees with duplicate distributions.")
    for dup in duplicates:
        tree1_index, tree2_index, distribution = dup
        print(f"\nDuplicate pair: Tree {tree1_index+1} and Tree {tree2_index+1}")
        print(f"Distribution for both trees: {distribution}\n")
        visualize_trees(all_trees_with_6_leaves, [tree1_index, tree2_index])
else:
    print("No duplicate distributions found among the generated trees.")
