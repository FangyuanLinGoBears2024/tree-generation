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

    # We are interested in all permutations of the first 6 leaves
    selected_leaves = leaves[:4]
    all_permutations = list(permutations(selected_leaves))

    # To store the counts of all length sequences
    length_sequences = Counter()

    # Compute the lengths for each permutation
    for perm in all_permutations:
        # Compute the lengths of the subtrees for k = 2 to 6
        subtree_lengths = []
        for k in range(2, 5):
            length = compute_subtree_length(tree, perm[:k])
            subtree_lengths.append(length)
        
        # Convert to tuple and count it
        length_sequences[tuple(subtree_lengths)] += 1

    # Compute the uniform probability distribution
    total_permutations = len(all_permutations)
    joint_prob_dist = {seq: count for seq, count in length_sequences.items()}

    return joint_prob_dist

def compare_all_tree_pairs(trees):
    # Compare every pair of trees to check if they have identical distributions
    num_trees = len(trees)
    duplicate_pairs = []

    for i in range(num_trees):
        for j in range(i + 1, num_trees):
            print(f"Comparing Tree {i+1} and Tree {j+1}...")
            dist_i = compute_distribution_over_permutations(trees[i])
            dist_j = compute_distribution_over_permutations(trees[j])
            print(f"FIRST distribution: {dist_i}\n")
            print(f"SECOND distribution : {dist_j}\n")
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
trees_order_4 = generate_non_isomorphic_trees(4)
trees_order_5 = generate_non_isomorphic_trees(5)
trees_order_6 = generate_non_isomorphic_trees(6)
trees_order_7 = generate_non_isomorphic_trees(7)
trees_order_8 = generate_non_isomorphic_trees(8)

# Filter trees with exactly 6 leaves and no vertex of degree 2
trees_4_with_4_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_4, 4)
trees_5_with_4_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_5, 4)
trees_6_with_4_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_6, 4)
trees_7_with_4_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_7, 4)
trees_8_with_4_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_8, 4)

# Combine all trees to be checked
all_trees_with_4_leaves = (
    trees_4_with_4_leaves +
    trees_5_with_4_leaves +
    trees_6_with_4_leaves +
    trees_7_with_4_leaves +
    trees_8_with_4_leaves
)

# Compare all pairs of trees
duplicates = compare_all_tree_pairs(all_trees_with_4_leaves)
visualize_trees(all_trees_with_4_leaves, [0,1])
if duplicates:
    print(f"\nFound {len(duplicates)} pairs of trees with duplicate distributions.")
    for dup in duplicates:
        tree1_index, tree2_index, distribution = dup
        print(f"\nDuplicate pair: Tree {tree1_index+1} and Tree {tree2_index+1}")
        print(f"Distribution for both trees: {distribution}\n")
        visualize_trees(all_trees_with_4_leaves, [tree1_index, tree2_index])
else:
    print("No duplicate distributions found among the generated trees.")
