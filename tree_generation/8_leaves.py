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

    # We are interested in all permutations of the first 7 leaves
    selected_leaves = leaves[:8]
    all_permutations = list(permutations(selected_leaves))

    # To store the counts of all length sequences
    length_sequences = Counter()

    # Compute the lengths for each permutation
    for perm in all_permutations:
        # Compute the lengths of the subtrees for k = 2 to 7
        subtree_lengths = []
        for k in range(2, 9):
            length = compute_subtree_length(tree, perm[:k])
            subtree_lengths.append(length)
        
        # Convert to tuple and count it
        length_sequences[tuple(subtree_lengths)] += 1

    # Compute the uniform probability distribution
    total_permutations = len(all_permutations)
    joint_prob_dist = {seq: count / total_permutations for seq, count in length_sequences.items()}

    return joint_prob_dist

def find_first_duplicate_pair(trees):
    # Precompute distributions for all trees
    distributions = []
    for i, tree in enumerate(trees):
        print(f"Computing distribution for Tree {i+1}...")
        dist = compute_distribution_over_permutations(tree)
        print(f"Distribution: {dist}\n")
        distributions.append(dist)
    
    # Compare every pair of trees to check if they have identical distributions
    num_trees = len(trees)

    for i in range(num_trees):
        for j in range(i + 1, num_trees):
            print(f"Comparing Tree {i+1} and Tree {j+1}...")
            if distributions[i] == distributions[j]:
                print(f"Duplicate found between Tree {i+1} and Tree {j+1} with distribution:")
                print(distributions[i])
                visualize_trees(trees, [i, j])
                return (i, j, distributions[i])  # Return immediately after finding a duplicate

    print("No duplicate distributions found among the generated trees.")
    return None

def visualize_trees(trees, tree_indices):
    # Visualize the specified trees
    for index in tree_indices:
        plt.figure(figsize=(4, 4))
        nx.draw(trees[index], with_labels=True, node_color='lightblue', font_weight='bold')
        plt.title(f"Tree {index+1}")
        plt.show()

# Generate all non-isomorphic trees for orders 7 to 15
trees_order_7 = generate_non_isomorphic_trees(7)
trees_order_8 = generate_non_isomorphic_trees(8)
trees_order_9 = generate_non_isomorphic_trees(9)
trees_order_10 = generate_non_isomorphic_trees(10)
trees_order_11 = generate_non_isomorphic_trees(11)
trees_order_12 = generate_non_isomorphic_trees(12)
trees_order_13 = generate_non_isomorphic_trees(13)
trees_order_14 = generate_non_isomorphic_trees(14)
trees_order_15 = generate_non_isomorphic_trees(15)

# Filter trees with exactly 7 leaves and no vertex of degree 2
trees_7_with_8_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_7, 8)
trees_8_with_8_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_8, 8)
trees_9_with_8_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_9, 8)
trees_10_with_8_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_10, 8)
trees_11_with_8_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_11, 8)
trees_12_with_8_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_12, 8)
trees_13_with_8_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_13, 8)
trees_14_with_8_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_14, 8)
trees_15_with_8_leaves = filter_trees_with_leaves_and_degree_constraint(trees_order_15, 8)

# Print the number of trees generated for each order
filtered_trees_count = {
    7: len(trees_7_with_8_leaves),
    8: len(trees_8_with_8_leaves),
    9: len(trees_9_with_8_leaves),
    10: len(trees_10_with_8_leaves),
    11: len(trees_11_with_8_leaves),
    12: len(trees_12_with_8_leaves),
    13: len(trees_13_with_8_leaves),
    14: len(trees_14_with_8_leaves),
    15: len(trees_15_with_8_leaves),
}

print("Number of non-isomorphic trees with 7 leaves for each order:")
for order, count in filtered_trees_count.items():
    print(f"Order {order}: {count} trees")

# Combine all trees to be checked
all_trees_with_8_leaves = (
    trees_7_with_8_leaves +
    trees_8_with_8_leaves +
    trees_9_with_8_leaves +
    trees_10_with_8_leaves +
    trees_11_with_8_leaves +
    trees_12_with_8_leaves +
    trees_13_with_8_leaves +
    trees_14_with_8_leaves +
    trees_15_with_8_leaves
)

# Find the first duplicate pair and stop immediately
first_duplicate = find_first_duplicate_pair(all_trees_with_8_leaves)

if first_duplicate:
    tree1_index, tree2_index, distribution = first_duplicate
    print(f"\nFirst duplicate pair found: Tree {tree1_index+1} and Tree {tree2_index+1}")
    print(f"Distribution for both trees: {distribution}\n")
else:
    print("No duplicates found.")
