import numpy as np
from cutnorm import compute_cutnorm, tools

print("Erdos Renyi Random Graph Cutnorms...")
# Generate Erdos Renyi Random Graph
n = 100
p = 0.5
erdos_renyi_a = tools.sbm.erdos_renyi(n, p)
erdos_renyi_b = tools.sbm.erdos_renyi(n, p)

# Compute l1 norm
normalized_diff = (erdos_renyi_a - erdos_renyi_b) / n**2
l1 = np.linalg.norm(normalized_diff.flatten(), ord=1)

# Compute cutnorm
cutn_round, cutn_sdp, info = compute_cutnorm(erdos_renyi_a, erdos_renyi_b)

print("l1 norm: ", l1)
print("cutnorm rounded: ", cutn_round)
print("cutnorm sdp: ", cutn_sdp)

print("\nSynthetic graphs with ones in center...")
n = 1000
a = np.ones((n, n))
a_small = np.ones((n // 2, n // 2))
b = np.ones((n, n))
b[n // 3:-n // 3, n // 3:-n // 3] = 0

print("Opt obj function val: " + str(np.sum(a - b) / np.sum(a)))

# Same-dim, not weighted
objf_lower, objf_upper, info = compute_cutnorm(a, b)
print("Same Dim, not weighted")
print(objf_lower, objf_upper)

# Diff-dim, not weighted
objf_lower, objf_upper, info = compute_cutnorm(a_small, b)
print("Diff Dim, not weighted")
print(objf_lower, objf_upper)

# Diff-dim Weighted
objf_lower, objf_upper, info = compute_cutnorm(a_small, b,
                                               np.ones(n // 2) / (n // 2),
                                               np.ones(n) / n)
print("Diff Dim, weighted")
print(objf_lower, objf_upper)

# Same-dim, not weighted, logn
objf_lower, objf_upper, info = compute_cutnorm(
    a, b, logn_lowrank=True, extra_info=True)
print("Same Dim, not weighted, lowrank")
print(objf_lower, objf_upper)
print("Debug Keys: ", info.keys())
