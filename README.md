# the-network-untangling-problem
Implementations for experiments included in the paper:

Polina Rozenshtein, Nikolaj Tatti and Aristides Gionis. "The network-untangling problem: From interactions to activity timelines."

= list of used packages =
* https://pypi.python.org/pypi/networkx/1.11

= Case k=1 =
inner_point.py implements algorithm Inner
budget.py implements algorithm Budget
baseline.py implements a greedy baseline

= experiments =
* quick_test_synthetic.py runs all three algorithms on synthetic dataset.
* test_iterations_inner.py tests how the solution by Inner Point algorithm evolves during iterations with re-initialization
* test_overlap.py tests both algorithms on the synthetic dataset with varying overlap parameter


