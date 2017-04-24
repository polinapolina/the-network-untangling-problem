# the-network-untangling-problem
Implementations for experiments included in the paper:

Polina Rozenshtein, Nikolaj Tatti and Aristides Gionis. "The network-untangling problem: From interactions to activity timelines."

= list of used packages =
* https://pypi.python.org/pypi/networkx/1.11

= experiments =

* test_budget.py runs Budget algorithm on the synthetic dataset with 100 nodes and 10000 interactions
* test_inner_point.py runs Inner Point algorithm on the synthetic dataset with 100 nodes and 10000 interactions
* test_iterations_inner.py tests how the solution by Inner Point algorithm evolves during iterations with re-initialization
* test_overlap.py tests both algorithms on the synthetic dataset with varying overlap parameter


