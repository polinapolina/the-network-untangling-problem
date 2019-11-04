# the-network-untangling-problem
Implementations for experiments included in the paper:

Polina Rozenshtein, Nikolaj Tatti and Aristides Gionis. "The network-untangling problem: From interactions to activity timelines."

= list of used packages =
* https://pypi.python.org/pypi/networkx/1.11

= Case k=1 =

* `inner_point.py` implements algorithm Inner. For details see function `runInner`.
* `budget.py` implements algorithm Budget. For details see function `runBudget`.
* `baseline.py` implements a greedy baseline. For details see function `baseline`.

= experiments =
* `quick_test_synthetic.py` runs any of 3 algorithms on a synthetic dataset and reports statistics of the solution.
```
usage: quick_test_synthetic.py [-h] [--intlen INTLEN] [--overlap OVERLAP]
                               [--nnodes NNODES]
                               algorithm

positional arguments:
  algorithm          baseline, budget, or inner

optional arguments:
  -h, --help         show this help message and exit
  --intlen INTLEN    length of each active interval (default: 100)
  --overlap OVERLAP  activity intervals overlap parameter, between 0 and 1
                     (default: 0.2)
  --nnodes NNODES    number of nodes in the graph (default: 100)
```

* `test_iterations_inner.py` tests how the solution by `Inner Point` algorithm evolves during iterations with re-initialization
```
usage: test_iterations_inner.py [-h] [--intlen INTLEN] [--overlap OVERLAP]
                                [--nnodes NNODES]

optional arguments:
  --intlen INTLEN    length of each active interval (default: 100)
  --overlap OVERLAP  activity intervals overlap parameter, between 0 and 1
                     (default: 0.5)
  --nnodes NNODES    number of nodes in the graph (default: 100)
```
* `test_overlap.py` tests both algorithms on the synthetic dataset with varying overlap parameter
```
usage: test_overlap.py [-h] [--intlen INTLEN] [--nnodes NNODES]

optional arguments:
  --intlen INTLEN  length of each active interval (default: 100)
  --nnodes NNODES  number of nodes in the graph (default: 100)
```


