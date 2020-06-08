## the-network-untangling-problem
Implementations for the experiments included in the paper:

Polina Rozenshtein, Nikolaj Tatti and Aristides Gionis. "The network-untangling problem: From interactions to activity timelines."

= list of used packages =
* https://pypi.python.org/pypi/networkx/1.11

### `single_interval` directory for case k=1

* `inner_point.py` implements algorithm Inner. For details see function `runInner`.
* `budget.py` implements algorithm Budget. For details see function `runBudget`.
* `baseline.py` implements a greedy baseline. For details see function `baseline`.

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

#### `other_experiments`

* `maximal_convergence.py` tests how the solution by `Inner` algorithm evolves during iterations with re-initialization (iterations of `Maximal` algorithm). The results here are reported based on a single run on a synthetic dataset.
```
usage: maximal_convergence.py [-h] [--intlen INTLEN] [--overlap OVERLAP]
                                [--nnodes NNODES]

optional arguments:
  --intlen INTLEN    length of each active interval (default: 100)
  --overlap OVERLAP  activity intervals overlap parameter, between 0 and 1
                     (default: 0.5)
  --nnodes NNODES    number of nodes in the graph (default: 100)
```
* `different_overlap.py` tests both algorithms on the synthetic dataset with varying overlap parameter. The results here are reported based on a single run on a synthetic dataset.
```
usage: test_overlap.py [-h] [--intlen INTLEN] [--nnodes NNODES]

optional arguments:
  --intlen INTLEN  length of each active interval (default: 100)
  --nnodes NNODES  number of nodes in the graph (default: 100)
```

### `multiple_intervals` directory for case k>1

* `k_inner_point.py` implements algorithm k-Inner. For details see function `runKInner`.
* `k_budget.py` implements algorithm k-Budget. For details see function `runKBudget`.
* `k_baseline.py` implements a greedy k-Baseline. For details see function `kbaseline`.

* `k_quick_test_synthetic.py` runs any of 3 algorithms on a synthetic dataset and reports statistics of the solution.
```
usage: k_quick_test_synthetic.py [-h] [-k K] [--intlen INTLEN]
                                 [--overlap OVERLAP] [--nnodes NNODES]
                                 algorithm

positional arguments:
  algorithm          baseline, budget, or inner

optional arguments:
  -h, --help         show this help message and exit
  -k K               number of intervals (default: 10)
  --intlen INTLEN    length of each active interval (default: 10)
  --overlap OVERLAP  activity intervals overlap parameter, between 0 and 1
                     (default: 0.5)
  --nnodes NNODES    number of nodes in the graph (default: 10)
```

#### `other_experiments`
* `effect_of_k.py` runs all three algorithms (`k-Inner`, `k-Budget`, `k-Baseline`) on a synthesic dataset with different `k`
* `k_budget_convergence.py` runs `k-Budget` algorithm and reports its performance on iterative search of inactive points 
* `k_inner_convergence.py` runs `k-Inner` algorithm and reports its performance on iterative search of active points 
* `sensitivity_to_initialization.py` runs `k-Inner` and `k-Budget` algorithm for different percent of randomly selected initial active/inactive points.

All results here are reported based on a single run on a synthetic dataset.


Please see `--help` for arguments 

