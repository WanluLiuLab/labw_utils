# Readme for `assess_interval_engine`

Here contains code that assesses data structure that indexes intervals. We perform this research because we need to find overlapping intervals from two sets of intervals.

Conclusion: Use the Numpy one. Pandas and NumExpr implementations are slow. Although `intervaltree` implementation is ultra-fast in small data size, it becomes incredibly slow if data size exceeds expectations (cannot even finish reading it!).
