# SVGSplines
A Spatially Variable Gene(SVG) detection method that improves SPARK-X. This repo also contains various benchmark results  showing comparison of SPARK-X to SVGSplines.

Software implementation
Two R functions are made available via GitHub:
• spark.x, which reproduces the output of SPARK-X (with option = mixture) from a regression-based perspective; and
• svg.spline, which implements the proposed additive spline approach.

The spark.x function implements two tests: the score test (method = "score", the default) and the Wald test (method = "wald"). Its purpose is to demonstrate that SPARK-X is indeed the score test for a regression model. It is not heavily optimized for computation efficiency and memory use, although it is often more than an order of magnitude faster than SPARK-X when
working on the simulated SVGbench data.
