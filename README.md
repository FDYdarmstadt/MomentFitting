# MomentFitting
This software package contains a basic implementation of the hierarchical
moment-fitting approach first presented in the following publication:
Müller, Kummer, Oberlack: Highly accurate surface and volume integration on
implicit domains by means of moment-fitting.
DOI: 10.1002/nme.4569
URL: http://onlinelibrary.wiley.com/doi/10.1002/nme.4569/abstract

The implementation considers the construction of quadrature rules for a
single quadrilateral cell K = [-1, 1] x [-1, 1] which is cut by the zero
iso-contour of a level set function.

Use the functions 'getSurfaceRule' and 'getVolumeRule' to obtain a set of
quadrature nodes and weights for the integration over the zero iso-contour and
the positive sub-volume (w.r.t. the level set function), respectively. Use the
Matlab test runner 'runtests' to execute the basic verification tests in
'surfaceTest.m' and 'volumeTest.m', which also serve as usage examples.

If you have any questions or comments, please contact me via
mueller@fdy.tu-darmstadt.de
