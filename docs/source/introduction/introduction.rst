.. _sec_introduction:

============
Introduction
============

This is the documentation for slime, the local ancestry simulator.

.. note:: This documentation is incomplete and under development.

Under the hood, each ``slime`` simulation consists of 3 components:

#. A forward-in-time simulation of admixture using SLiM.
#. A backwards-in-time simulation of ancestral population history using msprime.
#. Extraction of local ancestry from the resulting tree sequence.