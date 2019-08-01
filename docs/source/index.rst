.. slime documentation master file, created by
   sphinx-quickstart on Tue Jul 30 16:04:48 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to slime's documentation!
=================================

.. toctree::
   :maxdepth: 2
   :caption: Contents 1:
   introduction
   recenthistory

Introduction
============
Under the hood, each ``slime`` simulation consists of 3 components:

#. A forward-in-time simulation of admixture using SLiM. (:ref:`sec_recenthistory`)
#. A backwards-in-time simulation of ancestral population history using msprime
#. Extraction of local ancestry from the resulting tree sequence


Table of contents
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Installation
============

Step 1. Simulation of admixture
===============================
You must have a SLiM script describing the recent history of admixture.
In slime, we provide a script generator that allows you to generate these scripts for
a flexible range of situations in Python.

:ref:`sec_recenthistory`

Quickstart - a simple example
=============================
 
We'll run through a basic ``slime`` simulation involving a simple two-population
example of admixture.

*******************************************
Step 1: Specify recent history of admixture
*******************************************

Suppose our admixture involves a single-pulse of admixture between two ancestral
populations 10 generations in the past. 
We'll specify configurations for each of these ancestral populations as an
``msprime.PopulationConfiguration`` object.
(Please see the `msprime documentation <`tskit documentation <https://tskit.readthedocs.io/en/stable>`_>`_ for more information about these).

	>>> ref0_config = msprime.PopulationConfiguration(sample_size=0, initial_size=10, growth_rate=0)
	>>> ref1_config = msprime.PopulationConfiguration(sample_size=0, initial_size=15, growth_rate=0)

We'll also specify a configuration for the admixed population, as well as an admixture proportion from each of the reference populations.

	>>> adm_config = msprime.PopulationConfiguration(sample_size=5, initial_size=10, growth_rate = 0)
	>>> adm_props = [0.3, 0.7]

By default, ``slime`` will assume a single-pulse of admixture between all inputted populations, but we can change this later if we wish.
We'll need to decide some other parameters: in particular, a recombination rate (presumed constant at the moment), a sequence length, and the number of generations ago at which the admixture occurred.

	>>> rho = 0.1
	>>> length = 10
	>>> gens = 15

This is the minimal input required to run a ``slime`` simulation.

To generate a SLiM script with this information, we'll create a ``slime.RecentHistory`` object containing all of these parameters as input:

	>>> script = slime.RecentHistory(final_gen=gens, chrom_length=length,
            reference_configs=[ref0_config, ref1_config], adm_configs=adm_config,
            prop=adm_props)
