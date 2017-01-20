Forward-time Simulation
=======================

Forward time simulation for discrete generations is provided through the :class:`Population` class.
Populations can undergo a generation of random mating with :meth:`Population.advance_generation`. 
The size of the population and increase or decrease as necessary.
For each individual in the next generation, a new individual is formed by selecting simulated gametes from a randomly selected mother and father.
Initial individuals are added to the population with :meth:`Population.founder_individual`. 

When advancing generations the individuals from the previous generation are removed, but the objects themselves are not deleted.
This lets you track the genealogies of individuals back to the initial population. 
If you are simulating large datasets (individuals or variants) it may be necessary to clear the genotypes of the previous generation, or to remove the old individuals entirely.


Simulating Phenotypes
---------------------

Pydigree can simulate phenotypes under a quantitative genetic models.
Genetic effects on traits are specified by :class:`QuantitativeTraitArchitecture` objects, which are collections of :class:`QuantitativeGeneticEffect` objects.
Phenotypes are realized for individuals with the :meth:`Individual.predict_phenotype` method.

Traits can fixed to specific parameters. 
Narrow-sense heritability (:math:`h^2`) can be specified, and a corresponding amound of random noise added when traits are realized. 
Traits can also be rescaled to have a user-specified mean and variance.

Each :class:`QuantitativeGeneticEffect` instance contains a reference to a locus, as well as two parameters describing the effect, :attr:`a` and :attr:`k`. 
The parameter :attr:`a` describes the additive effect at a locus, *i.e.* each copy of a minor allele will add :attr:`a` to the breeding value.
:attr:`k` sets the deviance from additivity (*i.e.* dominance effect). 
Default value is 0, so that the effect of the genotype Aa is half that of aa. 
When :attr:`k` is not 0, the effect of the heterozygote is :math:`(1+k)a`. 

::

    effect:      0                        (1+k)a            2a 
                 |---------------------------|--------------|
    genotype:  A1/A1                       A1/A2          A2/A2 

When :attr:`k` = 1, the allele effect is purely dominant, and when :attr:`k` = 0, it is purely recessive.

Specifying the effects this way has a couple of advantages. 
First, computing individual breeding values is straightforward.
Second, since other quantitative genetic effects rely on allele frequencies, this allows trait variance to vary over the course of a simulation (if :math:`h^2` is not fixed).


Dichotomous Traits
^^^^^^^^^^^^^^^^^^
Quantitative traits can be transformed into binary traits by applying a threshold to the trait. 
When the architecture realizes the trait for the :class:`Individual` it will be set to a bool indicating that the threshold has been passed.

Fully-penetrant Mendelian traits can be considered as special quantiative traits that only involve a single locus. 
They can be modeled as a :class:`QuantativeGeneticEffect` that adds some breeding value, and any individual with a non-zero breeding value is affected.  

Chromosome Pools
----------------

