# Interactions and impacts of keystone species   
**J. J. Borrelli**

## Introduction

  

- Talk about the microbiome and microbial community ecology  

In the human gut there is a thriving microbial community with as many as 500 species coexisting. Recent advances in metagenomic sequencing have allowed us to catalogue these species and describe the variation in community structure and composition across individuals and within individuals over time. Several factors likely impact the composition of the community. Variation in diet has been shown to influence microbial community composition. This is likely because different microbes are able to more efficiently use different resources (macronutrients). Addtionally, either by hitchhiking on food or other means of transfer, new microbes may invade the community. For example, _Eschericha coli_ may be introduced to the gut via ingestion of contaminated or undercooked foods. The host enviroment can also influence which species are able to coexist in the gut community either through an immune response or mediated through some kind of niche-selection. Finally, the interactions among microbes may both set the boundaries for community composition and drive the response of the community to the other external impacts. 

Microbes can compete with one another directly for limited resources in the gut (both food and space). Alternatively they may interact via the production of metabolites. These compounds can be either beneficial or detrimental to the growth of other microbial populations. Thus, in microbial communities we may expect to see all five major interaction types: competition (-,-), mutualism (+,+), parasitism/predation (+,-), amensalism (-,0), and commensalism (+,0).  

An understanding of how microbes interact may allow us to 

- Talk about microbial community stability  



- Talk about identifying keystone species in microbial communities  

Keystone species are those that have a disproportional impact on the communtiy relative to thier abundance. The textbook case of a keystone species is the sea star _Pisaster_ in tidal communities of the Pacific Northwest. _Pisaster_ feeds on highly competitive mussels, clearing precious real-estate in the inter-tidal zone for multiple sedentary species to colonize. Removing _Pisaster_ from these communities leads to the competitive exclusion of most of these species, and the domination of the community by a single mussel. Other examples of keystone species include otters, whose predatory interaction with sea urchins helps stabilize the kelp forest, and [daphnia?].

Most studies involving keystone species are focused on a single interaction type, predator-prey. For example, Paine (1980) highlights that _Pisaster_ is effective at reducing competition among mussels because it is a generalist that consumes prey at a range of sizes. The impact that different species may have on their community, however, may be mediated by different types of interaction including competition and mutualism. Moreover, keystone interactions of different types may be variable in the form of the impact to the community. 

Typically our understanding of keystone species' impacts have come from species removal experiments. In this paper we describe an _in silico_ species removal experiment. We determine the impact that each species has on a simulated community  


## Methods

Simulations began by generating a species pool whose interactions were known. Interactions among species were assigned using the Erdos-Renyi model for random networks with 200 nodes and a connectance of 0.2. Each species was assigned a self-damping term drawn from a beta distributions (a = 1.1, b = 5) scaled between -5 and 0. Off-diagonal interaction strengths were drawn from a normal distribution based on the interactions derived by Stein et al. (2013). This gave a normal distribution with a mean of -0.07 and standard deviation of 0.354. Growth rates for all species were positive and drawn from a uniform distribution between 0.1 and 1. 

Each individual community was created by sampling 50 species from the pool. All interactions defined above were assumed to be universal, so individual communities represented subsets of the initial species pool. Dynamics of each community were simulated using the Lotka-Volterra model of species interactions,

dN/dt = rN + \sum(a_{ij}*N_j)  eq(1)   

where _r_ is the species specific growth rate, _aij_ is the effect of species _j_ on species _i_. The simulation was run for 1000 time steps, which was long enough for most communities to reach equilibrium (constant abundances) or steady-state (stable attractor). Note the difference is that at steady state the abundances of populations may be fluctuating, but they will remain that way unless perturbed. 

In order to identify which species were important to the stability of the community we systematically removed each species (one at a time) and measured the changes in the community. The starting point for species removals were the equilibrial/steady-state communities. Following the removal of a species, the resulting community dynamics were simulated for 1000 time steps using the same model and parameters as the initial community. At the end of each simulation, the impact on the community was measured using five metrics: (1) the mean change in abundance, (2) mean coefficient of variation in the first 50 time steps following removal, (3) mean coefficient of variation in the last 200 time steps following removal, (4) persistence (the fraction of species with positive abundance), and (5) eigenvalue of the resulting Jacobian matrix with the largest real part.

The types and strengths of the interactions each species participated in were identified for every community. 

## Results 

Equilibrium local communities ranged from 16 to 34 species (median = 24), with connectances between 0.18 and 0.3 (median = 0.23).     


## Discussion



