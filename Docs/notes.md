# Notes

## Abundance Interactions

Coyte _et al._ found that microbial communities with competitive interactions were more likely to be stable than those with mutualistic interactions

But microbial communities have abundance distributions that have a few abundant microbes and many rare taxa

__Are there differences in the interactions among abundant and rare taxa?__ 

Based on the data from Stein _et al._ the equilibrium community consists of 5 relative abundant taxa:  
	
	- undefined genus of Enterobacteriaceae
	- Clostridium difficile
	- Blautia
	- undefined genus of unclassified Mollicutes
	- Other
	
and 6 relatively less abundant taxa: 

	- Coprobacillus
	- Akkermansia
	- unclassified Lachnospiraceae
	- Barnesiella
	- undefined genus of Lachnospiraceae
	- Enterococcus  

![communityABUNDANCE](https://raw.githubusercontent.com/jjborrelli/microbial-dyn/master/Images/eqMICROBE1.png =650x)	

Decomposing the interaction matrix (_M_) into two submatrices comprised of only abundant OTUs (_M_{abund, abund}) or only rare OTUs (_M_{rare, rare}) I can assess what the interactions of these two submodules are. In these data there is a difference between the two submodules in their constituent interaction types. 


| Interaction	| _M_ ABUND | _M_ RARE 	|
|:------------:|:---------:|:---------:|
| Competition	| 2	| 6	|
| Mutualism		| 4	| 3	| 
| Predation		| 4	| 6	|

The interactions among the abundant OTUs are characterized primarily by mutualisms (40%) and predation/parasitism (40%), while the interactions among rare OTUs are characterized by competition (40%) and predation (40%). 

Thus the Coyte _et al._ finding that microbial communities should be dominated by competitive interactions may be supplemented by the idea that these competitive interactions may be primarily between taxa that are relatively rare. These competitive and rare taxa may then be what allows the dominant OTUs to be stable despite their mutualistic interactions. 

 