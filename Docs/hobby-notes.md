
__Why do we care about stability?__

Human Health implications

There is a widespread assumption among microbial biologists studying the human gut microbiome that these communities are relatively stable over long periods of time. Evidence supporting this assumption comes primarily from assessments of the microbial community composition at several time points and noting that the same dominant OTUs appear at each time point. The assumption of stability has been translated into the idea that a stable microbial community is important for maintaining the health of the human host. For example, the invasion by the pathogenic bacteria _Clostridium dificile_ disrupts the stable community and can lead to negative health outcomes. Microbiome composition has been linked to diabetes (REFERENCE), obesity (REFERENCE), other things (REFERENCE). 

We want to know what features of the microbial community confer stability. 


__How can we assess stability?__

_Dynamics and Loval Stability_

Need to figure out whether we can use relative abundance measures of OTUs to determine several stability measures. In particular we want to be able to parameterize Lotka-Volterra models of species interactions to assess whether the community is locally stable (in the sense of the eigenvalues) and whether the community has permanence (species maintain positive abundance based on projections from the model). By the approaches that typically require absolute abundances to realtive abundance data we will be able to unlock a tremendous amount of data that has been collected. 

$$  
	Lotka-Volterra Model Equation  
$$  

_Keystone Species_

Discovering the identity of keystone taxa in the gut microbial community is imperative to understanding the factors underlying the stability of the microbiome. Prior studies have used networks of co-occurrence and correlations to assign 'keystoneness' to different OTUs. A keystone species is one that has an disproportional impact on species coexistence relative to its own abundance. We want to use a novel approach, from the field of economics, called market basket analysis (MBA) to identify the relationships among microbial taxa. 

MBA is a data mining technique used to rules of association based on the idea that if a customer buys a group of items then that makes them more or less likely to buy other items. For example, a rule might state that a person buying bananas is also highly likely to purchase apples and strawberries. This type of analysis requires multiple observations of presence and absence information. These observations are typically across individuals, but can also be temporally explicit (i.e., the same person over multiple time points).     

Applying market basket analysis to microbiome data would allow us to determine association rules for different taxa. This analysis could be broadened then, to incorporate abundance information, and could then inform us, for example, if species A is present then species B will be present and abundant, but species C will be present and rare. While MBA is based on co-occurence, it is a better approach because it is not limited to pairwise relationships. Such an approach has not been used before in microbiome research, but has the potential to inform our understanding of how different microbial taxa can impact the whole community. 

A great deal of data that has been collected on the human gut microbiome is ideal for MBA. 