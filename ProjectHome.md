# Introduction #

This project curates a dataset of proteins that are distributed differentially to apical and basolateral regions of the [epithelial](http://en.wikipedia.org/wiki/Epithelium) [membrane](http://en.wikipedia.org/wiki/Cell_membrane). The challenge is to collect a sufficiently large number of protein sequences with a known location to allow [data mining](http://en.wikipedia.org/wiki/Data_mining) algorithms to find patterns from which the location of [en.wikipedia.org/wiki/Cell\_biology biologically] yet uncharacterised proteins can be predicted.

This is not easy, since the sorting happens when the protein is already integrated in the membrane, i.e. it is [no longer interpreted a sequence](http://en.wikipedia.org/wiki/Protein_folding). Multiple [molecular interactions](http://en.wikipedia.org/wiki/Protein%E2%80%93protein_interaction) are contributing to the targeting. Even after the the protein was first delivered, by a process called [transcytosis](http://en.wikipedia.org/wiki/Transcytosis). This is also how we get our food from the guts engulfed on one side transported across all the interior to the other side of the cell towards the blood. These transport processes are relevant for various diseases, e.g. for the entering of a virus into a cell or for its assembly. And it also just feels good to learn more about epithelia, alone since [more than 80% of human cancers are of epithelial origin](http://www.cancerhelp.org.uk/about-cancer/what-is-cancer/cells/types-of-cells-and-cancer).

This project is risky, i.e. I do not necessarily expect a success until the 3D structures of membrane proteins are available. This means that I can only invest spare time into this collection, or now and then ask a student to contribute and test another algorithm. I am now asking for some micro-support through [Flattr](http://flattr.com), which may help in various ways. If the Flattr community is helping, it will of course be acknowledged in the paper that will eventually accompany this effort.

# Implementation #

An XML file is prepared to store the essentials of >150 papers that describe the location of at least one epithelial protein. It is available in this google code repository in the folder 01\_literature . This effort yet offers a set of >100 transmembrane proteins that are expressed in the MDCK cell line. The preferred location is specified and if known, also experimental constraints for the signals and their effects are summarised.

# Funding #

There is none, except yours, possibly:

[![](http://api.flattr.com/button/button-compact-static-100x17.png)](http://flattr.com/thing/50250/Trying-to-predict-the-sorting-of-membrane-proteins-in-epithelia)