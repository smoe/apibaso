#!/usr/bin/make

.SUFFIXES: .fasta .spexs .hmmtop .spexs_patterns .spexs_output .spexs_incidence_random .spexs_pattern_similarity .tsv .pattern_similarity .pattern_incidence_on_random


#
# The following folders are created:
#
#  comparisons of apical and basolateral:
#
#       Sequences as a whole in both directions
#
#	Only
#		inner residues
#		outer residues
#		transmembranous residues
#		
#  distributions of patterns
#

target: all

NUMPATTERNS=80
LITERATURE=01_literature_protein_data
SEQUENCES=02_retrieved_sequences
SUBSEQUENCES=03_derived_protein_subsequences
CONSERVED=04_conserved_sequences
COMPARISONS=05_comparisons
DISTRIBUTIONS=06_distributions
ANALYSES=07_analyses

APIBASOFILTERPL=src/apibaso_filter.pl
FASTA2SPEXS=./bin/fasta2specs
APIBASOXML=$(LITERATURE)/apibaso.xml
PATTERNDISTRIBUTIONPL=src/patterndistribution.pl
HMMTOPSEQUENCESPL=src/hmmtop2sequences.pl
DATATRANSFORMING=src/data_transformation_for_SVMs_of_Kai.R

DOMAIN=in

%.spexs: %.fasta $(FASTA2SPEXS)
	$(FASTA2SPEXS) $< > $@

%.spexs_patterns: %.spexs_output
	sort +1  $< |head -$(NUMPATTENRS) | cut -f1 -d\  > $@

%.spexs_incidence_random: %.spexs_output
	src/evaluate_pattern.pl $< > $@

#%.spexs_pattern_similarity: %.spexs_incidence_random src/kai_alles_gegen_alles/doit
#	src/kai_alles_gegen_alles/doit $< $@

%.pattern_incidence_on_random: %.tsv
	src/evaluate_pattern.pl $< > $@

%.pattern_similarity: %.pattern_incidence_on_random src/kai_alles_gegen_alles/doit
	src/kai_alles_gegen_alles/doit $< $@

src/kai_alles_gegen_alles/doit: src/kai_alles_gegen_alles/doit.c
		gcc -o $@ $< -Isrc/kai_alles_gegen_alles/ -lm src/kai_alles_gegen_alles/utils.c src/kai_alles_gegen_alles/DiskUtils.c -D__USE_SSE__ -msse2

# The HMMTOP results require FASTA files as input
%.hmmtop: %.fasta
	hmmtop -if=$< -of=$@

# Shuffling of sequences with respective EMBOSS tool
$(SUBSEQUENCES)/apical.shuffled.complete.fasta: $(SUBSEQUENCES)/apical.complete.fasta
	shuffleseq $< $@

$(SUBSEQUENCES)/basolateral.shuffled.complete.fasta: $(SUBSEQUENCES)/basolateral.complete.fasta
	shuffleseq $< $@

# Aim: lists of patterns as derived from spexs
all:	
	@echo "Possible targets are:"
	@echo " retrieval      - fetches reference sequences"
	@echo " orthologues    - fetches reference sequences and their orthologues"
	@echo " subsequences   - separates external, internal and tm domains"
	@echo " hmmtop"
	@echo " comparisons"
	@echo " distributions"
	@echo " analyses"

#
# Retrieval of sequences across species ... 
#

retrieval: retrieval_stamp

$(SEQUENCES)/reference/apical.complete.fasta $(SEQUENCES)/reference/basolateral.complete.fasta: $(LITERATURE)/apibaso.xml
	echo Fetching sequences as specified in apibaso.xml
	if [ ! -d $(SEQUENCES)/reference ] ; then mkdir -p $(SEQUENCES)/reference; fi
	src/apibaso_filter.pl --xmlfile $(APIBASOXML) --destdir $(SEQUENCES)/reference

retrieval_stamp: $(SEQUENCES)/reference/apical.complete.fasta \
		 $(SEQUENCES)/reference/basolateral.complete.fasta
	touch $@


orthologues: orthologues_stamp

orthologues_stamp: $(LITERATURE)/apibaso.xml
	echo Fetching sequences as specified in apibaso.xml
	echo        !! and their orthologues !!
	src/apibaso.pl --xmlfile $(APIBASOXML) --destdir $(SEQUENCES)
	touch $@

subsequences: \
	$(SUBSEQUENCES)/apical.tm.fasta \
 	$(SUBSEQUENCES)/apical.tmin.fasta \
 	$(SUBSEQUENCES)/apical.in.fasta \
 	$(SUBSEQUENCES)/apical.out.fasta \
 	$(SUBSEQUENCES)/basolateral.tm.fasta \
 	$(SUBSEQUENCES)/basolateral.tmin.fasta \
 	$(SUBSEQUENCES)/basolateral.in.fasta \
 	$(SUBSEQUENCES)/basolateral.out.fasta

	$(SUBSEQUENCES)/apical.tm.fasta \
	$(SUBSEQUENCES)/apical.tmin.fasta \
	$(SUBSEQUENCES)/apical.in.fasta \
	$(SUBSEQUENCES)/apical.out.fasta \
	$(SUBSEQUENCES)/apical.tm.shuffled \
	$(SUBSEQUENCES)/apical.tmin.shuffled \
	$(SUBSEQUENCES)/apical.in.shuffled \

$(SUBSEQUENCES)/apical.out.shuffled: \
	$(SUBSEQUENCES)/apical.complete.hmmtop $(SUBSEQUENCES)/apical.shuffled.complete.fasta $(HMMTOPSEQUENCESPL)
	perl $(HMMTOPSEQUENCESPL) --hmmtop $(SUBSEQUENCES)/apical.complete.hmmtop --sequences $(SUBSEQUENCES)/apical.complete.fasta \
	    --corename=$(SUBSEQUENCES)/apical
	perl $(HMMTOPSEQUENCESPL) --hmmtop $(SUBSEQUENCES)/apical.complete.hmmtop --sequences $(SUBSEQUENCES)/apical.shuffled.complete.fasta \
	    --corename=$(SUBSEQUENCES)/apical.shuffled

	$(SUBSEQUENCES)/basolateral.tm.fasta \
	$(SUBSEQUENCES)/basolateral.tmin.fasta \
	$(SUBSEQUENCES)/basolateral.in.fasta \
	$(SUBSEQUENCES)/basolateral.out.fasta \
	$(SUBSEQUENCES)/basolateral.tm.shuffled \
	$(SUBSEQUENCES)/basolateral.tmin.shuffled \
	$(SUBSEQUENCES)/basolateral.in.shuffled \
	$(SUBSEQUENCES)/basolateral.out.shuffled: \
 	$(SUBSEQUENCES)/basolateral.complete.hmmtop $(SUBSEQUENCES)/basolateral.shuffled.complete.fasta $(HMMTOPSEQUENCESPL)
	perl $(HMMTOPSEQUENCESPL) --hmmtop $(SUBSEQUENCES)/basolateral.complete.hmmtop --sequences $(SUBSEQUENCES)/basolateral.complete.fasta \
	    --corename=$(SUBSEQUENCES)/basolateral
	perl $(HMMTOPSEQUENCESPL) --hmmtop $(SUBSEQUENCES)/basolateral.complete.hmmtop --sequences $(SUBSEQUENCES)/basolateral.shuffled.complete.fasta\
	    --corename=$(SUBSEQUENCES)/basolateral.shuffled

	$(SUBSEQUENCES)/apical.complete.fasta $(SUBSEQUENCES)/basolateral.complete.fasta: $(APIBASOXML) $(APIBASOFILTERPL)
	@[ -d $(SUBSEQUENCES) ] || mkdir $(SUBSEQUENCES)
	$(APIBASO@[ -d $(SUBSEQUENCES) ] || mkdir $(SUBSEQUENCES)FILTERPL) --xmlfile $(APIBASOXML)
	mv apical.complete.fasta basolateral.complete.fasta $(SUBSEQUENCES)


#
# hmmtopology for each .fasta entry in orthologues executed
#

hmmtop: hmmtop.stamp

hmmtop.stamp:  orthologues_stamp src/tmpredict.sh
	echo Topology prediction for each .fasta file with HMMTOP
	echo   apical
	src/tmpredict.sh $(SEQUENCES)/apical $(SEQUENCES)/apical_hmmtopology
	echo   basolateral
	src/tmpredict.sh $(SEQUENCES)/basolateral $(SEQUENCES)/basolateral_hmmtopology
	echo Completed topology prediction
	touch hmmtop.stamp

#
# take all internal sequences (.fasta format) and match them with
# orthologues (output .aln format)
#

align: 
	$(MAKE) DOMAIN=in align_impl
	$(MAKE) DOMAIN=out align_impl
	$(MAKE) DOMAIN=tm align_impl
	$(MAKE) DOMAIN=complete align_impl

align_impl align_$(DOMAIN).stamp: hmmtop.stamp src/align_tcoffe.sh
	echo Alignment with t_coffee of all .$(DOMAIN).fasta files
	echo    apical protein sequences
	src/align_tcoffe.sh $(SEQUENCES)/apical_hmmtopology $(SUBSEQUENCES)/apical_$(DOMAIN)_conserved $(DOMAIN)
	echo    basolateral protein sequences
	src/align_tcoffe.sh $(SEQUENCES)/basolateral_hmmtopology $(SUBSEQUENCES)/basolateral_$(DOMAIN)_conserved $(DOMAIN)
	echo Alignment completed
	touch align_$(DOMAIN).stamp


#
# search for conserved sequence string in all intern subsequences
# and save them as .fasta files in CONSERVED folder
#
#

search: 
	$(MAKE) DOMAIN=in search_impl
	$(MAKE) DOMAIN=out search_impl
	$(MAKE) DOMAIN=tm search_impl
	$(MAKE) DOMAIN=complete search_impl
	
search_impl search_$(DOMAIN).stamp: align_$(DOMAIN).stamp src/search.sh
	echo Starting: Fetching .aln files from alignment
	@[ -d $(CONSERVED) ] || mkdir $(CONSERVED)
	src/search.sh $(SUBSEQUENCES)/apical_$(DOMAIN)_conserved $(CONSERVED)/apical_$(DOMAIN)_conserved
	src/search.sh $(SUBSEQUENCES)/basolateral_$(DOMAIN)_conserved $(CONSERVED)/basolateral_$(DOMAIN)_conserved
	echo Completed: Fetching .aln files from alignment
	touch search_$(DOMAIN).stamp

mapping:
	$(MAKE) DOMAIN=in mapping_impl
	$(MAKE) DOMAIN=out mapping_impl
	$(MAKE) DOMAIN=tm mapping_impl
	$(MAKE) DOMAIN=complete mapping_impl
	
mapping_impl mapping_$(DOMAIN).stamp: src/permute.pl search_$(DOMAIN).stamp
	
	echo Mapping of .fasta files created by search

	@[ -d $(CONSERVED)/apical_$(DOMAIN)_conserved/main ] || mkdir $(CONSERVED)/apical_$(DOMAIN)_conserved/main
	@[ -d $(CONSERVED)/apical_$(DOMAIN)_conserved/test ] || mkdir $(CONSERVED)/apical_$(DOMAIN)_conserved/test
	src/permute.pl --inputfolder $(CONSERVED)/apical_$(DOMAIN)_conserved
	
	@[ -d $(CONSERVED)/basolateral_$(DOMAIN)_conserved/main ] || mkdir $(CONSERVED)/basolateral_$(DOMAIN)_conserved/main
	@[ -d $(CONSERVED)/basolateral_$(DOMAIN)_conserved/test ] || mkdir $(CONSERVED)/basolateral_$(DOMAIN)_conserved/test	
	src/permute.pl --inputfolder $(CONSERVED)/basolateral_$(DOMAIN)_conserved
	
	echo mapping job for basolateral as well as apical completed 
	touch mapping_$(DOMAIN).stamp

#
# take .fasta files and convert them into spexs files
#
#

conversion_to_spexs: 
	$(MAKE) DOMAIN=in conversion_impl
	$(MAKE) DOMAIN=out conversion_impl
	$(MAKE) DOMAIN=tm conversion_impl
	$(MAKE) DOMAIN=complete conversion_impl

conversion_impl conversion_$(DOMAIN).stamp: \
	    mapping_$(DOMAIN).stamp \
	    $(CONSERVED)/apical_$(DOMAIN)_conserved_main.spexs \
	    $(CONSERVED)/basolateral_$(DOMAIN)_conserved_main.spexs
	touch conversion_$(DOMAIN)_main.stamp

$(CONSERVED)/basolateral_$(DOMAIN)_conserved_main.fasta: \
	    $(wildcard $(CONSERVED)/basolateral_$(DOMAIN)_conserved/main/*.fasta)
	cat $(CONSERVED)/basolateral_$(DOMAIN)_conserved/main/*.fasta > $(CONSERVED)/basolateral_$(DOMAIN)_conserved_main.fasta
	cat $(CONSERVED)/basolateral_$(DOMAIN)_conserved/test/*.fasta > $(CONSERVED)/basolateral_$(DOMAIN)_conserved_test.fasta

$(CONSERVED)/apical_$(DOMAIN)_conserved_main.fasta: \
	    $(wildcard $(CONSERVED)/apical_$(DOMAIN)_conserved/main/*.fasta)
	cat $(CONSERVED)/apical_$(DOMAIN)_conserved/main/*.fasta > $(CONSERVED)/apical_$(DOMAIN)_conserved_main.fasta
	cat $(CONSERVED)/apical_$(DOMAIN)_conserved/test/*.fasta > $(CONSERVED)/apical_$(DOMAIN)_conserved_test.fasta



#
#comparison of .spexs files containing the internal sequences 
#



compare_internal: compare_internal_1 compare_internal_2

compare_internal_1: compare_internal_1.stamp
# 
# # 	$(MAKE) DOMAIN=in compare_impl
# 	$(MAKE) DOMAIN=out compare_impl
# 	$(MAKE) DOMAIN=tm compare_impl
# # 	$(MAKE) DOMAIN=complete compare_impl

compare_internal_1.stamp: \
	    $(COMPARISONS)/basolateral_in.vs.all_conserved_main.spexs_output \
	    $(COMPARISONS)/apical_out.vs.all_conserved_main.spexs_output \
	    $(COMPARISONS)/apical_tm.vs.all_conserved_main.spexs_output

$(COMPARISONS)/basolateral_in.vs.all_conserved_main.spexs_output: \
	    $(CONSERVED)/basolateral_in_conserved_main.spexs \
	    $(CONSERVED)/apical_in_conserved_main.spexs
	@[ -d $(COMPARISONS) ] || mkdir $(COMPARISONS)
	$(SPEXS) $(SPEXS_OPTIONS) \
	-f $(CONSERVED)/basolateral_in_conserved_main.spexs \
	-f $(CONSERVED)/basolateral_out_conserved_main.spexs \
	-f $(CONSERVED)/basolateral_tm_conserved_main.spexs \
	-f $(CONSERVED)/apical_complete_conserved_main.spexs \
	> $(COMPARISONS)/basolateral_in.vs.all_conserved_main.spexs_output

	
$(COMPARISONS)/apical_out.vs.all_conserved_main.spexs_output: \
	    $(CONSERVED)/basolateral_in_conserved_main.spexs \
	    $(CONSERVED)/apical_in_conserved_main.spexs
	@[ -d $(COMPARISONS) ] || mkdir $(COMPARISONS)
	$(SPEXS) $(SPEXS_OPTIONS) \
	-f $(CONSERVED)/apical_out_conserved_main.spexs \
	-f $(CONSERVED)/apical_in_conserved_main.spexs \
	-f $(CONSERVED)/basolateral_complete_conserved_main.spexs \
	> $(COMPARISONS)/apical_out.vs.all_conserved_main.spexs_output

$(COMPARISONS)/apical_tm.vs.all_conserved_main.spexs_output: $(CONSERVED)/basolateral_in_conserved_main.spexs $(CONSERVED)/apical_in_conserved_main.spexs
	@[ -d $(COMPARISONS) ] || mkdir $(COMPARISONS)
	$(SPEXS) $(SPEXS_OPTIONS) \
	-f $(CONSERVED)/apical_tm_conserved_main.spexs \
	-f $(CONSERVED)/apical_in_conserved_main.spexs \
	-f $(CONSERVED)/basolateral_complete_conserved_main.spexs \
	> $(COMPARISONS)/apical_tm.vs.all_conserved_main.spexs_output
	touch compare_internal_1.stamp



compare_internal_2:
	$(MAKE) DOMAIN=in compare_impl
	$(MAKE) DOMAIN=out compare_impl
	$(MAKE) DOMAIN=tm compare_impl
	$(MAKE) DOMAIN=complete compare_impl

compare_impl compare_$(DOMAIN).stamp: $(COMPARISONS)/basolateral_$(DOMAIN).vs.apical_$(DOMAIN)_conserved_main.spexs_output \
				$(COMPARISONS)/apical_$(DOMAIN).vs.basolateral_$(DOMAIN)_conserved.spexs_output
	touch compare_$(DOMAIN).stamp

$(COMPARISONS)/basolateral_$(DOMAIN).vs.apical_$(DOMAIN)_conserved_main.spexs_output: $(CONSERVED)/basolateral_$(DOMAIN)_conserved_main.spexs $(CONSERVED)/apical_$(DOMAIN)_conserved_main.spexs
	@[ -d $(COMPARISONS) ] || mkdir $(COMPARISONS)
	$(SPEXS) $(SPEXS_OPTIONS) \
	-f $(CONSERVED)/basolateral_$(DOMAIN)_conserved_main.spexs \
	-f $(CONSERVED)/apical_$(DOMAIN)_conserved_main.spexs \
	> $(COMPARISONS)/basolateral_$(DOMAIN).vs.apical_$(DOMAIN)_conserved_main.spexs_output

$(COMPARISONS)/apical_$(DOMAIN).vs.basolateral_$(DOMAIN)_conserved.spexs_output: $(CONSERVED)/basolateral_$(DOMAIN)_conserved_main.spexs $(CONSERVED)/apical_$(DOMAIN)_conserved_main.spexs
	@[ -d $(COMPARISONS) ] || mkdir $(COMPARISONS)
	$(SPEXS) $(SPEXS_OPTIONS) \
	-f $(CONSERVED)/apical_$(DOMAIN)_conserved_main.spexs \
	-f $(CONSERVED)/basolateral_$(DOMAIN)_conserved_main.spexs \
	> $(COMPARISONS)/apical_$(DOMAIN).vs.basolateral_$(DOMAIN)_conserved_main.spexs_output




# $(COMPARISONS)/basolateral_$(DOMAIN).vs.apical_$(DOMAIN)_conserved_main.spexs_output: $(CONSERVED)/basolateral_$(DOMAIN)_conserved_main.spexs $(CONSERVED)/apical_$(DOMAIN)_conserved_main.spexs
# 	@[ -d $(COMPARISONS) ] || mkdir $(COMPARISONS)
# 	$(SPEXS) $(SPEXS_OPTIONS) \
# 	-f $(CONSERVED)/basolateral_$(DOMAIN)_conserved_main.spexs \
# 	-f $(CONSERVED)/apical_$(DOMAIN)_conserved_main.spexs \
# 	> $(COMPARISONS)/basolateral_$(DOMAIN).vs.apical_$(DOMAIN)_conserved_main.spexs_output
# 
# $(COMPARISONS)/apical_$(DOMAIN).vs.basolateral_$(DOMAIN)_conserved.spexs_output: $(CONSERVED)/basolateral_$(DOMAIN)_conserved_main.spexs $(CONSERVED)/apical_$(DOMAIN)_conserved_main.spexs
# 	@[ -d $(COMPARISONS) ] || mkdir $(COMPARISONS)
# 	$(SPEXS) $(SPEXS_OPTIONS) \
# 	-f $(CONSERVED)/apical_$(DOMAIN)_conserved_main.spexs \
# 	-f $(CONSERVED)/basolateral_$(DOMAIN)_conserved_main.spexs \
# 	> $(COMPARISONS)/apical_$(DOMAIN).vs.basolateral_$(DOMAIN)_conserved_main.spexs_output




#
# Describe distribution of patterns in the respective sequences (test and main data set) 
# 
#





distributions:
	$(MAKE) DOMAIN=in distribution_subdomains_impl
	$(MAKE) DOMAIN=out distribution_subdomains_impl
	$(MAKE) DOMAIN=tm distribution_subdomains_impl
	$(MAKE) DOMAIN=complete distribution_subdomains_impl

distribution_subdomains_impl distribution_$(DOMAIN).stamp: compare_internal_1.stamp $(COMPARISONS)/basolateral_in.vs.all_conserved_main.spexs_output $(PATTERNDISTRIBUTIONPL) $(CONSERVED)/apical_$(DOMAIN)_conserved_main.fasta $(CONSERVED)/basolateral_$(DOMAIN)_conserved_main.fasta 

#
# start with specialised patterns from comparison_1
#
#



	@[ -d $(DISTRIBUTIONS) ] || mkdir $(DISTRIBUTIONS)
	for side in apical basolateral; do \
		$(PATTERNDISTRIBUTIONPL) \
		 --patterns $(COMPARISONS)/valid_apical.spexs_output \
		 --sequences $(CONSERVED)/$${side}_$(DOMAIN)_conserved_main.fasta \
		 > $(DISTRIBUTIONS)/valid_apical.patterns_on_$${side}.sequences_$(DOMAIN)_main.distribution.tsv; \
	done
	for side in apical basolateral; do \
		$(PATTERNDISTRIBUTIONPL) \
		 --patterns $(COMPARISONS)/valid_apical.spexs_output \
		 --sequences $(CONSERVED)/$${side}_$(DOMAIN)_conserved_test.fasta \
		 > $(DISTRIBUTIONS)/valid_apical.patterns_on_$${side}.sequences_$(DOMAIN)_test.distribution.tsv; \
	done
	for side in apical basolateral; do \
		$(PATTERNDISTRIBUTIONPL) \
		 --patterns $(COMPARISONS)/valid_basolateral.spexs_output \
		 --sequences $(CONSERVED)/$${side}_$(DOMAIN)_conserved_main.fasta \
		 > $(DISTRIBUTIONS)/valid_basolateral.patterns_on_$${side}.sequences_$(DOMAIN)_main.distribution.tsv; \
	done
	for side in apical basolateral; do \
		$(PATTERNDISTRIBUTIONPL) \
		 --patterns $(COMPARISONS)/valid_basolateral.spexs_output \
		 --sequences $(CONSERVED)/$${side}_$(DOMAIN)_conserved_test.fasta \
		 > $(DISTRIBUTIONS)/valid_basolateral.patterns_on_$${side}.sequences_$(DOMAIN)_test.distribution.tsv; \
	done
# 	for side in apical basolateral; do \
# 		$(PATTERNDISTRIBUTIONPL) \
# 		 --patterns $(COMPARISONS)/apical_tm.vs.all_conserved_main.spexs_output \
# 		 --sequences $(CONSERVED)/$${side}_$(DOMAIN)_conserved_main.fasta \
# 		 > $(DISTRIBUTIONS)/apical_tm.patterns_on_$${side}.sequences_$(DOMAIN)_main.distribution.tsv; \
# 	done
# 	for side in apical basolateral; do \
# 		$(PATTERNDISTRIBUTIONPL) \
# 		 --patterns $(COMPARISONS)/apical_tm.vs.all_conserved_main.spexs_output \
# 		 --sequences $(CONSERVED)/$${side}_$(DOMAIN)_conserved_test.fasta \
# 		 > $(DISTRIBUTIONS)/apical_tm.patterns_on_$${side}.sequences_$(DOMAIN)_test.distribution.tsv; \
# 	done
# 
# #
# # continue with the genereal set from comparison_2
# #
# #
# 
# 
# 	for side in apical basolateral; do \
# 		$(PATTERNDISTRIBUTIONPL) \
# 		 --patterns $(COMPARISONS)/basolateral_$(DOMAIN).vs.apical_$(DOMAIN)_conserved_main.spexs_output \
# 		 --sequences $(CONSERVED)/$${side}_$(DOMAIN)_conserved_main.fasta \
# 		 > $(DISTRIBUTIONS)/basolateral.patterns_on_$${side}.sequences_$(DOMAIN)_main.distribution.tsv; \
# 	done
# 	for side in apical basolateral; do \
# 		$(PATTERNDISTRIBUTIONPL) \
# 		 --patterns $(COMPARISONS)/basolateral_$(DOMAIN).vs.apical_$(DOMAIN)_conserved_main.spexs_output \
# 		 --sequences $(CONSERVED)/$${side}_$(DOMAIN)_conserved_test.fasta \
# 		 > $(DISTRIBUTIONS)/basolateral.patterns_on_$${side}.sequences_$(DOMAIN)_test.distribution.tsv; \
# 	done
# 	for side in apical basolateral; do \
# 		$(PATTERNDISTRIBUTIONPL) \
# 		 --patterns $(COMPARISONS)/apical_$(DOMAIN).vs.basolateral_$(DOMAIN)_conserved_main.spexs_output \
# 		 --sequences $(CONSERVED)/$${side}_$(DOMAIN)_conserved_main.fasta \
# 		 > $(DISTRIBUTIONS)/apical.patterns_on_$${side}.sequences_$(DOMAIN)_main.distribution.tsv; \
# 	done
# 	for side in apical basolateral; do \
# 		$(PATTERNDISTRIBUTIONPL) \
# 		 --patterns $(COMPARISONS)/apical_$(DOMAIN).vs.basolateral_$(DOMAIN)_conserved_main.spexs_output \
# 		 --sequences $(CONSERVED)/$${side}_$(DOMAIN)_conserved_test.fasta \
# 		 > $(DISTRIBUTIONS)/apical.patterns_on_$${side}.sequences_$(DOMAIN)_test.distribution.tsv; \
# 	done
	touch distribution_$(DOMAIN).stamp

#
#
#	distribution of test patterns in the test set
#
#

# distributions_test:
# 	$(MAKE) DOMAIN=in distribution_test_subdomains_impl
# 	$(MAKE) DOMAIN=out distribution_test_subdomains_impl
# 	$(MAKE) DOMAIN=tm distribution_test_subdomains_impl
# 	$(MAKE) DOMAIN=complete distribution_test_subdomains_impl
# 
# distribution_test_subdomains_impl: distribution_test_$(DOMAIN).stamp
# distribution_$(DOMAIN).stamp: compare_$(DOMAIN).stamp $(COMPARISONS)/basolateral_$(DOMAIN).vs.apical_$(DOMAIN)_conserved_test.spexs_output $(PATTERNDISTRIBUTIONPL) $(CONSERVED)/apical_$(DOMAIN)_conserved_test.fasta $(CONSERVED)/basolateral_$(DOMAIN)_conserved_test.fasta
# 
# 	@[ -d $(DISTRIBUTIONS) ] || mkdir $(DISTRIBUTIONS)
# 		$(PATTERNDISTRIBUTIONPL) \
# 		 --patterns $(COMPARISONS)/basolateral_$(DOMAIN).vs.apical_$(DOMAIN)_conserved_test.spexs_output \
# 		 --sequences $(CONSERVED)/basolateral_$(DOMAIN)_conserved_test.fasta \
# 		 > $(DISTRIBUTIONS)/basolateral_test.patterns_$(DOMAIN).distribution.tsv; \
# 	done
# 	for side in apical basolateral; do \
# 		$(PATTERNDISTRIBUTIONPL) \
# 		 --patterns $(COMPARISONS)/apical_$(DOMAIN).vs.basolateral_$(DOMAIN)_conserved_main.spexs_output \
# 		 --sequences $(CONSERVED)/$${side}_$(DOMAIN)_conserved_main.fasta \
# 		 > $(DISTRIBUTIONS)/apical_test.patterns_on_$${side}.sequences.$(DOMAIN).distribution.tsv; \
# 	done
# 	touch distribution_test_$(DOMAIN).stamp



alignment_2: align_time_stamp
align_time_stamp_2: orthologues_stamp 

	echo Make alignment of sequences derived from orthologues
	echo These are from several species, not only the human
	src/alignment.sh $(SEQUENCES)/apical $(SUBSEQUENCES)/apical_aln
	echo apical wird alignt
	src/alignment.sh $(SEQUENCES)/basolateral $(SUBSEQUENCES)/basolateral_aln 
	echo basolateral wird alignt

	touch align_time_stamp
	echo align_time_stamp touched

#
# Harvest of patterns with SPEXS
#



SPEXS_OPTIONS=-cf spexs.charclass -showratio 6 -goodratio 3 -minpos 10 -wildcards 3 -max_gap_nr 8 -maxgrp 10 -maxlevel 4
SPEXS_SUFFIX=cf_spexs.charclass__showratio_6__goodratio_3__minpos_nr__wildcards_3__max_gap_nr_8__maxgrp_10__maxlevel_4

SPEXS=spexs

comparisons: comparison_time_stamp

	# Conversion of alignments to FASTA files (today, only the consensus strings)

conversion: conversion_stamp
conversion_stamp: align_time_stamp
	echo Fetching .aln files from alignment
	
	src/search.sh $(SUBSEQUENCES)/apical_aln $(CONSERVED)/apical_conserved
	src/search.sh $(SUBSEQUENCES)/basolateral_aln $(CONSERVED)/basolateral_conserved

	echo jawohl David das ging klar! Weiter Junge!
	touch conversion_stamp

#
# Conversion of FASTA files to SPEXS input files 
#
#
#


#$(COMPARISONS)/basolateral.vs.apical_conserved.spexs_output:
#			$(COMPARISONS)/basolateral_conserved.spexs $(COMPARISONS)/apical_conserved.spexs	

compare:
	@[ -d $(COMPARISONS) ] || mkdir $(COMPARISONS)
	$(SPEXS) $(SPEXS_OPTIONS) \
	-f $(CONSERVED)/basolateral_conserved.spexs \
	-f $(CONSERVED)/apical_conserved.spexs \
	> $(COMPARISONS)/basolateral.vs.apical_conserved.spexs_output
	
comparison_time_stamp: $(COMPARISONS)/basolateral.vs.apical_conserved.spexs_output
	touch comparison_time_stamp

# comparisons:
# 	$(MAKE) DOMAIN=in comparison_subdomains_impl
# 	$(MAKE) DOMAIN=out comparison_subdomains_impl
# 	$(MAKE) DOMAIN=tm comparison_subdomains_impl
# 	$(MAKE) DOMAIN=tmin comparison_subdomains_impl
# 	$(MAKE) DOMAIN=complete comparison_subdomains_impl

#comparison_subdomains_evaluation: \
# $(COMPARISONS)/apical.vs.basolateral.$(DOMAIN).spexs_pattern_similarity \
# $(COMPARISONS)/basolateral.vs.apical.$(DOMAIN).spexs_pattern_similarity

comparison_subdomains_impl: \
 $(COMPARISONS)/basolateral.vs.apical.$(DOMAIN).spexs_output \
 $(COMPARISONS)/apical.vs.basolateral.$(DOMAIN).spexs_output

$(COMPARISONS)/basolateral.vs.apical.$(DOMAIN).spexs_output: \
 $(SUBSEQUENCES)/basolateral.$(DOMAIN).spexs \
 $(SUBSEQUENCES)/apical.$(DOMAIN).spexs spexs.charclass
	@[ -d $(COMPARISONS) ] || mkdir $(COMPARISONS)
	$(SPEXS) $(SPEXS_OPTIONS) \
	    -f $(SUBSEQUENCES)/basolateral.$(DOMAIN).spexs \
	    -f $(SUBSEQUENCES)/apical.$(DOMAIN).spexs \
	    > $(COMPARISONS)/basolateral.vs.apical.$(DOMAIN).spexs_output

$(COMPARISONS)/apical.vs.basolateral.$(DOMAIN).spexs_output: \
$(SUBSEQUENCES)/apical.$(DOMAIN).spexs $(SUBSEQUENCES)/basolateral.$(DOMAIN).spexs spexs.charclass
	@[ -d $(COMPARISONS) ] || mkdir $(COMPARISONS)
	$(SPEXS) $(SPEXS_OPTIONS) \
	    -f $(SUBSEQUENCES)/apical.$(DOMAIN).spexs \
	    -f $(SUBSEQUENCES)/basolateral.$(DOMAIN).spexs \
	    > $(COMPARISONS)/apical.vs.basolateral.$(DOMAIN).spexs_output

#
# INVESTIGATION OF PATTERNS
#
distributions2:  
	$(MAKE) DOMAIN=in distributions_subdomains_impl
	$(MAKE) DOMAIN=out distributions_subdomains_impl
	$(MAKE) DOMAIN=tm distributions_subdomains_impl
	$(MAKE) DOMAIN=tmin distributions_subdomains_impl
	$(MAKE) DOMAIN=complete distributions_subdomains_impl

distributions_subdomains_impl2: $(PATTERNDISTRIBUTIONPL) \
 		baso_$(DOMAIN)_distributions api_$(DOMAIN)_distributions

#BASO_SPEXS_OUTPUT_FILE=$(COMPARISONS)/basolateral.$(DOMAIN).$(SPEXS_SUFFIX).$(DOMAIN).spexs_output

baso_$(DOMAIN)_distributions: \
 $(DISTRIBUTIONS)/basolateral.patterns_on_basolateral.sequences.$(DOMAIN).distribution.tsv \
 $(DISTRIBUTIONS)/basolateral.patterns_on_apical.sequences.$(DOMAIN).distribution.tsv

$(DISTRIBUTIONS)/basolateral.patterns_on_basolateral.sequences.$(DOMAIN).distribution.tsv \
$(DISTRIBUTIONS)/basolateral.patterns_on_apical.sequences.$(DOMAIN).distribution.tsv \
$(DISTRIBUTIONS)/basolateral.patterns_on_basolateral.sequences.shuffled.$(DOMAIN).distribution.tsv \
$(DISTRIBUTIONS)/basolateral.patterns_on_apical.sequences.shuffled.$(DOMAIN).distribution.tsv: \
 $(SUBSEQUENCES)/basolateral.$(DOMAIN).fasta \
 $(SUBSEQUENCES)/apical.$(DOMAIN).fasta \
 $(COMPARISONS)/basolateral.vs.apical.$(DOMAIN).spexs_output
	@[ -d $(DISTRIBUTIONS) ] || mkdir $(DISTRIBUTIONS)
	for side in apical basolateral; do \
		$(PATTERNDISTRIBUTIONPL) \
		 --patterns $(COMPARISONS)/basolateral.vs.apical.$(DOMAIN).spexs_output \
		 --sequences $(SUBSEQUENCES)/$$side.$(DOMAIN).fasta \
		 > $(DISTRIBUTIONS)/basolateral.patterns_on_$$side.sequences.$(DOMAIN).distribution.tsv; \
		$(PATTERNDISTRIBUTIONPL) \
		 --patterns $(COMPARISONS)/basolateral.vs.apical.$(DOMAIN).spexs_output \
		 --sequences $(SUBSEQUENCES)/$$side.shuffled.$(DOMAIN).fasta \
		 > $(DISTRIBUTIONS)/basolateral.patterns_on_$$side.sequences.shuffled.$(DOMAIN).distribution.tsv; \
	done
	#$(PATTERNDISTRIBUTIONPL) \
	# --patterns $(COMPARISONS)/basolateral.vs.apical.$(DOMAIN).spexs_output \
	# --sequences $(SUBSEQUENCES)/basolateral.$(DOMAIN).fasta \
	# > $(DISTRIBUTIONS)/basolateral.patterns_on_basolateral.sequences_$(DOMAIN).distribution.tsv
	#$(PATTERNDISTRIBUTIONPL) \
	# --patterns $(COMPARISONS)/basolateral.vs.apical.$(DOMAIN).spexs_output \
	# --sequences $(SUBSEQUENCES)/apical.$(DOMAIN).fasta \
	# > $(DISTRIBUTIONS)/basolateral.patterns_on_apical.sequences_$(DOMAIN).distribution.tsv
	
#API_SPEXS_OUTPUT_FILE=$(COMPARISONS)/apical.$(DOMAIN).$(SPEXS_SUFFIX).$(DOMAIN).spexs_output

api_$(DOMAIN)_distributions: \
 $(DISTRIBUTIONS)/apical.patterns_on_apical.sequences.$(DOMAIN).distribution.tsv \
 $(DISTRIBUTIONS)/apical.patterns_on_basolateral.sequences.$(DOMAIN).distribution.tsv

$(DISTRIBUTIONS)/apical.patterns_on_apical.sequences.$(DOMAIN).distribution.tsv \
$(DISTRIBUTIONS)/apical.patterns_on_basolateral.sequences.$(DOMAIN).distribution.tsv \
$(DISTRIBUTIONS)/apical.patterns_on_apical.sequences.shuffled.$(DOMAIN).distribution.tsv \
$(DISTRIBUTIONS)/apical.patterns_on_basolateral.sequences.shuffled.$(DOMAIN).distribution.tsv: \
 $(SUBSEQUENCES)/apical.$(DOMAIN).fasta \
 $(SUBSEQUENCES)/apical.$(DOMAIN).shuffled \
 $(SUBSEQUENCES)/basolateral.$(DOMAIN).fasta \
 $(SUBSEQUENCES)/basolateral.shuffled.$(DOMAIN).fasta \
 $(COMPARISONS)/apical.vs.basolateral.$(DOMAIN).spexs_output
	@[ -d $(DISTRIBUTIONS) ] || mkdir $(DISTRIBUTIONS)
	for side in apical basolateral; do \
		$(PATTERNDISTRIBUTIONPL) \
		 --patterns $(COMPARISONS)/apical.vs.basolateral.$(DOMAIN).spexs_output \
		 --sequences $(SUBSEQUENCES)/$$side.$(DOMAIN).fasta \
		 > $(DISTRIBUTIONS)/apical.patterns_on_$$side.sequences.$(DOMAIN).distribution.tsv; \
		$(PATTERNDISTRIBUTIONPL) \
		 --patterns $(COMPARISONS)/apical.vs.basolateral.$(DOMAIN).spexs_output \
		 --sequences $(SUBSEQUENCES)/$$side.shuffled.$(DOMAIN).fasta \
		 > $(DISTRIBUTIONS)/apical.patterns_on_$$side.sequences.shuffled.$(DOMAIN).distribution.tsv; \
	done
#	$(PATTERNDISTRIBUTIONPL) \
#	 --patterns $(COMPARISONS)/apical.vs.basolateral.$(DOMAIN).spexs_output \
#	 --sequences $(SUBSEQUENCES)/basolateral.$(DOMAIN).fasta \
#	 > $(DISTRIBUTIONS)/apical.patterns_on_basolateral.sequences_$(DOMAIN).distribution.tsv
#	$(PATTERNDISTRIBUTIONPL) \
#	 --patterns $(COMPARISONS)/apical.vs.basolateral.$(DOMAIN).spexs_output \
#	 --sequences $(SUBSEQUENCES)/apical.$(DOMAIN).fasta \
#	> $(DISTRIBUTIONS)/apical.patterns_on_apical.sequences_$(DOMAIN).distribution.tsv

# Analyses
analysis analyses:
	@[ -d $(ANALYSES) ] || mkdir $(ANALYSES)
	$(MAKE) DOMAIN=in analyses_subdomains_evaluation
	$(MAKE) DOMAIN=out analyses_subdomains_evaluation
	$(MAKE) DOMAIN=tm analyses_subdomains_evaluation
	$(MAKE) DOMAIN=tmin analyses_subdomains_evaluation
	$(MAKE) DOMAIN=complete analyses_subdomains_evaluation
	wc -l $(ANALYSES)/*

analyses_subdomains_evaluation: analyses_subdomains_impl
	$(MAKE) $(ANALYSES)/analyses_pattern_on_sequences_$(DOMAIN).pattern_similarity

analyses_subdomains_impl: $(DATATRANSFORMING) \
 $(DISTRIBUTIONS)/basolateral.patterns_on_basolateral.sequences.$(DOMAIN).distribution.tsv \
 $(DISTRIBUTIONS)/basolateral.patterns_on_apical.sequences.$(DOMAIN).distribution.tsv \
 $(DISTRIBUTIONS)/basolateral.patterns_on_basolateral.sequences.shuffled.$(DOMAIN).distribution.tsv \
 $(DISTRIBUTIONS)/basolateral.patterns_on_apical.sequences.shuffled.$(DOMAIN).distribution.tsv \
 $(DISTRIBUTIONS)/apical.patterns_on_apical.sequences.$(DOMAIN).distribution.tsv \
 $(DISTRIBUTIONS)/apical.patterns_on_basolateral.sequences.$(DOMAIN).distribution.tsv \
 $(DISTRIBUTIONS)/apical.patterns_on_apical.sequences.shuffled.$(DOMAIN).distribution.tsv \
 $(DISTRIBUTIONS)/apical.patterns_on_basolateral.sequences.shuffled.$(DOMAIN).distribution.tsv
	@[ -d $(ANALYSES) ] || mkdir $(ANALYSES)
	cat $(DATATRANSFORMING) \
	| sed -e 's/THIS_IS_SUBSTITUTED_BY_THE_REAL_DOMAIN/$(DOMAIN)/' \
	| sed -e 's/THIS_IS_SUBSTITUTED_BY_THE_REAL_ANALYSISDIR/$(ANALYSES)/' \
	| tee /tmp/script.apibaso.R | R --vanilla --slave

#
# Janitor
#

clean:
	find . -name "~*" -o -name "*~" | xargs -r rm

distclean:
	rm -f $(SUBSEQUENCES)/* $(COMPARISONS)/* $(DISTRIBUTIONS)/* $(ANALYSES)/*

backup: clean
	tmpname=~/dontbackup; \
	localdirname=`pwd`; \
	localdirname=`basename $$localdirname`; \
	cd ..; \
	find $$localdirname/src/kai_matlab > $$tmpname; \
	for i in "in" "out" "complete" "tm" "tmin"; do \
		[ -d $$localdirname/src/kai_matlab_$$i ] && find $$localdirname/src/kai_matlab_$$i >> $$tmpname; \
	done; \
	f=~/apibaso_`date  "+%F_%T"|tr -d ":"`.tar.gz; \
	GZIP=-9 tar  -X $$tmpname -czf $$f \
	    $$localdirname/$(APIBASOXML) \
	    $$localdirname/$(COMPARISONS) \
	    $$localdirname/$(DISTRIBUTIONS) \
	    $$localdirname/$(ANALYSES) \
	    $$localdirname/*.charclass \
	    $$localdirname/bin \
	    $$localdirname/README \
	    $$localdirname/COPYRIGHT \
	    $$localdirname/src \
	    $$localdirname/Makefile ; \
	rm -f $$tmpname; \
	wc -c $$f

commit:
	svn commit

.PHONY: clean subsequences comparisons distributions backup all \
	analyses analysis analyses_subdomains_impl \
	api_$(DOMAIN)_distributions baso_$(DOMAIN)_distributions \
	commit retrieval orthologues conversion_to_spexs


.PRECIOUS: $(COMPARISONS)/apical.vs.basolateral.$(DOMAIN).spexs_output $(COMPARISONS)/basolateral.vs.apical.$(DOMAIN).spexs_output $(COMPARISONS)/basolateral.vs.apical.$(DOMAIN).spexs_output 
