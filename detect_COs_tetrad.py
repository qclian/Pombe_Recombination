#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#========================================================
#Program	:	detect_COs_tetrad.py
#Contact	:	Qichao Lian [qlian@mpipz.mpg.de]
#Date		:	01.08.2023
#Version	:	1.0
#========================================================


import argparse, hues, datetime, os, sys
from collections import OrderedDict
from scipy.special import comb, perm
import heapq
from decimal import *
import numpy as np
import math
import decimal

startTime = datetime.datetime.now()


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, \
	description="""This script is used to detect COs in tetrad.
Version:	v1.0
Author:		qclian [qlian@mpipz.mpg.de]
""")

Input = parser.add_argument_group("input files")
Input.add_argument('-r',	dest="chr_len",			help='chromosome length files', type=str, required=True)
Input.add_argument('-c', 	dest="centr_reg",		help='centromere regions file',	type=str, required=True)
Input.add_argument('-i', 	dest="tetrad2spores",	help='tetrad-spores file',	type=str, required=True)
Input.add_argument('-m',	dest="marker_list",		help='markers between parents', type=str, required=True)
Input.add_argument('-g',	dest="genotype_folder",	help='genotyped SNP folder (vcf format)', type=str, required=True)

SWargs = parser.add_argument_group("CO & GC detection")
SWargs.add_argument('-mrn',	dest="min_reads_num",	help='minimum reads number (default: 2)', type=int, required=False, default=2)
SWargs.add_argument('-maf',	dest="min_allele_freq",	help='minimum allele frequency (default: 0.8)', type=int, required=False, default=0.8)

Output = parser.add_argument_group("output files")
Output.add_argument('-o',	dest="output_prefix",	help='output file prefix name', type=str, required=True)

args = parser.parse_args()



def main():

	# Input files
	Chromosome_len = ReadChrLen(args.chr_len)
	Centromere = ReadCentroReg(args.centr_reg)

	Tetrad2spore = ReadTetradSporesFile(args.tetrad2spores)
	ParentalMarkers = ReadParentalMarkers(args.marker_list)


	TetradMarkerDepth, TetradMarkerRaw, TetradMarkerInferred, TetradGenoRaw, TetradGenoInferred = \
	TetradGenotypeProfiling(Tetrad2spore, ParentalMarkers, args.genotype_folder, args.min_reads_num, args.min_allele_freq)


	# Identify recombination events
	CandidateCOEvents, CandidateGCEvents = IdentifyRecombination(TetradMarkerInferred)

	# Combine recombination events
	CombinedRecombinationEvents = CombineRecombination(CandidateCOEvents, CandidateGCEvents)

	# Refine recombination events
	RefinedRecombinationEvents = RefineRecombination(ParentalMarkers, TetradGenoInferred, CombinedRecombinationEvents)

	RefinedRecombinationEvents_output = args.output_prefix + ".RefinedRecombinationEvents.txt"
	OutputFile(RefinedRecombinationEvents, RefinedRecombinationEvents_output)



def RefineRecombination(ParentalMarkers, TetradGenoInferred, CombinedRecombinationEvents):
	print()
	hues.info("Refine candidate recombination events")

	RefinedRecombinationEvents = OrderedDict()

	pre_tetrad_id = ''
	pre_chr_id = ''
	for event_key, combinedRecombinationEvent in CombinedRecombinationEvents.items():
		event_tetrad_id = combinedRecombinationEvent[0]
		event_index_tetrad = combinedRecombinationEvent[1]
		event_index_chr = combinedRecombinationEvent[2]
		event_type = combinedRecombinationEvent[3]
		event_chr_id = combinedRecombinationEvent[4]
		event_start = int(combinedRecombinationEvent[5])
		event_stop = int(combinedRecombinationEvent[6])
		event_pos = combinedRecombinationEvent[7]
		event_size = combinedRecombinationEvent[8]
		event_affected_spores = combinedRecombinationEvent[9]
		event_relation = combinedRecombinationEvent[10]

		event_start_refine = event_start
		key = event_chr_id + "_" + str(event_start)
		pre_marker_pos = ''
		for marker_key in ParentalMarkers:
			cur_marker_chr_id, cur_marker_pos = marker_key.split("_")

			marker_inferred_key = event_tetrad_id + "," + cur_marker_chr_id + "," + cur_marker_pos
			if marker_inferred_key not in TetradGenoInferred:
				continue
			
			if marker_key == key and pre_marker_pos != '' and event_type == "GC":
				event_start_refine = pre_marker_pos

			pre_marker_pos = int(cur_marker_pos)

		event_stop_refine = event_stop
		key = event_chr_id + "_" + str(event_stop)
		pre_marker_pos = ''
		for marker_key in ParentalMarkers:
			cur_marker_chr_id, cur_marker_pos = marker_key.split("_")

			marker_inferred_key = event_tetrad_id + "," + cur_marker_chr_id + "," + cur_marker_pos
			if marker_inferred_key not in TetradGenoInferred:
				continue

			if pre_marker_pos == event_stop and cur_marker_chr_id == event_chr_id and event_type == "GC":
				event_stop_refine = cur_marker_pos

			pre_marker_pos = int(cur_marker_pos)


		event_pos_refine = (float(str(event_start_refine)) + float(str(event_stop_refine))) / 2.0
		event_pos_refine = Decimal(event_pos_refine).quantize(Decimal("0.1"), rounding = "ROUND_HALF_UP")
		event_size_refine = float(str(event_stop_refine)) - float(str(event_start_refine)) + 1

		if event_type == "GC":
			event_pos_refine = event_pos

		refinedRecombinationEvent = [event_tetrad_id, event_index_tetrad, event_index_chr, event_type, event_chr_id, \
		event_start_refine, event_stop_refine, event_pos_refine, event_size_refine, event_affected_spores, event_relation]
		RefinedRecombinationEvents[event_key] = refinedRecombinationEvent

	return RefinedRecombinationEvents


def CombineRecombination(CandidateCOEvents, CandidateGCEvents):
	print()
	hues.info("Combine candidate recombination events")

	CombinedRecombinationEvents = OrderedDict()
	events_count = 0
	events_index_tetrad = 0
	events_index_chr = 0
	CandidateGCEvents_tmp = CandidateGCEvents.copy()

	pre_tetrad_id = ''
	pre_chr_id = ''
	CO_overlapped_GC = False
	CO_printed = False

	for co_key, candidateCOEvent in CandidateCOEvents.items():
		co_tetrad_id = candidateCOEvent[0]
		co_chr_id = candidateCOEvent[4]
		co_start = int(candidateCOEvent[5])
		co_stop = int(candidateCOEvent[6])
		co_pos = candidateCOEvent[7]
		co_size = candidateCOEvent[8]
		co_affected_spores = candidateCOEvent[9]

		if CO_printed:
			pass
		else:
			events_count += 1
		
		if pre_tetrad_id == '':
			pre_tetrad_id = co_tetrad_id
			pre_chr_id = co_chr_id
			events_index_tetrad = 1
			events_index_chr = 1
		else:
			if co_tetrad_id == pre_tetrad_id:
				if CO_printed:
					pass
				else:
					events_index_tetrad += 1
				if co_chr_id == pre_chr_id:
					if CO_printed:
						pass
					else:
						events_index_chr += 1
				else:

					if CO_printed:
						pass
					else:
						events_index_chr = 1

					for gc_key, candidateGCEvent in CandidateGCEvents.items():
						gc_tetrad_id = candidateGCEvent[0]
						gc_chr_id = candidateGCEvent[4]
						gc_start = int(candidateGCEvent[5])
						gc_stop = int(candidateGCEvent[6])
						gc_pos = candidateGCEvent[7]
						gc_size = candidateGCEvent[8]
						gc_affected_spores = candidateGCEvent[9]

						if gc_tetrad_id == pre_tetrad_id and gc_chr_id == pre_chr_id and gc_key in CandidateGCEvents_tmp:
							combinedRecombinationEvent = [gc_tetrad_id, events_index_tetrad, events_index_chr, "GC", gc_chr_id, \
							gc_start, gc_stop, gc_pos, gc_size, gc_affected_spores, "NO"]
							CombinedRecombinationEvents[gc_tetrad_id + "," + str(events_count)] = combinedRecombinationEvent
							events_count += 1
							events_index_tetrad += 1
							events_index_chr += 1
						else:
							pass

					events_index_chr = 1
					pre_chr_id = co_chr_id
			else:

				if CO_printed:
					pass
				else:
					events_index_tetrad += 1
					events_index_chr += 1

				for gc_key, candidateGCEvent in CandidateGCEvents.items():
					gc_tetrad_id = candidateGCEvent[0]
					gc_chr_id = candidateGCEvent[4]
					gc_start = int(candidateGCEvent[5])
					gc_stop = int(candidateGCEvent[6])
					gc_pos = candidateGCEvent[7]
					gc_size = candidateGCEvent[8]
					gc_affected_spores = candidateGCEvent[9]

					if gc_tetrad_id == pre_tetrad_id and gc_chr_id == pre_chr_id and gc_key in CandidateGCEvents_tmp:
						combinedRecombinationEvent = [gc_tetrad_id, events_index_tetrad, events_index_chr, "GC", gc_chr_id, \
						gc_start, gc_stop, gc_pos, gc_size, gc_affected_spores, "NO"]
						CombinedRecombinationEvents[gc_tetrad_id + "," + str(events_count)] = combinedRecombinationEvent
						events_count += 1
						events_index_tetrad += 1
						events_index_chr += 1
					else:
						pass

				events_index_tetrad = 1
				events_index_chr = 1
				pre_tetrad_id = co_tetrad_id
				pre_chr_id = co_chr_id

		CO_overlapped_GC = False
		CO_printed = False
		for gc_key, candidateGCEvent in CandidateGCEvents.items():
			gc_tetrad_id = candidateGCEvent[0]
			gc_chr_id = candidateGCEvent[4]
			gc_start = int(candidateGCEvent[5])
			gc_stop = int(candidateGCEvent[6])
			gc_pos = candidateGCEvent[7]
			gc_size = candidateGCEvent[8]
			gc_affected_spores = candidateGCEvent[9]

			if gc_tetrad_id == co_tetrad_id and gc_chr_id == co_chr_id:
				if gc_stop < co_start and gc_key in CandidateGCEvents_tmp:
					combinedRecombinationEvent = [gc_tetrad_id, events_index_tetrad, events_index_chr, "GC", gc_chr_id, \
					gc_start, gc_stop, gc_pos, gc_size, gc_affected_spores, "NO"]
					CombinedRecombinationEvents[gc_tetrad_id + "," + str(events_count)] = combinedRecombinationEvent
					del CandidateGCEvents_tmp[gc_key]
					events_count += 1
					events_index_tetrad += 1
					events_index_chr += 1

				elif gc_start < co_start and gc_stop > co_start and gc_stop < co_stop and gc_key in CandidateGCEvents_tmp:
					CO_overlapped_GC = True

					combinedRecombinationEvent = [gc_tetrad_id, events_index_tetrad, events_index_chr, "GC", \
					gc_chr_id, gc_start, gc_stop, gc_pos, gc_size, gc_affected_spores, "Yes"]
					CombinedRecombinationEvents[gc_tetrad_id + "," + str(events_count)] = combinedRecombinationEvent
					del CandidateGCEvents_tmp[gc_key]
					events_count += 1
					events_index_tetrad += 1
					events_index_chr += 1

					combinedRecombinationEvent = [co_tetrad_id, events_index_tetrad, events_index_chr, "CO", \
					co_chr_id, co_start, co_stop, co_pos, co_size, co_affected_spores, "Yes"]
					CombinedRecombinationEvents[co_tetrad_id + "," + str(events_count)] = combinedRecombinationEvent
					events_count += 1
					events_index_tetrad += 1
					events_index_chr += 1
					CO_printed = True

				elif gc_start > co_start and gc_stop < co_stop and gc_stop > co_start and gc_key in CandidateGCEvents_tmp:
					CO_overlapped_GC = True

					if CO_printed:
						combinedRecombinationEvent = [gc_tetrad_id, events_index_tetrad, events_index_chr, "GC", \
						gc_chr_id, gc_start, gc_stop, gc_pos, gc_size, gc_affected_spores, "Yes"]
						CombinedRecombinationEvents[gc_tetrad_id + "," + str(events_count)] = combinedRecombinationEvent
						del CandidateGCEvents_tmp[gc_key]
						events_count += 1
						events_index_tetrad += 1
						events_index_chr += 1
					else:
						combinedRecombinationEvent = [co_tetrad_id, events_index_tetrad, events_index_chr, "CO", \
						co_chr_id, co_start, co_stop, co_pos, co_size, co_affected_spores, "Yes"]
						CombinedRecombinationEvents[co_tetrad_id + "," + str(events_count)] = combinedRecombinationEvent
						events_count += 1
						events_index_tetrad += 1
						events_index_chr += 1
						CO_printed = True

						combinedRecombinationEvent = [gc_tetrad_id, events_index_tetrad, events_index_chr, "GC", \
						gc_chr_id, gc_start, gc_stop, gc_pos, gc_size, gc_affected_spores, "Yes"]
						CombinedRecombinationEvents[gc_tetrad_id + "," + str(events_count)] = combinedRecombinationEvent
						del CandidateGCEvents_tmp[gc_key]
						events_count += 1
						events_index_tetrad += 1
						events_index_chr += 1

				elif gc_start > co_start and gc_start < co_stop and gc_stop > co_stop and gc_key in CandidateGCEvents_tmp:
					CO_overlapped_GC = True

					combinedRecombinationEvent = [co_tetrad_id, events_index_tetrad, events_index_chr, "CO", \
					co_chr_id, co_start, co_stop, co_pos, co_size, co_affected_spores, "Yes"]
					CombinedRecombinationEvents[co_tetrad_id + "," + str(events_count)] = combinedRecombinationEvent
					events_count += 1
					events_index_tetrad += 1
					events_index_chr += 1
					CO_printed = True

					combinedRecombinationEvent = [gc_tetrad_id, events_index_tetrad, events_index_chr, "GC", \
					gc_chr_id, gc_start, gc_stop, gc_pos, gc_size, gc_affected_spores, "Yes"]
					CombinedRecombinationEvents[gc_tetrad_id + "," + str(events_count)] = combinedRecombinationEvent
					del CandidateGCEvents_tmp[gc_key]
					events_count += 1
					events_index_tetrad += 1
					events_index_chr += 1
				else:
					continue
			else:
				continue

		if CO_overlapped_GC:
			pass
		else:
			combinedRecombinationEvent = [co_tetrad_id, events_index_tetrad, events_index_chr, "CO", co_chr_id, \
			co_start, co_stop, co_pos, co_size, co_affected_spores, "NO"]
			CombinedRecombinationEvents[co_tetrad_id + "," + str(events_count)] = combinedRecombinationEvent

	for gc_key, candidateGCEvent in CandidateGCEvents.items():
		gc_tetrad_id = candidateGCEvent[0]
		gc_chr_id = candidateGCEvent[4]
		gc_start = int(candidateGCEvent[5])
		gc_stop = int(candidateGCEvent[6])
		gc_pos = candidateGCEvent[7]
		gc_size = candidateGCEvent[8]
		gc_affected_spores = candidateGCEvent[9]

		if gc_tetrad_id == pre_tetrad_id and gc_chr_id == pre_chr_id and gc_key in CandidateGCEvents_tmp:
			events_count += 1
			events_index_tetrad += 1
			events_index_chr += 1
			combinedRecombinationEvent = [gc_tetrad_id, events_index_tetrad, events_index_chr, "GC", gc_chr_id, \
			gc_start, gc_stop, gc_pos, gc_size, gc_affected_spores, "NO"]
			CombinedRecombinationEvents[gc_tetrad_id + "," + str(events_count)] = combinedRecombinationEvent
		else:
			pass

	return CombinedRecombinationEvents


def IdentifyRecombination(TetradMarkerInferred):
	print()
	hues.info("Identify candidate recombination events")

	CandidateCOEvents = OrderedDict()
	events_count = 0
	events_index_tetrad = 0
	events_index_chr = 0

	hues.info("Identify candidate recombination events -- COs")
	pre_tetrad_id = ''
	pre_chr_id = ''
	pre_pos = ''
	pre_spores_genotype = ''
	for key, tetradMarkerInferred in TetradMarkerInferred.items():
		cur_tetrad_id, cur_chr_id, cur_pos = key.split(",")
		cur_pos = int(cur_pos)
		cur_spores_genotype = tetradMarkerInferred[3]
		spores_segregation = tetradMarkerInferred[4]

		if spores_segregation == "2:2:0":
			if pre_tetrad_id == '':
				pre_tetrad_id = cur_tetrad_id
				pre_chr_id = cur_chr_id
				pre_pos = cur_pos
				pre_spores_genotype = cur_spores_genotype
				events_index_tetrad = 0
				events_index_chr = 0
			else:
				if cur_tetrad_id == pre_tetrad_id:
					if cur_chr_id == pre_chr_id:
						if cur_spores_genotype == pre_spores_genotype:
							pre_pos = cur_pos
						else:
							affected_spores = CompareSporesGenotype(pre_spores_genotype, cur_spores_genotype)
							co_pos = 0.0
							co_pos = (pre_pos + cur_pos)/2
							co_pos = Decimal(co_pos).quantize(Decimal("0.1"), rounding = "ROUND_HALF_UP")
							co_size = cur_pos - pre_pos + 1

							events_count += 1
							events_index_tetrad += 1
							events_index_chr += 1
							candidateCOEvent = [cur_tetrad_id, events_index_tetrad, events_index_chr, "CO", \
							cur_chr_id, pre_pos, cur_pos, str(co_pos), co_size, affected_spores]
							CandidateCOEvents[cur_tetrad_id + "," + str(events_count)] = candidateCOEvent

							pre_pos = cur_pos
							pre_spores_genotype = cur_spores_genotype
					else:
						pre_chr_id = cur_chr_id
						pre_pos = cur_pos
						pre_spores_genotype = cur_spores_genotype
						events_index_chr = 0
				else:
					pre_tetrad_id = cur_tetrad_id
					pre_chr_id = cur_chr_id
					pre_pos = cur_pos
					pre_spores_genotype = cur_spores_genotype
					events_index_tetrad = 0
					events_index_chr = 0
		elif spores_segregation == "3:1:0":
			continue
		elif spores_segregation == "1:3:0":
			continue
		elif spores_segregation == "4:0:0":
			continue
		elif spores_segregation == "0:4:0":
			continue
		else:
			hues.warn("weird segregation!")
	hues.info(str(events_count) + " COs identified")

	CandidateGCEvents = OrderedDict()
	events_count = 0
	events_index_tetrad = 0
	events_index_chr = 0

	hues.info("Identify candidate recombination events -- GCs")
	pre_tetrad_id = ''
	pre_chr_id = ''
	tract_start = ''
	tract_stop = ''
	pre_spores_genotype = ''
	pre_spores_segregation = ''
	for key, tetradMarkerInferred in TetradMarkerInferred.items():
		cur_tetrad_id, cur_chr_id, cur_pos = key.split(",")
		cur_pos = int(cur_pos)
		cur_spores_genotype = tetradMarkerInferred[3]
		cur_spores_segregation = tetradMarkerInferred[4]

		if cur_spores_segregation == "2:2:0":

			if pre_tetrad_id == '':
				pre_tetrad_id = cur_tetrad_id
				pre_chr_id = cur_chr_id
				tract_start = cur_pos
				tract_stop = cur_pos
				pre_spores_genotype = cur_spores_genotype
				pre_spores_segregation = cur_spores_segregation
				events_index_tetrad = 0
				events_index_chr = 0
			else:
				if cur_tetrad_id == pre_tetrad_id:
					if cur_chr_id == pre_chr_id:
						tract_start = cur_pos
						tract_stop = cur_pos
						pre_spores_genotype = cur_spores_genotype
						pre_spores_segregation = cur_spores_segregation
					else:
						pre_chr_id = cur_chr_id
						tract_start = cur_pos
						tract_stop = cur_pos
						pre_spores_genotype = cur_spores_genotype
						pre_spores_segregation = cur_spores_segregation
						events_index_chr = 0
				else:
					pre_tetrad_id = cur_tetrad_id
					pre_chr_id = cur_chr_id
					tract_start = cur_pos
					tract_stop = cur_pos
					pre_spores_genotype = cur_spores_genotype
					pre_spores_segregation = cur_spores_segregation
					events_index_tetrad = 0
					events_index_chr = 0

		elif cur_spores_segregation == "3:1:0" or cur_spores_segregation == "1:3:0" \
		or cur_spores_segregation == "4:0:0" or cur_spores_segregation == "0:4:0":
			if pre_tetrad_id == '':
				pre_tetrad_id = cur_tetrad_id
				pre_chr_id = cur_chr_id
				tract_start = cur_pos
				tract_stop = cur_pos
				pre_spores_genotype = cur_spores_genotype
				pre_spores_segregation = cur_spores_segregation
				events_index_tetrad = 0
				events_index_chr = 0
			else:
				if cur_tetrad_id == pre_tetrad_id:
					if cur_chr_id == pre_chr_id:
						if cur_spores_genotype == pre_spores_genotype:
							tract_stop = cur_pos

							del CandidateGCEvents[cur_tetrad_id + "," + str(events_count)]

							gc_pos = 0.0
							gc_pos = (tract_start + tract_stop)/2
							gc_pos = Decimal(gc_pos).quantize(Decimal("0.1"), rounding = "ROUND_HALF_UP")
							gc_size = tract_stop - tract_start + 1

							candidateGCEvent = [cur_tetrad_id, events_index_tetrad, events_index_chr, "GC", \
							cur_chr_id, tract_start, tract_stop, str(gc_pos), gc_size, affected_spores]
							CandidateGCEvents[cur_tetrad_id + "," + str(events_count)] = candidateGCEvent

						else:
							tract_start = cur_pos
							tract_stop = cur_pos

							affected_spores = CompareSporesGenotype(pre_spores_genotype, cur_spores_genotype)
							gc_pos = 0.0
							gc_pos = (tract_start + tract_stop)/2
							gc_pos = Decimal(gc_pos).quantize(Decimal("0.1"), rounding = "ROUND_HALF_UP")
							gc_size = tract_stop - tract_start + 1

							events_count += 1
							events_index_tetrad += 1
							events_index_chr += 1
							candidateGCEvent = [cur_tetrad_id, events_index_tetrad, events_index_chr, "GC", \
							cur_chr_id, tract_start, tract_stop, str(gc_pos), gc_size, affected_spores]
							CandidateGCEvents[cur_tetrad_id + "," + str(events_count)] = candidateGCEvent

							pre_spores_genotype = cur_spores_genotype
							pre_spores_segregation = cur_spores_segregation
					else:
						pre_chr_id = cur_chr_id
						tract_start = cur_pos
						tract_stop = cur_pos
						pre_spores_genotype = cur_spores_genotype
						pre_spores_segregation = cur_spores_segregation
						events_index_chr = 0
				else:
					pre_tetrad_id = cur_tetrad_id
					pre_chr_id = cur_chr_id
					tract_start = cur_pos
					tract_stop = cur_pos
					pre_spores_genotype = cur_spores_genotype
					pre_spores_segregation = cur_spores_segregation
					events_index_tetrad = 0
					events_index_chr = 0
		else:
			hues.warn("weird segregation!")
	hues.info(str(events_count) + " GCs identified")

	return CandidateCOEvents, CandidateGCEvents


def CompareSporesGenotype(geno1, geno2):

	geno1_spores = geno1.split(":")
	geno2_spores = geno2.split(":")
	affected_spores = ''

	for x in range(0,4):
		geno1_spore = geno1_spores[x]
		geno2_spore = geno2_spores[x]
		diff = 0
		
		if geno1_spore == geno2_spore:
			diff = 0
		else:
			diff = 1

		if x == 0:
			affected_spores = str(diff)
		else:
			affected_spores = affected_spores + ":" + str(diff)

	return affected_spores 


def OutputFile(output, output_file):
	CheckOutput(output_file)

	OUT = open(output_file, 'w')
	for key, value in output.items():		
		OUT_str = value[0]
		length = len(value)
		for i in range(1, length):
			OUT_str = OUT_str + "\t" + str(value[i])

		OUT.write(OUT_str + "\n")


def TetradGenotypeProfiling(Tetrad2spore, ParentalMarkers, SporeGenotypeFolder, min_reads_num, min_allele_freq):
	print()
	hues.info("Profilling tetrad genotype")

	TetradMarkerDepth = OrderedDict()
	TetradMarkerRaw = OrderedDict()
	TetradMarkerInferred = OrderedDict()
	TetradGenoRaw = OrderedDict()
	TetradGenoInferred = OrderedDict()

	for tetrad_id, spores in Tetrad2spore.items():
		hues.log(tetrad_id + "\t" + str(spores))

		spore_a = spores[0]
		spore_b = spores[1]
		spore_c = spores[2]
		spore_d = spores[3]

		spore_a_file = SporeGenotypeFolder + spore_a + ".geno.vcf"
		spore_b_file = SporeGenotypeFolder + spore_b + ".geno.vcf"
		spore_c_file = SporeGenotypeFolder + spore_c + ".geno.vcf"
		spore_d_file = SporeGenotypeFolder + spore_d + ".geno.vcf"

		spore_a_geno = ReadSporeGenotype(spore_a_file, ParentalMarkers, min_reads_num, min_allele_freq)
		spore_b_geno = ReadSporeGenotype(spore_b_file, ParentalMarkers, min_reads_num, min_allele_freq)
		spore_c_geno = ReadSporeGenotype(spore_c_file, ParentalMarkers, min_reads_num, min_allele_freq)
		spore_d_geno = ReadSporeGenotype(spore_d_file, ParentalMarkers, min_reads_num, min_allele_freq)

		for key, value in ParentalMarkers.items():
			marker_spore_a_geno_raw = spore_a_geno[key]
			marker_spore_b_geno_raw = spore_b_geno[key]
			marker_spore_c_geno_raw = spore_c_geno[key]
			marker_spore_d_geno_raw = spore_d_geno[key]
			chr_id, pos = key.split("_")
			print(key + "\t" + str(value) + 
				"\t" + str(marker_spore_a_geno_raw) + "\t" + str(marker_spore_b_geno_raw) + 
				"\t" + str(marker_spore_c_geno_raw) + "\t" + str(marker_spore_d_geno_raw))
			
			tetradMarkerDepth = [tetrad_id, chr_id, pos, \
			marker_spore_a_geno_raw[0], marker_spore_a_geno_raw[1], marker_spore_a_geno_raw[2], \
			marker_spore_b_geno_raw[0], marker_spore_b_geno_raw[1], marker_spore_b_geno_raw[2], \
			marker_spore_c_geno_raw[0], marker_spore_c_geno_raw[1], marker_spore_c_geno_raw[2], \
			marker_spore_d_geno_raw[0], marker_spore_d_geno_raw[1], marker_spore_d_geno_raw[2]]
			TetradMarkerDepth[tetrad_id + "_" + key] = tetradMarkerDepth
			print(str(tetradMarkerDepth))

			genotype_pattern_raw = marker_spore_a_geno_raw[6] + ":" + marker_spore_b_geno_raw[6] + ":" \
			+ marker_spore_c_geno_raw[6] + ":" + marker_spore_d_geno_raw[6]

			m_marker_num = f_marker_num = na_marker_num = 0
			if marker_spore_a_geno_raw[6] == "m":
				m_marker_num += 1
			elif marker_spore_a_geno_raw[6] == "f":
				f_marker_num += 1
			else:
				na_marker_num += 1

			if marker_spore_b_geno_raw[6] == "m":
				m_marker_num += 1
			elif marker_spore_b_geno_raw[6] == "f":
				f_marker_num += 1
			else:
				na_marker_num += 1

			if marker_spore_c_geno_raw[6] == "m":
				m_marker_num += 1
			elif marker_spore_c_geno_raw[6] == "f":
				f_marker_num += 1
			else:
				na_marker_num += 1

			if marker_spore_d_geno_raw[6] == "m":
				m_marker_num += 1
			elif marker_spore_d_geno_raw[6] == "f":
				f_marker_num += 1
			else:
				na_marker_num += 1
			segregation_pattern_raw = str(m_marker_num) + ":" + str(f_marker_num) + ":" + str(na_marker_num)
			
			tetradMarkerRaw = [tetrad_id, chr_id, pos, genotype_pattern_raw, segregation_pattern_raw]
			TetradMarkerRaw[tetrad_id + "," + chr_id + "," + pos] = tetradMarkerRaw

			tetradGenoRaw = [tetrad_id, chr_id, pos, \
			marker_spore_a_geno_raw[6], marker_spore_b_geno_raw[6], marker_spore_c_geno_raw[6], marker_spore_d_geno_raw[6]]
			TetradGenoRaw[tetrad_id + "," + chr_id + "," + pos] = tetradGenoRaw

			marker_spore_a_geno_infer = marker_spore_a_geno_raw[6]
			marker_spore_b_geno_infer = marker_spore_b_geno_raw[6]
			marker_spore_c_geno_infer = marker_spore_c_geno_raw[6]
			marker_spore_d_geno_infer = marker_spore_d_geno_raw[6]
			segreation_pattern_infer = segregation_pattern_raw
			if segregation_pattern_raw == "2:1:1":
				if marker_spore_a_geno_infer == "NA":
					marker_spore_a_geno_infer = "f"
				elif marker_spore_b_geno_infer == "NA":
					marker_spore_b_geno_infer = "f"
				elif marker_spore_c_geno_infer == "NA":
					marker_spore_c_geno_infer = "f"
				else:
					marker_spore_d_geno_infer = "f"
				segreation_pattern_infer = "2:2:0"
			elif segregation_pattern_raw == "1:2:1":
				if marker_spore_a_geno_infer == "NA":
					marker_spore_a_geno_infer = "m"
				elif marker_spore_b_geno_infer == "NA":
					marker_spore_b_geno_infer = "m"
				elif marker_spore_c_geno_infer == "NA":
					marker_spore_c_geno_infer = "m"
				else:
					marker_spore_d_geno_infer = "m"
				segreation_pattern_infer = "2:2:0"
			elif segregation_pattern_raw == "1:1:2" or segregation_pattern_raw == "0:2:2" or segregation_pattern_raw == "2:0:2" \
			or segregation_pattern_raw == "0:1:3" or segregation_pattern_raw == "1:0:3" or segregation_pattern_raw == "0:0:4":
				continue
			else:
				pass
			genotype_pattern_infer = marker_spore_a_geno_infer + ":" + marker_spore_b_geno_infer + ":" \
			+ marker_spore_c_geno_infer + ":" + marker_spore_d_geno_infer

			tetradMarkerInferred = [tetrad_id, chr_id, pos, genotype_pattern_infer, segreation_pattern_infer]
			TetradMarkerInferred[tetrad_id + "," + chr_id + "," + pos] = tetradMarkerInferred

			tetradGenoInferred = [tetrad_id, chr_id, pos, \
			marker_spore_a_geno_infer, marker_spore_b_geno_infer, marker_spore_c_geno_infer, marker_spore_d_geno_infer]
			TetradGenoInferred[tetrad_id + "," + chr_id + "," + pos] = tetradGenoInferred

	return TetradMarkerDepth, TetradMarkerRaw, TetradMarkerInferred, TetradGenoRaw, TetradGenoInferred


def ReadSporeGenotype(input_file, ParentalMarkers, min_reads_num, min_allele_freq):
	CheckInput(input_file)

	SporeGenoFull = OrderedDict()
	SporeGenoFlt = OrderedDict()
	marker_num_flt = 0

	with open(input_file, 'r') as IN:
		for line in IN:
			if not line.startswith("I"):
				continue

			lines = line.strip().split("\t")

			chr_id = lines[0]
			pos = lines[1]
			key = chr_id + "_" + str(pos)

			infos = lines[9].split(":")
			ref_depth = int(infos[1].split(",")[0])
			alt_depth = int(infos[1].split(",")[1])
			depth = int(infos[2])

			allele_freq = 0.0
			ratio = 0.0
			geno = ""
			if (depth == 0) or (ref_depth == 0 and alt_depth == 0):
				SporeGenoFlt[key] = [ref_depth, alt_depth, depth, 'NA']
			else:
				if ref_depth >= alt_depth:
					ratio = ref_depth / (ref_depth + alt_depth)
					geno = 0
				else:
					ratio = alt_depth / (ref_depth + alt_depth)
					geno = 1
				allele_freq = Decimal(ratio).quantize(Decimal("0.01"), rounding = "ROUND_HALF_UP")

				if ((depth > min_reads_num) or (alt_depth >= min_reads_num)) and (allele_freq > min_allele_freq):
					SporeGenoFlt[key] = [ref_depth, alt_depth, depth, geno]
					marker_num_flt += 1
				else:
					SporeGenoFlt[key] = [ref_depth, alt_depth, depth, 'NA']

	hues.log(input_file + " Input file loaded!")
	hues.log(str(marker_num_flt) + " markers loaded!")

	marker_num_flt = 0
	for key, value in ParentalMarkers.items():
		if key in SporeGenoFlt:
			sporeGenoFull = SporeGenoFlt[key]

			ref_depth = sporeGenoFull[0]
			alt_depth = sporeGenoFull[1]
			depth = sporeGenoFull[2]
			geno = sporeGenoFull[3]
			p1_geno = value[0]
			p2_geno = value[1]

			phase = ''
			if str(geno) == str(p1_geno):
				phase = 'm'
				marker_num_flt += 1
			elif str(geno) == str(p2_geno):
				phase = 'f'
				marker_num_flt += 1
			else:
				print(key + "\t" + str(value) + "\t" + str(sporeGenoFull))
				phase = 'NA'

			sporeGenoFull = [ref_depth, alt_depth, depth, geno, value[0], value[1], phase]
		else:
			sporeGenoFull = [0, 0, 0, 'NA', value[0], value[1], 'NA']
		SporeGenoFull[key] = sporeGenoFull
	hues.log(str(marker_num_flt) + " markers kept!")

	return SporeGenoFull


def ReadTetradSporesFile(input_file):
	CheckInput(input_file)
	Tetrad2spore = OrderedDict()
	tetrad_num = 0

	with open(input_file, 'r') as IN:
		for line in IN:
			lines = line.strip().split("\t")

			if len(lines) < 3:
				hues.error("Input file:\t" + input_file + " format is not correct!")
				sys.exit()
			else:
				spore_name = lines[0]
				spore_id = lines[1]
				tetrad_id = lines[3]

				if tetrad_id in Tetrad2spore:
					spores = Tetrad2spore[tetrad_id]
					spores.append(spore_id)
					Tetrad2spore[tetrad_id] = spores
				else:
					Tetrad2spore[tetrad_id] = [spore_id]
					tetrad_num += 1

	hues.log(input_file + " Input file loaded!")
	hues.log(str(tetrad_num) + " tetrads loaded!")
	return(Tetrad2spore)


def ReadParentalMarkers(input_file):
	CheckInput(input_file)

	ParentalMarkers = OrderedDict()
	marker_num = 0

	with open(input_file, 'r') as IN:
		for line in IN:
			lines = line.strip().split("\t")

			chr_id = lines[0]
			pos = lines[1]
			p1_geno = lines[5]
			p2_geno = lines[6]

			key = chr_id + "_" + str(pos)
			ParentalMarkers[key] = [p1_geno, p2_geno]

			marker_num += 1

	hues.log(input_file + " Input file loaded!")
	hues.log(str(marker_num) + " markers loaded!")
	return ParentalMarkers


def ReadCentroReg(input_centro_reg):
	CheckInput(input_centro_reg)
	centromere = OrderedDict()

	with open(input_centro_reg, 'r') as IN:
		for line in IN:
			if line.startswith("Name"):
				continue

			lines = line.strip().split("\t")

			chr_id = lines[0]
			region = [int(lines[1]), int(lines[2])]

			centromere[chr_id] = region

	hues.log(str(len(centromere)) + " Chromosome centromere regions loaded!")
	return centromere


def ReadChrLen(input_chr_len):
	CheckInput(input_chr_len)
	chr_len = OrderedDict()

	with open(input_chr_len, 'r') as IN:
		for line in IN:
			if line.startswith("Name"):
				continue

			lines = line.strip().split("\t")
			chr_len[lines[0]] = lines[1]

	hues.log(str(len(chr_len)) + " Chromosome length loaded!")
	return chr_len


def CheckInput(input_file):
	if os.path.exists(input_file):
		hues.log("Found input file:\t" + input_file)

		if os.path.getsize(input_file):
			hues.log("Input file:\t" + input_file + " is not empty!")
		else:
			hues.error("Input file:\t" + input_file + " is empty!")
			sys.exit()

	else:
		hues.error("Not found input file:\t" + input_file)
		sys.exit()


def CheckOutput(output_file):
	if os.path.exists(output_file):
		hues.warn("Found output file:\t" + output_file)
		hues.warn("Remove\t" + output_file + "\t...")
		os.remove(output_file)



if __name__ == "__main__":
	main()


stopTime = datetime.datetime.now()
hues.success("This script has run " + str((stopTime - startTime).seconds) + "s!")

