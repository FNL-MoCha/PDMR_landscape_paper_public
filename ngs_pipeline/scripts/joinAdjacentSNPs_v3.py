#!/usr/bin/python3
"""Find adjacent SNPs/MNPs in a VCF, and look in the BAM file to see if these
should be joined into a larger variant.

This script can only handle up to 2 internal compound ALT alleles for every
stretch of adjacent SNPs. This script is not designed to test all possible
combinations of internal joined SNPs for SNP stretches longer than 3.
This script also cannot handle multiallelic variants or indels.

"ALERT" will be added to VCF INFO columns for cases that are too complex.


##########################################################################

 module load samtools
 python joinAdjacentSNPs.py -v 128128/20170910/128128~338-R~L43~WES/calls/128128~338-R~L43~WES.HC_DNASeq.raw.vcf -o test.vcf 1 128128/20170910/128128~338-R~L43~WES/128128~338-R~L43~WES.bwa.final.bam /data/MoCha/patidarr/ref/ucsc.hg19.2bit



"""
from __future__ import division
import copy
import operator
import subprocess
import sys


# We're splitting multiallelics now so should not find any multiallelic VCF lines anymore.
# If the adjacent variants contain some multiallelics, should add alert.
# If the adjacent variants are only multiallelics, should not add an alert.
def checkMultiallelic(variants):
    seenPos = []
    countSplit = 0
    countComma = 0
    for var in variants:
        if not [var[0], var[1]] in seenPos:
            seenPos.append([var[0], var[1]])
        if 'SPLITMULTIALLELIC' in var[7]:
            countSplit += 1
        elif ',' in var[4]:
            countComma += 1
    if countSplit == len(variants) and len(seenPos) == 1:
        # This set of variants came from a single original multiallelic variant. don't add ALERT to this set.
        multiallelic = 'All'
    elif countSplit > 0:
        # Not everything in this set came from a single original multiallelic variant.
        multiallelic = 'Some'
    elif countComma > 0:
        multiallelic = 'Some'
    else:
        multiallelic = 'None'
    return multiallelic


def checkIndels(variants):
    indelsPresent = False
    for var in variants:
        for alt in var[4].split(','):
            if len(var[3]) != len(alt):
                indelsPresent = True
                break
    return indelsPresent


def getReference(chrom, start, end, ref2bit):
    # For twoBit2Fa, coords need to be 0-based, left closed, right open.
    fixStart = int(start) - 1
    fixEnd = int(end) - 1
    coords = "{}:{}-{}".format(chrom, fixStart, fixEnd)
    step_get_ref = subprocess.Popen(
        ['twoBitToFa', '-noMask', ref2bit + ":" + coords, 'stdout'],
        stdout=subprocess.PIPE
    )
    step_remove_header = subprocess.Popen(
        ['grep', '-v', '>'],
        stdin=step_get_ref.stdout,
        stdout=subprocess.PIPE
    )
    step_get_ref.stdout.close()
    # Stuff from twoBitToFa should only ever be a few bases, since we are
    # looking for things that are "adjacent" so there shouldn't be multiple
    # lines in output, but set this up so that multiple lines won't break
    # things.
    ref = "".join(step_remove_header.communicate()[0].rstrip('\n\r').splitlines())
    print("Intervening genomic sequence: %s" % ref, file=sys.stderr)
    return ref


def getJoinedAlleles(variants, ref2bit):
    # Initialize joined alleles
    refJoined = variants[0][3]
    altJoined = variants[0][4]
    for i in range(1, len(variants)):
        # Is there space between this and the previous variant?
        # If there is, add it to the joined alleles.
        currChrom = variants[i][0]
        currPos = int(variants[i][1])
        prevPos = int(variants[i - 1][1])
        prevLength = len(variants[i - 1][3])
        prevEnd = prevPos + prevLength - 1
        if currPos - prevEnd > 1:
            intervening_genomic_seq = getReference(currChrom, prevEnd + 1, currPos, ref2bit)
            refJoined += intervening_genomic_seq
            altJoined += intervening_genomic_seq
        refJoined += variants[i][3]
        altJoined += variants[i][4]
    return refJoined, altJoined


def correctVcfLine(variants, newRef, newAlt):
    if len(variants) > 1:
        correctedLine = copy.copy(variants[0])
        correctedLine[3] = newRef
        correctedLine[4] = newAlt
        correctedLine[7] = correctedLine[7].replace('SNP', 'MNP').replace('snp', 'mnp') + ';JOINED'
    else:
        correctedLine = variants[0]
    return correctedLine


def alertSet(variants):
    alertVars = []
    for var in variants:
        alertLine = copy.copy(var)
        alertLine[7] += ";ALERT"
        alertVars.append("\t".join(map(str, alertLine)))
    return alertVars


def examineAdjacent(adjacentVar, inBam, ref2bit):
    # Minimum number of reads to support an allele
    min_support = 5
    print("Examining adjacent variants:\n%s" % "\n".join(map(str, adjacentVar)), file=sys.stderr)
    correctedVars = []
    if len(adjacentVar) > 3:
        print("WARNING: This is a stretch of %s adjacent variants. This script is not designed to test all possible combinations of internal joined SNPs for stretches longer than 3." % len(adjacentVar), file=sys.stderr)
    if checkIndels(adjacentVar) is True:
        print("WARNING: Variant stretch includes indels.", file=sys.stderr)
        correctedVars = alertSet(adjacentVar)
    if checkMultiallelic(adjacentVar) == 'Some':
        print("WARNING: Variant stretch includes multiallelic variants.", file=sys.stderr)
        correctedVars = alertSet(adjacentVar)
    elif checkMultiallelic(adjacentVar) == 'All':
        print("WARNING: This variant set comes from a single original multiallelic variant and does not represent true adjacent SNPs/MNPs. Leaving variants as is.", file=sys.stderr)
        # correctedVars might have already been set if there are indels, so need to reset.
        correctedVars = []
        for var in adjacentVar:
            correctedVars.append("\t".join(map(str, var)))
    if len(correctedVars) == 0:
        refJoined, altJoined = getJoinedAlleles(adjacentVar, ref2bit)
        print("Joined REF allele is %s" % refJoined, file=sys.stderr)
        print("Joined ALT allele is %s" % altJoined, file=sys.stderr)
        # coordinates of joined region
        chrom = adjacentVar[0][0]
        start = int(adjacentVar[0][1])
        end = int(adjacentVar[-1][1]) + len(adjacentVar[-1][4]) - 1
        coords = chrom + ":" + str(start) + "-" + str(end)

        # get counts of alleles from BAM file.
        observed_alleles = {}
        print("Samtools view -f 2 -F 1792 %s" % coords, file=sys.stderr)
        # don't count duplicates.
        step_samtools = subprocess.Popen(
            ['samtools', 'view', '-f', '2', '-F', '1792', inBam, coords],
            stdout=subprocess.PIPE
        )
        step_cut = subprocess.Popen(['cut', '-f', '4,6,10'], stdin=step_samtools.stdout, stdout=subprocess.PIPE)
        step_samtools.stdout.close()
        bamData = step_cut.communicate()[0].rstrip('\n\r').split('\n')
        if len(bamData[0]) > 0:
            print("Total reads at this position: %s" % len(bamData), file=sys.stderr)
            fullMatchCount = 0
            spanFullLocusCount = 0
            for read in iter(bamData):
                pos, cigar, seq = read.split("\t")
                fullMatchCigar = str(len(seq)) + "M"
                if cigar == fullMatchCigar:
                    fullMatchCount += 1
                    # VCF is 1-based but python is 0-based
                    left = start - int(pos)
                    right = end - int(pos) + 1
                    allele = seq[left:right]
                    if len(allele) == len(altJoined):
                        spanFullLocusCount += 1
                        if allele in observed_alleles:
                            observed_alleles[allele] += 1
                        else:
                            observed_alleles[allele] = 1
            print("Reads with full match CIGAR (%s): %s" % (fullMatchCigar, fullMatchCount), file=sys.stderr)
            print("Reads with full match CIGAR that span the full locus: %s (%.3f%% of all reads)" % (spanFullLocusCount, 100.0 * spanFullLocusCount / len(bamData)), file=sys.stderr)
        # Sort alleles by abundance. This produces a list of duples. For every duple, the first element is the allele
        # and the second is the count.
        observed_alleles_sorted = sorted(observed_alleles.items(), key=operator.itemgetter(1), reverse=True)

        # Require at least 40% usable reads
        if (len(observed_alleles.keys()) >= 1
            and observed_alleles_sorted[0][1] >= min_support
            and 1.0 * spanFullLocusCount / len(bamData) >= 0.4):
            # For now, script handles at most 2 major alleles
            major_alleles = [observed_alleles_sorted[0][0]]
            # major_allele_read_counts = int(observed_alleles_sorted[0][1])
            if (len(observed_alleles_sorted) > 1
                and (observed_alleles_sorted[1][1] >= min_support
                     or observed_alleles_sorted[1][1] >= 0.6 * observed_alleles_sorted[0][1])):
                major_alleles.append(observed_alleles_sorted[1][0])
                # major_allele_read_counts += int(observed_alleles_sorted[1][1])
            # Convert each major allele into a binary array showing wheter it
            # matches REF or ALT at each variant position
            abBinary = []
            for abAllele in major_alleles:
                convBin = []
                for k in range(len(refJoined)):
                    if refJoined[k:k + 1] == altJoined[k:k + 1]:
                        # This is some intervening genomic sequence between the adjacent SNPs.
                        continue
                    elif abAllele[k:k + 1] == refJoined[k:k + 1]:
                        convBin.append(0)
                    else:
                        convBin.append(1)
                abBinary.append(convBin)
            print("Major alleles: %s" % major_alleles, file=sys.stderr)
            print("Major alleles binary arrays: %s" % abBinary, file=sys.stderr)

            # How many reads aren't from the major allele(s)/the reference allele?
            # Only look at variant positions, not intervening genomic sequence
            nonref_minor_allele_read_counts = 0
            for y in range(len(major_alleles), len(observed_alleles_sorted)):
                currAllele = observed_alleles_sorted[y][0]
                nonRef = 0
                for k in range(len(refJoined)):
                    if refJoined[k:k + 1] == altJoined[k:k + 1]:
                        # This is some intervening genomic sequence between the adjacent SNPs.
                        continue
                    elif refJoined[k:k + 1] == currAllele[k:k + 1]:
                        continue
                    else:
                        nonRef += 1
                if nonRef > 0:
                    nonref_minor_allele_read_counts += int(observed_alleles_sorted[y][1])
            nonref_minor_allele_percent = 100 * float(nonref_minor_allele_read_counts) / float(spanFullLocusCount)

            for obs in observed_alleles_sorted:
                print("Found allele %s, count %s" % (obs[0], obs[1]), file=sys.stderr)
            print("Reads from non-reference minor alleles are %.3f%% of total" % nonref_minor_allele_percent, file=sys.stderr)

            if nonref_minor_allele_percent <= 5:
                if ((len(major_alleles) == 2 and altJoined in major_alleles and refJoined in major_alleles)
                        or (len(major_alleles) == 1 and altJoined in major_alleles)
                        or (len(major_alleles) == 2 and [0] * len(adjacentVar) in abBinary and [1] * len(adjacentVar) in abBinary)
                        or (len(major_alleles) == 1 and [1] * len(adjacentVar) in abBinary)):
                    # Join the variants
                    print("Joining variants:\n%s" % "\n".join(map(str, adjacentVar)), file=sys.stderr)
                    joinedVars = correctVcfLine(adjacentVar, refJoined, altJoined)
                    correctedVars.append("\t".join(map(str, joinedVars)))
                else:
                    if len(adjacentVar) > 1:
                        correctedVars = alertSet(adjacentVar)
            else:
                if len(adjacentVar) > 1:
                    correctedVars = alertSet(adjacentVar)
        else:
            if len(adjacentVar) > 1:
                correctedVars = alertSet(adjacentVar)
    return correctedVars


def readVcf(vcfFile):
    vcfData = open(vcfFile).readlines()
    vcfLines = []
    for line in vcfData:
        if line.startswith('#'):
            vcfLines.append(line.strip().split('\t'))
        else:
            vcfLines.append(line.strip().split('\t'))
    return vcfLines


def writeVcf(vcfLines, outputFile):
    with open(outputFile, 'w') as vcfOut:
        for line in vcfLines:
            print("\t".join(map(str, line)), file=vcfOut)


def main():
    # Parsing command-line arguments
    #print (len(sys.argv))
    if len(sys.argv) != 6:
        print("Usage: python joinAdjacentSNPs.py <vcf> <output_vcf> <min_adjacent> <bam> <ref2bit>")
        sys.exit(1)

    vcfFile = sys.argv[1]
    outputFile = sys.argv[2]
    minAdjacent = int(sys.argv[3])
    bamFile = sys.argv[4]
    ref2bit = sys.argv[5]

    vcfLines = readVcf(vcfFile)
    correctedVcfLines = []
    for i in range(len(vcfLines)):
        line = vcfLines[i]
        if line[0].startswith('#'):
            correctedVcfLines.append(line)
        else:
            currChrom = line[0]
            currPos = int(line[1])
            currRef = line[3]
            currAlt = line[4]
            if len(currRef) == 1 and len(currAlt) == 1:
                # SNP
                adjacentVars = [line]
                j = i + 1
                while j < len(vcfLines) and int(vcfLines[j][1]) - currPos == 1:
                    adjacentVars.append(vcfLines[j])
                    currPos = int(vcfLines[j][1])
                    j += 1
                if len(adjacentVars) >= minAdjacent:
                    correctedVars = examineAdjacent(adjacentVars, bamFile, ref2bit)
                    if len(correctedVars) > 0:
                        correctedVcfLines.extend(correctedVars)
                    else:
                        correctedVcfLines.extend(adjacentVars)
                else:
                    correctedVcfLines.extend(adjacentVars)
            else:
                # Not a SNP
                correctedVcfLines.append(line)

    writeVcf(correctedVcfLines, outputFile)


if __name__ == "__main__":
    main()

