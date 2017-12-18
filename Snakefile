from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

FTP = FTPRemoteProvider()


# load config file
configfile: 'config.yml'

###########################
# ALL rule
###########################

rule all:
    input:
        # positives
        expand("results/{db}/positives/design/final/{filter}/minSamples_{minsamples}.overlap_{overlap}.{transcriptDB}_minDistanceToTSS_{size}.probeLength.{length}.{ext}",
                db=['RoadMap'], filter=['enhancer'],
                minsamples=config['params']['positives']['minsamples'],
                overlap=config['params']['positives']['overlap'],
                transcriptDB=config['params']['transcript_DBs'], size=config['params']['minDistanceToTSS'],
                length=config['params']['probe_length'],
                ext=['fa', 'bed']),
        expand("results/{db}/negatives/design/final/{filter}/{group}.minsamplesWithinGroup_{minsamples}.maxOther_{other}.groupOverlap_{overlap}.{transcriptDB}_minDistanceToTSS_{size}.probeLength.{length}.{ext}",
                db=['RoadMap'], filter=['enhancer'],
                group=config['data']['RoadMap']['negative_groups'],
                overlap=config['params']['negatives']['overlap'],
                minsamples=config['params']['negatives']['minsamplesWithinGroup'],
                other=config['params']['negatives']['maxsamplesToOther'],
                transcriptDB=config['params']['transcript_DBs'], size=config['params']['minDistanceToTSS'],
                length=config['params']['probe_length'],
                ext=['fa', 'bed']),
###########################
# Download the DNase data
###########################

# imputed narrow peak data from RoadMap. Sample-IDs are defined in the config file
rule downloadRoadMapNarrowPeak:
    input:
    output:
        "data/RoadMap/imputedNarrowPeak/{id}.bed.{ext}"
    params:
        url="http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidatedImputed/narrowPeak/{id}-DNase.imputed.narrowPeak.bed.nPk.{ext}"
    shell:
        """
        wget {params.url} -O {output}
        """
# imputed signal data from RoadMap. Sample-IDs are defined in the config file
rule downloadRoadMapSignal:
    input:
    output:
        "data/RoadMap/imputedSignal/{id}.bigwig"
    params:
        url="http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/{id}-DNase.imputed.pval.signal.bigwig"
    shell:
        """
        wget {params.url} -O {output}
        """

####################
# positives
####################

include: "positives.rule"

####################
# negatives
####################

include: "negatives.rule"

#############################################
# select enhancers
# The regions have to have a certain distance apart from a TSS site
#############################################

def getRemote(wc):
    if (wc.database=="ensembl" and wc.release=="75"):
        return FTP.remote("ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz")
    elif (wc.database=="gencode" and wc.release == "27"):
        return FTP.remote("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.annotation.gtf.gz")
rule getTSS:
    input:
        getRemote
    output:
        "data/transcript_dbs/{database}/{database}.{release}.bed.gz"
    shell:
        """
        export LC_ALL=C;
        zcat {input} | awk -v OFS="\\t" '{{if ($3=="transcript") print $1,$2,$3,$4,$5,$7,$10}}' | \
        awk -v OFS="\\t" '{{if ($6=="+") print $1,$5-1,$5,substr($7,1,16),".",$6; if ($6=="-") print $1,$4-1,$4,substr($7,1,16),".",$6;}}' | \
        sed 's/"//g'  | sed 's/;//g' | sed 's/^chr//g' | egrep -v "^GL" | egrep -v "^HS" | egrep -v "^HG" | awk '{{print "chr"$0}}' | \
        sort -k 1,1 -k 2,2n | bgzip -c > {output};
        tabix {output};
        """

rule extendTSS:
    input:
        "data/transcript_dbs/{database}/{database}.{release}.bed.gz"
    output:
        out="data/transcript_dbs/{database}/extend/{database}.{release}.TSSextend_{size}.bed.gz",
        idx="data/transcript_dbs/{database}/extend/{database}.{release}.TSSextend_{size}.bed.gz.tbi"
    params:
        size="{size}"
    shell:
        """
        export LC_ALL=C;
        zcat {input} | awk -v OFS='\t' '{{ if ($6=="+") print $1,$2-1000*{params.size},$3,$4,$5,$6; else print $1,$2,$3+1000*{params.size},$4,$5,$6;}}' | \
        awk -v OFS='\t' '{{ if ($2<0) print $1,0,$3,$4,$5,$6; else print $0;}}' | sort -k 1,1 -k 2,2n | bgzip -c > {output.out};
        tabix {output.out};
        """

rule filterForEnhancers:
    input:
        filtered="results/{db}/{type}/filter/samples/{name}.bed.gz",
        idx="results/{db}/{type}/filter/samples/{name}.bed.gz.tbi",
        database="data/transcript_dbs/{database}/extend/{database}.{release}.TSSextend_{size}.bed.gz",
        db_idx="data/transcript_dbs/{database}/extend/{database}.{release}.TSSextend_{size}.bed.gz.tbi"
    output:
        filtered="results/{db}/{type}/filter/enhancer/{name}.{database}_{release}_minDistanceToTSS_{size}.bed.gz",
        idx="results/{db}/{type}/filter/enhancer/{name}.{database}_{release}_minDistanceToTSS_{size}.bed.gz.tbi"
    shell:
        """
        bedtools intersect -v -a {input.filtered} -b {input.database} | bgzip -c > {output.filtered};
        tabix {output.filtered};
        """




###########################
# get average signals
###########################

# to use the signal calls we have to convert bigwig files to wig and later to bed files
# bigwig to wig
rule bigwigToWig:
    input:
        bigwig="data/{db}/imputedSignal/{id}.bigwig"
    output:
        wig=temp("data/{db}/imputedSignal/{id}.wig")
    shell:
        """
        bigWigToWig {input.bigwig} {output.wig}
        """
# wig to bed
rule wigToBed:
    input:
        wig="data/{db}/imputedSignal/{id}.wig"
    output:
        bed="data/{db}/imputedSignal/{id}.bed.gz",
        idx="data/{db}/imputedSignal/{id}.bed.gz.tbi"
    shell:
        """
        export LC_ALL=C;
        wig2bed -d < {input.wig} | sort -k1,1 -k2,2n | bgzip -c >{output.bed};
        tabix {output.bed};
        """

def getSignalSamples(wc):
    return expand("data/{db}/imputedSignal/{id}.bed.gz",
                    db=wc.db, id=config['data'][wc.db]['samples'])
def getSignalSamplesIdx(wc):
    return expand("data/{db}/imputedSignal/{id}.bed.gz.tbi",
                    db=wc.db, id=config['data'][wc.db]['samples'])

# use the Sinal files of the samples to extract the center of the peak.
# use teh average of all centers as new center peak
rule extractAverageCenterOfPositiveRegions:
    input:
        regions="results/{db}/{type}/filter/{filter}/{name}.bed.gz",
        regions_idx="results/{db}/{type}/filter/{filter}/{name}.bed.gz.tbi",
        samples=getSignalSamples,
        samples_idx=getSignalSamplesIdx
    output:
        centers=temp("results/{db}/{type}/design/averageSignalCenter/{filter}/{name}.bed")
    run:
        import pysam
        with open(output.centers, 'w') as fout:
            with  pysam.TabixFile("%s" % input.regions) as regions:
                for region in regions.fetch():
                    line_split = region.split("\t")
                    chr=line_split[0]
                    start=int(line_split[1])+1
                    end=int(line_split[2])
                    counts = line_split[4:len(line_split)]
                    positions = []

                    for signal_file in input.samples:
                        maxPos = None
                        maxValue = -1
                        with  pysam.TabixFile("%s" % signal_file) as signals:
                             for signal in signals.fetch(chr,start,end):
                                 signal_split = signal.split("\t")
                                 value_start=int(signal_split[1])+1
                                 value_end=int(signal_split[2])
                                 value = float(signal_split[4])

                                 if (value > maxValue):
                                     maxValue=value
                                     maxPos = int(value_start + (value_end - value_start) / 2)
                        positions.append(maxPos)
                    positions =  [x for x in positions if x is not None]
                    fout.write("%s\t%d\t%d\t%s\t%d\n" % (chr, start-1, end, "\t".join(counts).strip(), int(sum(positions)/ len(positions))))

# convert it as bgzip and index it
rule bgzipExtractAverageCenterOfPositiveRegions:
    input:
        centers="results/{db}/{type}/design/{filter}/{name}.bed"
    output:
        centers="results/{db}/{type}/design/{filter}/{name}.bed.gz",
        idx="results/{db}/{type}/design/{filter}/{name}.bed.gz.tbi"
    shell:
        """
        cat {input.centers} | bgzip -c > {output.centers};
        tabix {output.centers}
        """

###########################
# create final probes
###########################


# using a bedfile with center (5th column) and the number of active counts (4th column)
# to create a fasta file using the probe_length
rule extractProbe:
    input:
        centers="results/{db}/{type}/design/averageSignalCenter/{filter}/{design}.bed.gz",
        idx="results/{db}/{type}/design/averageSignalCenter/{filter}/{design}.bed.gz.tbi",
        reference=config['files']['reference_genome']
    output:
        probes="results/{db}/{type}/design/final/{filter}/{design}.probeLength.{length}.fa",
        bed="results/{db}/{type}/design/final/{filter}/{design}.probeLength.{length}.bed"
    params:
        length="{length}"
    run:
        import pysam
        import math
        from pyfaidx import Fasta

        reference = Fasta(input['reference'])

        length=int(params.length)
        length=math.floor(length/2)

        with open(output.bed, 'w') as foutBed:
            with open(output.probes, 'w') as fout:
                with  pysam.TabixFile("%s" % input.centers) as regions:
                    for region in regions.fetch():
                        line_split = region.split("\t")
                        chr=line_split[0]
                        counts = line_split[3:len(line_split)-1]
                        # this is 1 based
                        center = int(line_split[-1])

                        # 1 based
                        start=center-length
                        # 1 based
                        end=center+length

                        sequence = reference[chr][start:end+1].seq

                        fout.write(">%s:%d-%d active_count:%s\n%s\n" % (chr, start, end, "_".join(counts), sequence))
                        foutBed.write("%s\t%d\t%d\t%s" % (chr, start-1,end, "\t".join(counts)))
