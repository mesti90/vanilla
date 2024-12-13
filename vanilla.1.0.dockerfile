# Dockerfile: Lightweight Version (Python Slim Base Image)
# Multi-stage build using python:3.11-slim

# Stage 1: Build stage (dependencies installation)
FROM python:3.11-slim AS builder

# Install system dependencies and third-party tools
RUN apt-get update && apt-get install -y --no-install-recommends \
gcc build-essential git wget unzip libstdc++6 libbz2-1.0 autoconf \
bzip2 gzip perl less libdatetime-perl libxml-simple-perl libdigest-md5-perl \
zlib1g-dev liblzma-dev libbz2-dev xz-utils curl libcurl4-gnutls-dev libssl-dev \
default-jdk-headless parallel libvcflib-dev libvcflib-tools

#Install Snippy
FROM builder AS snippy
	RUN git clone https://github.com/tseemann/snippy.git && sed -i '298a print $cfg "codon.Bacterial_and_Plant_Plastid: TTT/F, TTC/F, TTA/L, TTG/L+, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, TGT/C, TGC/C, TGA/*, TGG/W, CTT/L, CTC/L, CTA/L, CTG/L+, CCT/P, CCC/P, CCA/P, CCG/P, CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I+, ATC/I+, ATA/I+, ATG/M+, ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/R, AGG/R, GTT/V, GTC/V, GTA/V, GTG/V+, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, GGT/G, GGC/G, GGA/G, GGG/G\\n";' snippy/bin/snippy

#Install SnpEff
FROM builder AS snpeff
	RUN VERSION=v4_3t && \
	wget -O snpEff.zip https://sourceforge.net/projects/snpeff/files/snpEff_"$VERSION"_core.zip/download && \
	unzip snpEff.zip && \
	rm snpEff.zip 

FROM builder AS bwa
	RUN git clone https://github.com/lh3/bwa.git && cd bwa && make -j 8

FROM builder AS minimap
	RUN git clone https://github.com/lh3/minimap2 && cd minimap2 && make -j 8

FROM builder AS wgsim
	RUN git clone https://github.com/lh3/wgsim && cd wgsim && gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm

FROM builder AS samclip
	RUN wget https://raw.githubusercontent.com/tseemann/samclip/master/samclip -O samclip && chmod +x samclip

FROM builder AS any2fasta
	RUN wget https://raw.githubusercontent.com/tseemann/any2fasta/master/any2fasta -O any2fasta && chmod +x any2fasta

FROM builder AS vt
	RUN git clone https://github.com/atks/vt.git && cd vt && git submodule update --init --recursive && make -j 8

#Final Stage
FROM builder AS final
	
	#Install additional dependencies
	
	RUN apt-get -y --no-install-recommends install bison flex bioperl hmmer emboss python3-biopython python3-yaml freebayes samtools bamtools bedtools bcftools any2fasta samclip snp-sites seqtk mc && apt-get clean && apt-get autoclean && rm -rf /var/lib/apt/lists/* 
	
	#Copy binaries and scripts
	
	WORKDIR /usr/local/bin
	
	COPY --from=snippy /snippy/bin/* ./
	COPY --from=snippy /snippy/perl5/Snippy/Version.pm /usr/share/perl5/Snippy/Version.pm
	COPY --from=snpeff /snpEff/ /usr/local/snpEff/
	COPY --from=bwa /bwa/bwa ./
	COPY --from=minimap /minimap2/minimap2 ./
	COPY --from=wgsim /wgsim/wgsim ./
	COPY --from=wgsim /wgsim/wgsim_eval.pl ./
	COPY --from=vt /vt/vt ./
	
	#Set environment variables
	ENV PATH="/usr/bin:/usr/local/bin/:/usr/local/sbin:/sbin:/bin:/usr/sbin:/usr/local/snpEff/:/usr/local/snpEff/scripts:/usr/share/freebayes/scripts" \
	LANG=C.UTF-8 \
	LC_ALL=C.UTF-8
	
	ENTRYPOINT ["/bin/bash"]
