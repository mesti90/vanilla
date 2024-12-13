# Dockerfile: Lightweight Version (Python Slim Base Image)
# Multi-stage build using python:3.11-slim

# Stage 1: Build stage (dependencies installation)
FROM python:3.11-slim AS builder

# Install system dependencies and third-party tools
RUN apt-get update && apt-get install -y --no-install-recommends \
	gcc build-essential git wget unzip

#Install Snippy
FROM builder AS snippy
	RUN git clone https://github.com/tseemann/snippy.git

#Install SnpEff
FROM builder AS snpeff
	RUN VERSION=v4_3t && \
	wget -O snpEff.zip https://sourceforge.net/projects/snpeff/files/snpEff_"$VERSION"_core.zip/download && \
	unzip snpEff.zip && \
	rm snpEff.zip 

#Final Stage
FROM builder AS final
	
	#Install additional dependencies
	RUN apt-get update && apt-get install -y --no-install-recommends \
	emboss python3-biopython python3-yaml openjdk-17-jdk-headless && rm -rf /var/lib/apt/lists/*
	
	#Copy binaries and scripts
	COPY --from=snippy /snippy/bin/* /usr/local/bin
	COPY --from=snippy /snippy/perl5/Snippy/Version.pm /usr/share/perl5/Snippy/Version.pm
	COPY --from=snpeff /snpEff/ /usr/local/snpEff/
	
	#Set environment variables
	ENV PATH="/usr/bin:/usr/local/bin/:/usr/local/sbin:/sbin:/bin:/usr/sbin:/usr/local/snpEff/:/usr/local/snpEff/scripts" \
	LANG=C.UTF-8 \
	LC_ALL=C.UTF-8
