FROM openjdk:16-slim-buster

ARG PICARD_VERSION=2.23.8
ARG SAMTOOLS_VERSION=1.11
ARG HTSLIB_VERSION=1.11
ARG STAR_VERSION=2.7.6a
ARG BEDTOOLS_VERSION=2.29.2
ARG DROPSEQ_TOOLS_VERSION=2.4.0

WORKDIR /working

RUN apt-get update \
    && apt-get install -y \
		wget \
		make \
		g++ \
		bzip2 \
		libbz2-dev \
        zlib1g-dev \
        libncurses5-dev \
        liblzma-dev \
		libcurl4-openssl-dev \
		libssl-dev \
	&& mkdir -p /working/software/bin/ \
	&& cd /working/software/ \
	# pull Picard jar
	&& wget "https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar" \
	# install samtools
    && wget "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
    && mkdir /working/software/samtools-install \
    && tar -xjf "samtools-${SAMTOOLS_VERSION}.tar.bz2" \
        -C /working/software/samtools-install \
        --strip-components 1 \
    && cd /working/software/samtools-install \
    && ./configure --prefix=/working/software/samtools \
    && make \
    && make install \
    && ln -s /working/software/samtools/bin/* /usr/local/bin/ \
    && cd /working/software \
    && rm -r "samtools-${SAMTOOLS_VERSION}.tar.bz2" samtools-install/ \
	# install htslib
	&& wget "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" \
    && mkdir /working/software/htslib-install \
    && tar -xjf "htslib-${HTSLIB_VERSION}.tar.bz2" \
        -C /working/software/htslib-install \
        --strip-components 1 \
    && cd /working/software/htslib-install \
    && ./configure --prefix=/working/software/htslib \
    && make \
    && make install \
    && ln -s /working/software/htslib/bin/* /usr/local/bin/ \
    && cd /working/software \
    && rm -r "htslib-${HTSLIB_VERSION}.tar.bz2" htslib-install/ \
    # install STAR
    && wget "https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" \
	&& mkdir /working/software/star-install \
	&& tar -xzf "${STAR_VERSION}.tar.gz" \
	&& cd "STAR-${STAR_VERSION}/source" \
	&& make STAR \
	&& cp STAR /working/software/bin/ \
	&& cd /working/software \ 
	&& rm -rf "${STAR_VERSION}.tar.gz" "STAR-${STAR_VERSION}" \
    # install bedtools
    && wget -O /working/software/bedtools "https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools.static.binary" \
    && chmod a+x /working/software/bedtools \
    && ln -s /working/software/bedtools /working/software/bin/ \
    # install Drop-seq tools
    && wget "https://github.com/broadinstitute/Drop-seq/releases/download/v${DROPSEQ_TOOLS_VERSION}/Drop-seq_tools-${DROPSEQ_TOOLS_VERSION}.zip" \
    && unzip "Drop-seq_tools-${DROPSEQ_TOOLS_VERSION}.zip" \
    && ln -sf /working/software/Drop-seq_tools-${DROPSEQ_TOOLS_VERSION}/* /working/software/bin/ \
    && rm "Drop-seq_tools-${DROPSEQ_TOOLS_VERSION}.zip"

ENV PATH="/working/software/bin:$PATH"

CMD ["/bin/bash"]
