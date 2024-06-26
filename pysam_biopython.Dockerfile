FROM quay.io/biocontainers/pysam:0.22.1--py312hcfdcdd7_1

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    libbz2-dev \
    liblzma-dev \
    make \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir biopython
