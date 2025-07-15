# Multi-stage build for GBRS (Genome Reconstruction from RNA-Seq)
# Stage 1: Build stage for compiling dependencies
FROM ubuntu:22.04 AS builder

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Install system dependencies for building
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    bzip2 \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    git \
    unzip \
    python3 \
    python3-dev \
    python3-venv \
    python3-pip \
    libncurses5-dev \
    libncursesw5-dev \
    libssl-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libpng-dev \
    libfreetype6-dev \
    libjpeg-dev \
    libtiff-dev \
    libhdf5-dev \
    libblas-dev \
    liblapack-dev \
    gfortran \
    software-properties-common \
    && rm -rf /var/lib/apt/lists/*

# Install Python 3.12 from deadsnakes PPA (available for both AMD64 and ARM64)
RUN add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y python3.12 python3.12-dev python3.12-venv && \
    rm -rf /var/lib/apt/lists/*

# Create Python virtual environment
RUN python3.12 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Upgrade pip and install wheel
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# Copy source code and install Python dependencies
COPY . /tmp/gbrs
WORKDIR /tmp/gbrs
RUN pip install --no-cache-dir . && \
    rm -rf /tmp/gbrs

# Install SAMtools
RUN wget -q https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
    && tar -xjf samtools-1.18.tar.bz2 \
    && cd samtools-1.18 \
    && ./configure --prefix=/opt/samtools \
    && make \
    && make install \
    && cd .. \
    && rm -rf samtools-1.18*

# Install Bowtie (platform-specific)
ARG TARGETPLATFORM
RUN if [ "$TARGETPLATFORM" = "linux/amd64" ]; then \
        wget -q -O bowtie.zip https://github.com/BenLangmead/bowtie/releases/download/v1.3.1/bowtie-1.3.1-linux-x86_64.zip; \
    elif [ "$TARGETPLATFORM" = "linux/arm64" ]; then \
        wget -q -O bowtie.zip https://github.com/BenLangmead/bowtie/releases/download/v1.3.1/bowtie-1.3.1-linux-aarch64.zip; \
    else \
        echo "Unsupported platform: $TARGETPLATFORM"; exit 1; \
    fi \
    && unzip bowtie.zip -d /opt/ \
    && ln -s /opt/bowtie-1.3.1-linux-* /opt/bowtie \
    && rm bowtie.zip

# Stage 2: Runtime stage
FROM ubuntu:22.04 AS runtime

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV PATH="/opt/venv/bin:/opt/samtools/bin:/opt/bowtie:$PATH"

# Install runtime dependencies only
RUN apt-get update && apt-get install -y \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    libz1 \
    libbz2-1.0 \
    liblzma5 \
    libncurses5 \
    libncursesw5 \
    libssl3 \
    zlib1g \
    libcurl4 \
    libxml2 \
    libpng16-16 \
    libfreetype6 \
    libjpeg-turbo8 \
    libtiff5 \
    libhdf5-103-1 \
    libblas3 \
    liblapack3 \
    libgfortran5 \
    software-properties-common \
    && rm -rf /var/lib/apt/lists/*

# IMPORTANT: Install Python 3.12 in the runtime image so the venv works
# In a multi-stage build, the venv created in the builder points to the Python binary path in the builder image.
# If the runtime image does not have the same Python version at the same path, the venv will be broken ("bad interpreter").
RUN add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y python3.12 python3.12-venv && \
    rm -rf /var/lib/apt/lists/*

# Copy Python virtual environment from builder
COPY --from=builder /opt/venv /opt/venv

# Copy SAMtools and Bowtie from builder
COPY --from=builder /opt/samtools /opt/samtools
COPY --from=builder /opt/bowtie /opt/bowtie

# Create non-root user
RUN groupadd -r gbrs && useradd -r -g gbrs -m -s /bin/bash gbrs

# Fix permissions on virtual environment
RUN chown -R gbrs:gbrs /opt/venv

# Create application directory
RUN mkdir -p /app /data /output && \
    chown -R gbrs:gbrs /app /data /output

# Copy application code
COPY --chown=gbrs:gbrs . /app/gbrs

# Switch to non-root user
USER gbrs

# Set working directory
WORKDIR /app/gbrs

# Install the package (regular install, not editable)
RUN pip install --no-cache-dir . && \
    rm -rf /app/gbrs

# Create data and output directories
RUN mkdir -p /data /output

# Add health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import gbrs; print('GBRS is healthy')" || exit 1

# Default command
CMD ["gbrs", "--help"]

# Metadata
LABEL maintainer="Matthew Vincent <matt.vincent@jax.org>, Mike Lloyd <mike.lloyd@jax.org>"
LABEL description="GBRS: Genome Reconstruction from RNA-Seq"
LABEL version="1.1.0"
LABEL org.opencontainers.image.source="https://github.com/churchill-lab/gbrs"
LABEL org.opencontainers.image.vendor="Churchill Lab"
LABEL org.opencontainers.image.licenses="MIT"
