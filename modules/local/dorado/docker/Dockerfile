FROM nvidia/cuda:12.3.2-cudnn9-runtime-rockylinux9

RUN dnf install -y \
    python3 \
    python3-pip \
    python3-devel \
    gcc \
    wget \
    procps \
    && dnf clean all

RUN pip3 install --upgrade pip

# Specify the version of the ONT Dorado software to install
ENV DORADO_VERSION=1.3.0

WORKDIR /opt

RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-linux-x64.tar.gz \
    && tar -xzf dorado-${DORADO_VERSION}-linux-x64.tar.gz \
    && rm dorado-${DORADO_VERSION}-linux-x64.tar.gz \
    && mv dorado-${DORADO_VERSION}-linux-x64 dorado


ENV PATH="/opt/dorado/bin:${PATH}"

CMD ["dorado"]
