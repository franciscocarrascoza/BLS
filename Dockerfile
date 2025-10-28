FROM ubuntu:22.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    pkg-config \
    libxdrfile-dev \
    libtng-io-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/bls
COPY . .

RUN cmake -S . -B build -DBLS_USE_XDRFILE=ON -DBLS_USE_TNG=ON \
    && cmake --build build -j$(nproc) \
    && ctest --test-dir build

CMD ["./build/bls_analyze", "--help"]

