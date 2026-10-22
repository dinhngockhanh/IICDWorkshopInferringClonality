# Workshop image: R 4.5.1, PyClone-VI 0.2.0, and PhyloWGS on Python 2.7.18.
# uv.lock and renv.lock pin the Python 3 and R libraries.
# PhyloWGS cannot live in either lock; its commit and Python 2 packages are ARGs below.
FROM rocker/r-ver:4.5.1

ARG PHYLOWGS_COMMIT=2205be16e9ad68426e837229b593f47e36b3547e
ARG OPENSSL_VERSION=1.1.1w
ARG PYTHON27_VERSION=2.7.18
ARG NUMPY_VERSION=1.16.6
ARG SCIPY_VERSION=1.2.3
ARG ETE2_VERSION=2.3.10

COPY docker/apt-packages.txt /tmp/apt-packages.txt
RUN apt-get update \
    && xargs -a /tmp/apt-packages.txt apt-get install -y --no-install-recommends \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:0.10.7 /uv /uvx /bin/
ENV UV_PYTHON_INSTALL_DIR=/opt/uv/python \
    UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1
WORKDIR /opt/workshop
COPY pyproject.toml uv.lock /opt/workshop/
RUN uv python install 3.12 \
    && uv sync --frozen --no-dev
ENV PATH="/opt/workshop/.venv/bin:${PATH}"

RUN R -q -e 'install.packages("https://cloud.r-project.org/src/contrib/renv_1.3.1.tar.gz", repos = NULL)'
COPY renv.lock /opt/workshop/renv.lock
COPY docker/restore.R /opt/workshop/restore.R
RUN Rscript /opt/workshop/restore.R

RUN wget -O /tmp/openssl.tar.gz "https://www.openssl.org/source/old/1.1.1/openssl-${OPENSSL_VERSION}.tar.gz" \
    && tar -xzf /tmp/openssl.tar.gz -C /tmp \
    && cd "/tmp/openssl-${OPENSSL_VERSION}" \
    && ./config --prefix=/opt/openssl11 --openssldir=/opt/openssl11 shared -fPIC \
    && make -j"$(nproc)" \
    && make install_sw \
    && rm -rf /tmp/openssl*

COPY docker/enable_py27_ssl.py /tmp/enable_py27_ssl.py
RUN wget -O /tmp/python.tgz "https://www.python.org/ftp/python/${PYTHON27_VERSION}/Python-${PYTHON27_VERSION}.tgz" \
    && tar -xzf /tmp/python.tgz -C /tmp \
    && cd "/tmp/Python-${PYTHON27_VERSION}" \
    && python3 /tmp/enable_py27_ssl.py \
    && CFLAGS="-O2 -fcommon -std=gnu99 -Wno-error=implicit-function-declaration" \
       LDFLAGS="-Wl,-rpath,/opt/openssl11/lib" \
       ./configure --prefix=/opt/python27 --enable-unicode=ucs4 --with-ensurepip=no \
    && make -j"$(nproc)" \
    && make install \
    && rm -rf "/tmp/Python-${PYTHON27_VERSION}" /tmp/python.tgz
ENV LD_LIBRARY_PATH=/opt/openssl11/lib \
    PATH="/opt/python27/bin:/opt/workshop/.venv/bin:${PATH}" \
    PYTHON2=/opt/python27/bin/python

RUN /opt/python27/bin/python -c "import ssl,sys; print(ssl.OPENSSL_VERSION); print(sys.version)"
RUN wget -O /tmp/get-pip.py https://bootstrap.pypa.io/pip/2.7/get-pip.py \
    && /opt/python27/bin/python /tmp/get-pip.py
RUN CFLAGS="-fcommon -std=gnu99 -Wno-error=implicit-function-declaration" \
    CXXFLAGS="-fcommon" \
    /opt/python27/bin/pip install --no-cache-dir "numpy==${NUMPY_VERSION}"
# gfortran 10+ rejects the rank mismatches in SciPy 1.2.3's ARPACK sources.
RUN CFLAGS="-fcommon -std=gnu99 -Wno-error=implicit-function-declaration" \
    CXXFLAGS="-fcommon" \
    FFLAGS="-fcommon -fallow-argument-mismatch" \
    F77FLAGS="-fcommon -fallow-argument-mismatch" \
    F90FLAGS="-fcommon -fallow-argument-mismatch" \
    NPY_DISTUTILS_APPEND_FLAGS=1 \
    /opt/python27/bin/pip install --no-cache-dir "scipy==${SCIPY_VERSION}"
RUN /opt/python27/bin/pip install --no-cache-dir "ete2==${ETE2_VERSION}" \
    && /opt/python27/bin/python -c "import ssl, numpy, scipy, ete2; print(ssl.OPENSSL_VERSION); print(numpy.__version__); print(scipy.__version__)"

RUN git clone https://github.com/morrislab/phylowgs.git /opt/phylowgs \
    && git -C /opt/phylowgs checkout "${PHYLOWGS_COMMIT}" \
    && cd /opt/phylowgs \
    && g++ -o mh.o -O3 mh.cpp util.cpp $(gsl-config --cflags --libs) \
    && test -s /opt/phylowgs/mh.o

ENV PHYLOWGS_DIR=/opt/phylowgs \
    OMP_NUM_THREADS=1 \
    RGL_USE_NULL=TRUE \
    MPLBACKEND=Agg \
    PATH="/opt/workshop/.venv/bin:/opt/python27/bin:${PATH}"

WORKDIR /work
CMD ["bash", "src/run_all.sh"]
