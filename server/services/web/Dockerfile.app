FROM python:3.10-slim

# Install system dependencies, including wget, gcc, and Julia
RUN apt-get update && apt-get install -y --no-install-recommends gcc wget netcat-openbsd vim git

# Install Python dependencies
RUN pip install --upgrade pip flake8==6.0.0

# Set environment variables and directories
ENV HOME=/home/app
ENV APP_HOME=/home/app/web
RUN mkdir -p $APP_HOME

# Install Julia
WORKDIR /opt
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.1-linux-x86_64.tar.gz \
    && tar -zxvf julia-1.7.1-linux-x86_64.tar.gz \
    && rm julia-1.7.1-linux-x86_64.tar.gz \
    && echo "export PATH=/opt/julia-1.7.1/bin:$PATH" >> /etc/bash.bashrc

# Set PATH for Julia
ENV PATH="/opt/julia-1.7.1/bin:${PATH}"

# Newer Linux kernels refuse to load libraries that request an executable stack, as Julia's do
RUN apt-get update && apt-get install -y --no-install-recommends patchelf && rm -rf /var/lib/apt/lists/* \
 && find /opt/julia-1.7.1/lib -name "*.so*" -type f -exec patchelf --clear-execstack {} \;

# Install FlashWeave with Julia
RUN julia -e 'using Pkg; Pkg.add("FlashWeave"); Pkg.add("PyCall")'

# Install additional Python packages (manta, fuzzywuzzy)
RUN pip install fuzzywuzzy
RUN pip install git+https://github.com/hariszaf/manta.git@scipy-version

# To enable pip insatll . for microbetag
RUN apt-get update && apt-get install -y \
    build-essential \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Install app Python dependencies
COPY ./requirements.txt /tmp/app_requirements.txt
RUN pip wheel --no-cache-dir --no-deps --wheel-dir /wheels -r /tmp/app_requirements.txt
RUN pip install --no-cache /wheels/*

# microbetag, at the commit the server ran with as a git submodule (v1.0.4).
# Until the server moves to the package of this repository, it is installed from GitHub.
ARG MICROBETAG_REF=49e94a77bc944d563544440fc9ac1823e8ab204f
RUN git clone https://github.com/msysbio/microbetag.git /opt/microbetag \
 && git -C /opt/microbetag checkout ${MICROBETAG_REF}

# Install microbetag Python dependencies
RUN pip wheel --no-cache-dir --no-deps --wheel-dir /wheels -r /opt/microbetag/requirements/requirements.txt
RUN pip install --no-cache /wheels/*

# Copy project files to the appropriate directory
COPY . $APP_HOME

# Change to the app user -- uncomment in production
WORKDIR $APP_HOME

# v1.0.4's setup.py runs `pip install pyshorteners` from inside its own install, which current pip
# does not allow; pyshorteners is already installed above, so drop that custom install command.
RUN cd /opt/microbetag \
 && sed -i "s/'install': CustomInstallCommand,//" setup.py \
 && pip install .