FROM python:3.10.0-slim  

# Install system dependencies, including wget, gcc, and Julia
RUN apt-get update && apt-get install -y --no-install-recommends gcc wget netcat vim git

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

# Install microbetag Python dependencies
COPY ./microbetag/requirements/requirements.txt /tmp/mtg_requirements.txt
RUN pip wheel --no-cache-dir --no-deps --wheel-dir /wheels -r /tmp/mtg_requirements.txt
RUN pip install --no-cache /wheels/*

# Copy project files to the appropriate directory
COPY . $APP_HOME

# Change to the app user -- uncomment in production
WORKDIR $APP_HOME

RUN cd microbetag &&  pip install .