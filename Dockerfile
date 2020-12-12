FROM jupyter/scipy-notebook:c1b0cf6bf4d6
MAINTAINER Matthew Ware "matt.ware@raytheon.com"

# for binder environment
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

USER root
#RUN adduser --disabled-password \
#    --gecos "Default user" \
#   --uid ${NB_UID} \
#    ${NB_USER}

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

USER jovyan
RUN pip install --no-cache-dir notebook==5.*

# PyCall needs this
RUN conda install pyqt matplotlib

# 'right' version of pygsti
RUN git clone https://github.com/pyGSTio/pyGSTi.git --branch v0.9.3
RUN cd pyGSTi; pip install -e .

# install julia
WORKDIR /home/jovyan
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.3-linux-x86_64.tar.gz
RUN tar xf julia-0.6.3-linux-x86_64.tar.gz
RUN mv julia-d55cadc350/ julia-0.6.3/


RUN /home/jovyan/julia-0.6.3/bin/julia -e 'Pkg.add("Convex.jl",v"0.5.0");Pkg.add("Formatting");Pkg.add("IJulia");Pkg.add("PyPlot");Pkg.add("StatsBase");Pkg.add("NPZ");Pkg.add("Glob");Pkg.add("SchattenNorms");Pkg.add("Distributions");Pkg.add("Cliffords");'

RUN /home/jovyan/julia-0.6.3/bin/julia -e 'Pkg.clone("QuantumInfo");Pkg.checkout("QuantumInfo", "fix/0.6transforms");'

RUN export LD_LIBRARY_PATH=/opt/conda/lib:$LD_LIBRARY_PATH
RUN /home/jovyan/julia-0.6.3/bin/julia -e 'ENV["PYTHON"]="/opt/conda/bin/python";Pkg.build("PyCall");Pkg.build("PyPlot")'

# copy in the data
COPY --chown=jovyan:users data/ work/data/
COPY --chown=jovyan:users scripts/ work/scripts/
COPY --chown=jovyan:users notebooks/ work/notebooks/
COPY --chown=jovyan:users things-to-copy/ work/things-to-copy/

RUN chown -R jovyan work/
RUN chmod -R 777 work/
