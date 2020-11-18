FROM jupyter/scipy-notebook:c1b0cf6bf4d6
MAINTAINER Matthew Ware "matt.ware@raytheon.com"

#RUN /bin/bash -c "source activate python2"
USER jovyan

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


RUN /home/jovyan/julia-0.6.3/bin/julia -e 'Pkg.add("Convex.jl",v"0.5.0");Pkg.add("Formatting");Pkg.add("IJulia");Pkg.add("PyPlot");Pkg.add("QuantumInfo");Pkg.add("StatsBase");Pkg.add("NPZ");Pkg.add("Glob");Pkg.add("SchattenNorms");Pkg.add("Distributions");Pkg.add("Cliffords");'

RUN export LD_LIBRARY_PATH=/opt/conda/lib:$LD_LIBRARY_PATH
RUN /home/jovyan/julia-0.6.3/bin/julia -e 'ENV["PYTHON"]="/opt/conda/bin/python";Pkg.build("PyCall");Pkg.build("PyPlot")'

# copy in the data
COPY --chown=jovyan:users data/ work/data/
COPY --chown=jovyan:users scripts/ work/scripts/
COPY --chown=jovyan:users notebooks/ work/notebooks/
COPY --chown=jovyan:users things-to-copy/ work/things-to-copy/

RUN chown -R jovyan work/
RUN chmod -R 777 work/
