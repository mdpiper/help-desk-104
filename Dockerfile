FROM condaforge/miniforge3

WORKDIR /opt
COPY requirements.txt .
COPY run.py .
COPY build_bmi_fortran_heat ./build_bmi_fortran_heat

RUN conda update -n base -c conda-forge conda
RUN conda install -c conda-forge --file requirements.txt

RUN git clone https://github.com/csdms/babelizer
WORKDIR /opt/babelizer
RUN git checkout mdpiper/help-desk-104
RUN conda install -c conda-forge --file requirements.txt
RUN make install

# Persist CONDA_PREFIX through build. Needed to build pymt_heat.
ENV CONDA_PREFIX=/opt/conda

WORKDIR /opt/build_bmi_fortran_heat/bmi-example-fortran/_build
RUN cmake ../ -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} &&  make all install

WORKDIR /opt/build_bmi_fortran_heat
RUN babelize init babel_heat.toml
WORKDIR /opt/build_bmi_fortran_heat/pymt_heat
RUN pip install -e .

WORKDIR /opt
CMD [ "python", "run.py" ]

# docker build --tag help-desk-104 .
# docker run --rm help-desk-104
