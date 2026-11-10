# Noah-MP v5.2.1 for HydroTuring.
#
# Upstream's source, from the HRLDAS revision that carries v5.2.1 as its submodule,
# is compiled with HRLDAS's own makefiles and flags (with -O2 added), unmodified.
# ht_noahmp_driver.F90 links against the objects and only feeds them and reads them.
# hrldas.exe and the point-data converter stay in the image so that verify_driver.py
# can compare the driver with upstream's own entry point. Nothing is fetched at run
# time. The UCAR licences are copied beside the adapter; the image is never pushed.
#
# Both stages use the same bookworm-slim index, pinned by digest. The apt packages
# resolve from the bookworm archive at build time; README.md lists their versions.

FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171 AS build
ARG HRLDAS_COMMIT=cd96df470220f7d7133cdbccd5f9c5355cf173e2
# The tag v5.2.1, and the submodule pointer of HRLDAS_COMMIT.
ARG NOAHMP_COMMIT=17751dcd7a2442a3c57138f530f523daf08dc1c5
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates git make cpp gfortran libnetcdff-dev libnetcdf-dev \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src/hrldas && cd /src/hrldas \
 && git remote add origin https://github.com/NCAR/hrldas.git \
 && git fetch -q --depth 1 origin "${HRLDAS_COMMIT}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${HRLDAS_COMMIT}" \
 && git ls-tree HEAD noahmp | grep -q "${NOAHMP_COMMIT}"
RUN git init -q /src/hrldas/noahmp && cd /src/hrldas/noahmp \
 && git remote add origin https://github.com/NCAR/noahmp.git \
 && git fetch -q --depth 1 origin "${NOAHMP_COMMIT}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${NOAHMP_COMMIT}"
# HRLDAS's serial gfortran options, with the paths of this image and without
# Jasper, which only its GRIB converter needs.
RUN cd /src/hrldas/hrldas \
 && sed -e 's#^ F90FLAGS .*# F90FLAGS = -O2 -g -fconvert=big-endian -fbounds-check -fno-range-check -fallow-argument-mismatch#' \
        -e 's#^ NETCDFMOD .*# NETCDFMOD = -I/usr/include#' \
        -e 's#^ NETCDFLIB .*# NETCDFLIB = -lnetcdff -lnetcdf#' \
        -e 's#^ LIBJASPER .*# LIBJASPER =#' \
        -e 's#^ INCJASPER .*# INCJASPER =#' \
        arch/user_build_options.gfortran.serial > user_build_options
# The official targets in Makefile order, less the GRIB converter.
RUN cd /src/hrldas/hrldas/Utility_routines && make -j1
RUN cd /src/hrldas/noahmp/utility && make -j1
RUN cd /src/hrldas/noahmp/src && make -j1
RUN cd /src/hrldas/noahmp/drivers/hrldas && make -j1
RUN cd /src/hrldas/urban/wrf && make -j1
RUN cd /src/hrldas/hrldas/IO_code && make -j1
RUN cd /src/hrldas/hrldas/run && make -j1
COPY ht_noahmp_driver.F90 /src/ht_noahmp_driver.F90
RUN cd /src \
 && gfortran -O2 -g -ffree-line-length-none -fno-range-check \
      -Ihrldas/noahmp/src -Ihrldas/noahmp/utility -Ihrldas/noahmp/drivers/hrldas -Ihrldas/hrldas/Utility_routines \
      -c ht_noahmp_driver.F90 -o ht_noahmp_driver.o \
 && gfortran -o ht_noahmp_driver ht_noahmp_driver.o \
      hrldas/noahmp/drivers/hrldas/*.o $(ls hrldas/noahmp/src/*.o | grep -v NoahmpIOVarType) \
      hrldas/noahmp/utility/*.o hrldas/hrldas/Utility_routines/*.o \
      -lnetcdff -lnetcdf
RUN mkdir -p /src/point \
 && gfortran -I/usr/include -o /src/point/create_point_data.exe \
      /src/hrldas/hrldas/HRLDAS_forcing/run/examples/single_point/create_point_data.f90 -lnetcdff -lnetcdf

FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171
ENV PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1
RUN apt-get update && apt-get install -y --no-install-recommends \
      libgfortran5 libnetcdff7 python3 python3-numpy python3-netcdf4 \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/ht_noahmp_driver /opt/noahmp/ht_noahmp_driver
COPY --from=build /src/hrldas/hrldas/run/hrldas.exe /opt/noahmp/hrldas.exe
COPY --from=build /src/point/create_point_data.exe /opt/noahmp/create_point_data.exe
COPY --from=build /src/hrldas/hrldas/HRLDAS_forcing/run/examples/single_point/bondville.dat /opt/noahmp/bondville.dat
COPY --from=build /src/hrldas/noahmp/parameters/NoahmpTable.TBL /opt/noahmp/NoahmpTable.TBL
COPY --from=build /src/hrldas/hrldas/run/namelist.hrldas_example /opt/noahmp/namelist.hrldas_example
COPY --from=build /src/hrldas/noahmp/LICENSE.txt /model/LICENSE.noahmp.txt
COPY --from=build /src/hrldas/LICENSE.txt /model/LICENSE.hrldas.txt
RUN chmod -R a+rX /opt/noahmp /model

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py
COPY verify_driver.py /model/verify_driver.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
