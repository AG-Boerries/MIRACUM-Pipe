FROM debian:9.9

ADD . /opt/MIRACUM-Pipe

RUN apt-get update && apt-get upgrade -y

RUN /opt/MIRACUM-Pipe/tools/setup.debian.sh

# install tools
RUN /opt/MIRACUM-Pipe/tools/install.sh

# set PATH
RUN export MY_PATH=`cat /opt/MIRACUM-Pipe/tools/my_path`
ENV PATH="$MY_PATH:$PATH"

# set LD_LIBRARY_PATH
RUN export MY_LD_LIBRARY_PATH=`cat /opt/MIRACUM-Pipe/tools/my_ld_library_path`
ENV LD_LIBRARY_PATH="$MY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH"

RUN echo $PATH