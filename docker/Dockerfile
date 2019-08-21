FROM miracumpipe_tools:latest

COPY ./*.sh /opt/MIRACUM-Pipe/
COPY ./conf/default.yaml /opt/MIRACUM-Pipe/conf

# ENTRYPOINT ["/opt/MIRACUM-Pipe/interface.sh"]

# define volumes
VOLUME /opt/MIRACUM-Pipe/assets/input
VOLUME /opt/MIRACUM-Pipe/assets/output
VOLUME /opt/MIRACUM-Pipe/assets/references

VOLUME /opt/MIRACUM-Pipe/tools/annovar
VOLUME /opt/MIRACUM-Pipe/tools/gatk
