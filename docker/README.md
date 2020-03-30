# Docker

For development, we recommend to use the hierarchical docker stack, in order to save time while testing and building (the full production build can take about 3 hours)

This stack is based on

- Dockerfile_base
    Based on `debian:10` the linux environment is included into this file. All system tools are setup.
- Dockerfile_R
  All R packages are installed on top of `Dockerfile_base` which are required.
- Dockerfile_tools
  All tools which are used by our pipeline are downloaded and installed on top of `Dockerfile_R`
- Dockerfile_report
  This is intended if only scripts intended to change the report are changed. Tis is build upon `Dockerfile_tools`.
- Dockerfile_interface
  Finally only the shell interface data are changed, based on `Dockerfile_report`.
- Dockerfile
  This file is used for production build and is a composition of all other docker files and bases on `debian:10`.

## Use the dockerized version

For more information of who to use the docker container, see [MIRACUM-Pipe-docker](https://github.com/AG-Boerries/MIRACUM-Pipe-docker).

## Building behind a proxy server

If you are behind a proxy server and need to build the docker from scratch you need to add the following build arguments in the docker/build file.

- --build-arg http_proxy=http://proxy.server.de:port
- --build-arg https_proxy=https://proxy.server.de:port

Additionally, the proxy settings can also be added in the Dockerfiles as ENV parameters

- ENV http_proxy http://proxy.server.de:port
- ENV https_proxy https://proxy.server.de:port
