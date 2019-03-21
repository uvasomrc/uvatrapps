# `uvatrapps` Docker instructions

This repository contains Docker instructions for a Shiny server that runs [UVATRAPPS apps](https://github.com/uvasomrc/uvatrapps).

A built version of the container is hosted on [Docker Hub](https://hub.docker.com/r/somrc/consultr). 

To run the container locally:

```
docker pull somrc/uvatrapps
```

```
docker run -d -p 80:80 somrc/uvatrapps
```

Alternatively, you can build the container locally based on the `Dockerfile` in this repository:

```
git clone https://github.com/uvasomrc/uvatrapps.git
```

```
cd uvatrapps/docker
```

```
docker build -t uvatrapps .
```

```
docker run -d -p 80:80 uvatrapps
```
