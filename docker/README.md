# Docker + raycloudtools tutorial
---
#### Building and using raycloudtools with docker
#### author: Tiago de Conto
#### date: June 24, 2023
---

The provided dockerfile allows easy building. It is based on [Ubuntu22.04](https://hub.docker.com/_/ubuntu/).

Be sure to have the latest [Docker](https://docs.docker.com/engine/install/) installed in your system. To build the raycloudtools image, run:

```
docker build -f Dockerfile -t raycloudtools .
```

Or on Windows:

```
docker build -f Dockerfile.txt -t raycloudtools .
```

This image will download and compile raycloudtools and all its dependencies (including treetools, so you can process rayextract outputs). You can then run the `raycloudtools` image as a standalone executable. The basic command is `docker run -v your/local/datadir:/data raycloudtools` followed by the specific tool you want to use, as in: 

```
docker run -v your/local/datadir:/data raycloudtools rayimport
```

Here's an example of a pipeline to process a forest plot:

```
# create a local directory where to store your data and the raycloudtools outputs
mkdir mydata
cd mydata

# [optional] create a local environment variable with the arguments you 
# want to pass to docker - so you don't have to rewrite it everytime
CMD="run --rm -u $(id -u $USER) -v $PWD:/data raycloudtools"

# generate some sample data
docker $CMD raycreate forest 3
```

The command above will generate a `forest.ply` file in your local directory. From there, you can use `rayextract` to get other products:

```
docker $CMD rayextract terrain forest.ply

docker $CMD rayextract trunks forest.ply

docker $CMD rayextract forest forest.ply

docker $CMD rayextract trees forest.ply forest_mesh.ply

docker $CMD rayextract leaves forest.ply forest_trees.txt
```

The outputs from the above commands will be at your local directory. You can visually inspect the `ply` files using software such as [CloudCompare](https://www.danielgm.net/cc/).
