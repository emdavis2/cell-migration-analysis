BIOS 611 Project 
================
## Quick Summary:
The goal of this project is to analyze shape and motion data from movies of cells migrating and find correlations that can explain how the two are related and what unique behaviors are exhibited on different substrates the cells are moving on. 

## Full Description:
Cell migration is characterized by an elongation with a leading edge of the cell membrane that contains protrusions. These protrusions adhere to the substrate the cell is situated on, probe the physical properties of the substrate, and generate forces that ultimately cause the rear of the cell to retract, thus leading to overall cell motion. The lifetime and relative positioning of these focal adhesions affects the persistence of the cellular motion. Cell migration has typically been modeled as a persistent random walk, but this assumption fails to explain many aspects of the data, such as the non-Gaussian distribution of velocities and the non-singular exponential decay of the velocity autocorrelation function. Incorporating cell shape into motion behavior has yet to be done and may explain the deviations from a persistent random walk that we observe.  Looking at this data set, I hope to learn how whole cell cytoskeleton dynamics play a role in the persistence and other observable migration behaviors of cells through the correlation of cell shape and motion. Additionally, I hope to characterize the differences cells display while migrating on gel, a more physiologically relevant substrate, and glass. Gaining insights from the analysis of this data, I aim to eventually make a mathematical model that incorporates shape and motion and can describe the observed differences on gel and glass. 

## How to Run:
Clone this repositorty and navigate to the folder of this repository on your computer (where the Docker file is located) and build the docker container by running:

```
docker build -t python_latex .
```

Now, run this container with the following: 

```
docker run -v $(pwd):/home/ -it python_latex bash
```

Now your terminal should be interactive with the docker container we have just built. We can use the Makefile to make the figures. To generate the autocorrelation figures for the cells migrating on glass, use:

```
make figures/acf_figures_glass
```

To do the same for gel, run:

```
make figures/acf_figures_gel
```

To build the report, run:

```
make writeup.pdf
```