[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4323548.svg)](https://doi.org/10.5281/zenodo.4323548)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/matthewware/RGST-analysis/HEAD)

# RGST-analysis environment

Code for analyzing and plotting GST data from [arxiv:1803.01818](https://arxiv.org/abs/1803.01818)

## Docker
The `Docker/` folder contains a `Dockerfile` and instructions for building the analysis environment used in the tomographic reconstruction, analysis, and plotting of data in the `data/` folder. It is highly recommended that this image be used to run the analysis -- all the code was originally written in 2016-2017 and depends on rather outdated version of Python and Julia, as well as the assorted libraries used, and the image downloads all the appropriate versions. Building the image will take ~10 minutes depending on your machine and internet connection.

Note that the pygsti reconstruction can be very resource intensive! It can take ~2hrs when run serially on an 8 core machine (this corresponds to performing 14 different pygsti reconstructions). You can skip this step and go straight to plotting/playing with the operators using the supplied `*.npz` files containing the LGST results. There are more details in the analysis notebooks.

If you _would_ like to
do the full analysis, we recommend you make sure your Docker is configured to
use a significant amount of RAM (i.e. as much as you can offer. The default value on Mac OSx is 2 GB for example). Data files in
the `data/` directory were reconstructed using roughly 6 GB of RAM.

Once the docker image is built the `start_env.sh`, `stop_env.sh`, and
`copy_data.sh` utility scripts can be used to manage the environment.

### mybinder.org
The top-level Dockerfile is for the [mybinder.org](https://mybinder.org/) build. This can be launched from the badge above. The mybinder build will take some time. The amount of time could depend on how much is cached on the mybinder.org side. The compute resources on the mybinder.org side are modest so we recommend you not try to do the `pyGSTi` reconstruction or the calculation of the Hessian-based error bars on mybinder.org. These two are particularly resource intensive. The results of these calculations are precomputed and the results saved in `.npz` files. See the analysis notebook for more details.


## Analysis notebooks
  * The `work/Randomize.ipynb` notebook walks through the process of generating randomized GST experiments based on sequences generated by pyGSTi.
  * The `work/RGST-analysis.ipynb` contains all the code for analyzing the experimental data (after it has been collated into the format necessary for pyGSTi to run), as well as generating plots and error bars.

## Disclaimers

The code and scripts found here are provided as-is for the purposes of reproducing the results described in [arxiv:1803.01818](https://arxiv.org/abs/1803.01818).

## Acknowledgements

This work was funded by LPS/ARO grant W911NF-14-C-0048.  The content of this paper does not necessarily reflect the position or the policy of the Government,
and no official endorsement should be inferred.

## Citations

Please see the manuscript by Ware et al., [arxiv:1803.01818](https://arxiv.org/abs/1803.01818) (2018) for additional details,
and cite the arXiv version and the most recently published version of the manuscript when referring to this work.

```bibtex
@misc{ware2018,
      title={Experimental demonstration of Pauli-frame randomization on a superconducting qubit},
      author={Matthew Ware and Guilhem Ribeill and Diego Rist{\`e} and Colm A. Ryan and Blake Johnson and Marcus P. da Silva},
      year={2018},
      eprint={1803.01818},
      archivePrefix={arXiv},
      primaryClass={quant-ph}
}
```

## Contributors

Experiment generation scripts by [Marcus Silva](https://marcusps.github.io/). Data analysis by 
[Matthew Ware](https://matthewware.dev/) and Marcus Silva.
