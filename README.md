# Visium stitcher

The repository houses a Visium stitcher, allowing for a simple joining of adjacent tissue slices into a joint AnnData object for analysis. This is a two-step process, with the first step requiring aligning the images in Fiji. The resulting `.xml` file is subsequently used to transform the images and coordinates in Python.

Please refer to the [example notebook](notebooks/demo.ipynb) for exact Fiji instructions and a subsequent Python demonstration. The demo data is available at `ftp://ftp.sanger.ac.uk/pub/users/kp9/demo_data.tar.gz`.

## Installation

```bash
pip install git+https://github.com/Teichlab/visium_stitcher.git
```

## Citation

Please cite our [preprint](https://www.biorxiv.org/content/10.1101/2022.04.27.489800v1.abstract)
```
@article{zhang2022human,
  title={A human embryonic limb cell atlas resolved in space and time},
  author={Zhang, Bao and He, Peng and Lawrence, John E and Wang, Shuaiyu and Tuck, Elizabeth and Williams, Brian A and Roberts, Kenny and Kleshchevnikov, Vitalii and Mamanova, Lira and Bolt, Liam and others},
  journal={bioRxiv},
  pages={2022--04},
  year={2022},
  publisher={Cold Spring Harbor Laboratory}
}
```