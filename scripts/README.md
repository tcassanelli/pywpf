
# Scripts

This directory contains standalone code not part of the main package `pywpf`.

# Examples

To execute example data please download the files [here](https://drive.google.com/drive/folders/1SOJtwLRdiWdfYEAD3PQeZvM3XB_Mk48C), or access the Aqueye and Iqueye website [https://web.oapd.inaf.it/zampieri/aqueye-iqueye/](https://web.oapd.inaf.it/zampieri/aqueye-iqueye/) then click data and pca_data.

**Warning**: These files are large (1 GBs or more) and the scripts will take over one hour to be completed in a laptop computer.

To execute the software first transform the text files to numpy binaries:

```bash
python3 data2npy.py -i pathto/20091215-031945UTC_crab-NTT.baricentrizzato.tempo2.txt
python3 b0531+21.py -i pathto/20091215-031945UTC_crab-NTT.baricentrizzato.tempo2.npy
```

Then the data outputs will be a figure and the computed data sets from `pywpca`, in the `pywpf_out/` directory (`cd pathto/`).

```bash
├── pywpf_out
    ├── 20091215-031945UTC_crab-NTT.baricentrizzato.tempo2_run.pdf
    ├── 20091215-031945UTC_crab-NTT.baricentrizzato.tempo2-000
        ├── M200.npy
        ├── info.dat
```

If you have problems please submit an issue. Thanks!
