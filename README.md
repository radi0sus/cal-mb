# Cal-MB

A Python 3 script for easy calibration of ⁵⁷Fe-Mößbauer (MB) spectra from a ⁵⁷Fe sample. 

## External modules

`lmfit`
`numpy` 
`scipy` 
`matplotlib`

## Quick start

Start the script with:
```console
python3 cal-mb.py 57Fe_calib_raw_data.ws5
```
calculates the folding point `FP`, `v0` (channel where the velocity is zero),
`vmax` (maximum velocity), and the velocity / channel `f`.

The file should contain intensities from a multi-channel analyzer. The popular WissEl format
has the extension `.ws5`. In principle any raw data (not only WissEl) can be processed.

Terminal output:
```
-------------------------------------
Results for 57Fe_calib_raw_data.ws5 :
File modified on 17.12.2023 19:53:55
-------------------------------------
FP (channel) = 256.8019±0.0754
v₀ (channel) = 125.7526±0.0858
vmax /mm·s⁻¹ =  -4.7010±0.0151
f /mm·s⁻¹/c  =   0.0369±0.0001
-------------------------------------
```
You can use the first three parameters to fit MB spectra with [fit-mb](https://github.com/radi0sus/fit-mb). 

Start the script with:
```console
python3 cal-mb.py 57Fe_calib_raw_data.ws5 -s
```
calculates the same values as described above. In addition a `matplotlib` window is shown that 
summarizes the fit results.

<img src='examples\Figure_1.png' alt='Fit' width=600 align='center'> 

The results will not be saved. To keep the output you have to start the script with:
```console
python3 cal-mb.py 57Fe_calib_raw_data.ws5 > calib_from_today.txt
```
To keep the figure, you have to click the floppy symbol (similar to 💾) in the `matplotlib` window.
