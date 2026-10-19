# SOXS UV-VIS Wavelength Calibration Analysis

This README documents the key results of adjusting and measuring the SOXS UVVIS wavelength calibration accuracy. 

* See here for the analysis and subsequent method used to [Tightening the SOXS UVVIS wavelength solution.md](Tightening%20the%20SOXS%20UVVIS%20wavelength%20solution.md).
* See here for the method used to [quantify the wavelength calibration accuracy](Measuring%20the%20UVVIS%20Wavelength%20Calibration%20Accuracy.md) (Geza's modified notebook). 

## Method to Measure Wavelength Calibration Accuracy

1. **Measure line shifts.** For each line, find the flux minimum and subtract the reference minimum. For **stellar** lines, use the star's rest-frame spectrum. **Telluric** lines use the telluric model. Each line is assigned to its order (u/g/r/i).
2. **Barycentric correction.** locations of stellar lines are corrected for Earth's velocity and rotation. Telluric lines form in Earth's atmosphere, so they stay a fixed reference.
3. **Fit star velocity.** In the r and i orders, take the median telluric offset per exposure, then compute `vStar` from each stellar line: shift minus telluric offset, converted to velocity, plus the barycentric correction. Take the 3-sigma MAD-clipped median.
4. **Compute residuals.** Subtract the expected star-velocity shift from each line's shift. Telluric lines expect 0, so their residual equals their shift. The residual is the wavelength-solution error.
5. **Per-order correction.** Compute clipped stats per order. The **u and g orders anchor on stellar lines** because they have no telluric lines. The **r and i orders anchor on telluric lines**, because their stellar lines fit the velocity and their residuals are about 0 by construction. The correction is minus the median shift.


Here are the before-and-after results of the skyline wavelength-calibration enhancement. The ‘before’ results come from data reduced with soxpipe v0.17.4, and the ‘after’ results from v0.19.0.

For each reduction, we reduce the same 11 standard stars. First, we reduce them with the daily arcs (normal operations mode), then with the attached arcs. 

## BEFORE (v0.17.4) with Daily Arcs 

<img width="1584" height="1084" alt="image" src="https://github.com/user-attachments/assets/b9698feb-6c71-4c1d-84d5-eed62a0956a9" />

| Line type   | Order   |   N |   Median [Å] |   σ [Å] |
|:------------|:--------|----:|-------------:|--------:|
| telluric    | r       |  11 |       +0.421 |   0.401 |
| telluric    | i       |  33 |       +0.202 |   0.816 |
| stellar     | u       |  22 |       +0.050 |   0.430 |
| stellar     | g       |  25 |       -0.146 |   0.596 |
| stellar     | r       |  14 |       +0.294 |   0.408 |
| stellar     | i       |   3 |       +0.186 |   1.117 |

**Telluric lines in i-band**

<img width="1548" height="441" alt="image" src="https://github.com/user-attachments/assets/e441858a-4547-4ebd-9aac-fb7ec627fc94" />

**H-beta stellar line in g-band**

<img width="1542" height="745" alt="image" src="https://github.com/user-attachments/assets/fab0897d-fd65-460b-968a-7d7559aad760" />


## B. Attached Arcs

<img width="1587" height="1084" alt="image" src="https://github.com/user-attachments/assets/a43466b2-7fb2-4acd-b04d-96119ec451d5" />

| Line type   | Order   |   N |   Median [Å] |   σ [Å] |
|:------------|:--------|----:|-------------:|--------:|
| telluric    | r       |  28 |       +0.097 |   0.846 |
| telluric    | i       |  84 |       -0.148 |   1.382 |
| stellar     | u       |  56 |       +0.373 |   0.593 |
| stellar     | g       |  59 |       +0.053 |   1.022 |
| stellar     | r       |  31 |       +0.007 |   0.757 |
| stellar     | i       |   3 |       +0.198 |   0.188 |

<img width="1544" height="445" alt="image" src="https://github.com/user-attachments/assets/a115e5be-abba-4040-a5a1-4a51bd7a2b3f" />

<img width="1544" height="446" alt="image" src="https://github.com/user-attachments/assets/6a0d01a0-2201-431d-b36b-760c5602c249" />


## C. Nightly Arcs with Skyline Calibration

<img width="1585" height="1085" alt="image" src="https://github.com/user-attachments/assets/5fa57203-2f9e-4972-9464-cfb51bc3be38" />

| Line type   | Order   |   N |   Median [Å] |   σ [Å] |
|:------------|:--------|----:|-------------:|--------:|
| telluric    | r       |  28 |       +0.069 |   0.956 |
| telluric    | i       |  84 |       -0.050 |   1.088 |
| stellar     | u       |  56 |       +0.193 |   0.570 |
| stellar     | g       |  59 |       +0.094 |   0.749 |
| stellar     | r       |  31 |       +0.001 |   1.148 |
| stellar     | i       |   3 |       +0.069 |   0.315 |

<img width="1547" height="443" alt="image" src="https://github.com/user-attachments/assets/bfd5d163-00c5-429b-9f7a-3b067c5ff78c" />

<img width="1543" height="443" alt="image" src="https://github.com/user-attachments/assets/9ba109d8-0b3b-49cb-bc56-c65b39eae8f7" />

## D. Attached Arcs with Skyline Calibration

<img width="1584" height="1084" alt="image" src="https://github.com/user-attachments/assets/eddb279e-2fce-41a4-9b14-36394300a70f" />

| Line type   | Order   |   N |   Median [Å] |   σ [Å] |
|:------------|:--------|----:|-------------:|--------:|
| telluric    | r       |  28 |       +0.066 |   0.794 |
| telluric    | i       |  84 |       -0.079 |   1.765 |
| stellar     | u       |  56 |       +0.319 |   0.449 |
| stellar     | g       |  59 |       -0.036 |   0.768 |
| stellar     | r       |  31 |       -0.131 |   1.180 |
| stellar     | i       |   3 |       +0.467 |   0.398 |

<img width="1545" height="443" alt="image" src="https://github.com/user-attachments/assets/6d64da41-5134-4943-9b2e-94f913e9df42" />


<img width="1545" height="446" alt="image" src="https://github.com/user-attachments/assets/dd196ed7-7723-40f1-92a9-a5d912583615" />

