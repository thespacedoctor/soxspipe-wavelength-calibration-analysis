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


## BEFORE (v0.17.4) with Attached Arcs 

<img width="1585" height="1081" alt="image" src="https://github.com/user-attachments/assets/8e8f1d6c-539a-4e7f-9a08-2ae858d00536" />


| Line type   | Order   |   N |   Median [Å] |   σ [Å] |
|:------------|:--------|----:|-------------:|--------:|
| telluric    | r       |  10 |       -0.045 |   0.453 |
| telluric    | i       |  30 |       +0.042 |   1.127 |
| stellar     | u       |  20 |       +0.166 |   0.317 |
| stellar     | g       |  23 |       -0.074 |   0.362 |
| stellar     | r       |  13 |       -0.271 |   0.479 |
| stellar     | i       |   3 |       +0.129 |   0.267 |

**Telluric lines in i-band**

<img width="1534" height="749" alt="image" src="https://github.com/user-attachments/assets/dfab0436-c115-47a2-8158-4d2e7b6563ae" />

**H-beta stellar line in g-band**

<img width="1543" height="447" alt="image" src="https://github.com/user-attachments/assets/2da77ddd-4ecb-4816-8b88-eaed6987f70f" />



## AFTER (v0.19.0) with Daily Arcs 

<img width="1585" height="1084" alt="image" src="https://github.com/user-attachments/assets/cd03e0af-5c1c-4465-abc9-8097d6b5ee80" />


| Line type   | Order   |   N |   Median [Å] |   σ [Å] |
|:------------|:--------|----:|-------------:|--------:|
| telluric    | r       |  11 |       +0.035 |   0.221 |
| telluric    | i       |  30 |       +0.026 |   0.220 |
| stellar     | u       |  20 |       +0.293 |   0.274 |
| stellar     | g       |  24 |       +0.223 |   0.278 |
| stellar     | r       |  13 |       +0.026 |   0.273 |
| stellar     | i       |   3 |       +0.121 |   0.253 |

**Telluric lines in i-band**

<img width="1546" height="446" alt="image" src="https://github.com/user-attachments/assets/2986a5d4-fc8c-402a-be21-3a373085ac71" />

**H-beta stellar line in g-band**

<img width="1541" height="742" alt="image" src="https://github.com/user-attachments/assets/11b422be-dbff-4d55-8ba0-9f3c6aa85802" />


## AFTER (v0.19.0) with Attached Arcs 

<img width="1582" height="1080" alt="image" src="https://github.com/user-attachments/assets/e7bbc269-4b75-47cd-9615-dd58ce7f09e2" />

| Line type   | Order   |   N |   Median [Å] |   σ [Å] |
|:------------|:--------|----:|-------------:|--------:|
| telluric    | r       |  11 |       +0.014 |   0.244 |
| telluric    | i       |  31 |       +0.013 |   0.953 |
| stellar     | u       |  19 |       +0.231 |   0.350 |
| stellar     | g       |  23 |       +0.144 |   0.395 |
| stellar     | r       |  12 |       -0.015 |   0.329 |
| stellar     | i       |   3 |       +0.119 |   0.344 |

**Telluric lines in i-band**

<img width="1546" height="442" alt="image" src="https://github.com/user-attachments/assets/9767ef9b-e2c5-4e17-aa43-31bfa319e2ea" />

**H-beta stellar line in g-band**

<img width="1542" height="744" alt="image" src="https://github.com/user-attachments/assets/2e3e194d-0350-44e3-a8cb-56e6359cba22" />


